/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2020 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 *
 * Authors: Wolfgang Bangerth, 1999,
 *          Guido Kanschat, 2011
 *          Luca Heltai, 2021
 */
#include "poisson.h"

#include <deal.II/base/multithread_info.h>

#include <deal.II/grid/grid_refinement.h>

#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/trilinos_precondition.h>

#include <deal.II/numerics/error_estimator.h>

using namespace dealii;

template <int dim>
Poisson<dim>::Poisson()
  : dof_handler(triangulation)
  , solver_control("Solver control", 1000, 1e-12, 1e-12)
{
  add_parameter("Finite element degree", fe_degree);
  add_parameter("Mapping degree", mapping_degree);
  add_parameter("Number of global refinements", n_refinements);
  add_parameter("Output filename", output_filename);
  add_parameter("Forcing term expression", forcing_term_expression);
  add_parameter("Dirichlet boundary condition expression",
                dirichlet_boundary_conditions_expression);
  add_parameter("Coefficient expression", coefficient_expression);
  add_parameter("Exact solution expression", exact_solution_expression);
  add_parameter("Neumann boundary condition expression",
                neumann_boundary_conditions_expression);

  add_parameter("Local pre-refinement grid size expression",
                pre_refinement_expression);

  add_parameter("Dirichlet boundary ids", dirichlet_ids);
  add_parameter("Neumann boundary ids", neumann_ids);

  add_parameter("Problem constants", constants);
  add_parameter("Grid generator function", grid_generator_function);
  add_parameter("Grid generator arguments", grid_generator_arguments);
  add_parameter("Number of refinement cycles", n_refinement_cycles);

  add_parameter("Estimator type",
                estimator_type,
                "",
                this->prm,
                Patterns::Selection("exact|kelly|residual"));

  add_parameter("Marking strategy",
                marking_strategy,
                "",
                this->prm,
                Patterns::Selection("global|fixed_fraction|fixed_number"));

  add_parameter("Coarsening and refinement factors",
                coarsening_and_refinement_factors);

  add_parameter("Use direct solver", use_direct_solver);

  this->prm.enter_subsection("Error table");
  error_table.add_parameters(this->prm);
  this->prm.leave_subsection();
}


template <int dim>
void
Poisson<dim>::initialize(const std::string &filename)
{
  ParameterAcceptor::initialize(filename,
                                "last_used_parameters.prm",
                                ParameterHandler::Short);
}



template <int dim>
void
Poisson<dim>::parse_string(const std::string &parameters)
{
  ParameterAcceptor::prm.parse_input_from_string(parameters);
  ParameterAcceptor::parse_all_parameters();
}



template <int dim>
void
Poisson<dim>::make_grid()
{
  const auto vars = dim == 1 ? "x" : dim == 2 ? "x,y" : "x,y,z";
  pre_refinement.initialize(vars, pre_refinement_expression, constants);
  GridGenerator::generate_from_name_and_arguments(triangulation,
                                                  grid_generator_function,
                                                  grid_generator_arguments);

  for (unsigned int i = 0; i < n_refinements; ++i)
    {
      for (const auto &cell : triangulation.active_cell_iterators())
        if (pre_refinement.value(cell->center()) < cell->diameter())
          cell->set_refine_flag();
      triangulation.execute_coarsening_and_refinement();
    }

  std::cout << "Number of active cells: " << triangulation.n_active_cells()
            << std::endl;
}



template <int dim>
void
Poisson<dim>::refine_grid()
{
  // Cells have been marked in the mark() method.
  triangulation.execute_coarsening_and_refinement();
}



template <int dim>
void
Poisson<dim>::setup_system()
{
  if (!fe)
    {
      fe              = std::make_unique<FE_Q<dim>>(fe_degree);
      mapping         = std::make_unique<MappingQGeneric<dim>>(mapping_degree);
      const auto vars = dim == 1 ? "x" : dim == 2 ? "x,y" : "x,y,z";
      forcing_term.initialize(vars, forcing_term_expression, constants);
      coefficient.initialize(vars, coefficient_expression, constants);
      exact_solution.initialize(vars, exact_solution_expression, constants);

      dirichlet_boundary_condition.initialize(
        vars, dirichlet_boundary_conditions_expression, constants);

      neumann_boundary_condition.initialize(
        vars, neumann_boundary_conditions_expression, constants);
    }

  dof_handler.distribute_dofs(*fe);
  std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;
  constraints.clear();
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);

  for (const auto &id : dirichlet_ids)
    VectorTools::interpolate_boundary_values(
      *mapping, dof_handler, id, dirichlet_boundary_condition, constraints);
  constraints.close();


  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
  sparsity_pattern.copy_from(dsp);
  system_matrix.reinit(sparsity_pattern);
  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
  error_per_cell.reinit(triangulation.n_active_cells());
}



template <int dim>
void
Poisson<dim>::assemble_system()
{
  QGauss<dim>     quadrature_formula(fe->degree + 1);
  QGauss<dim - 1> face_quadrature_formula(fe->degree + 1);

  FEValues<dim> fe_values(*mapping,
                          *fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  FEFaceValues<dim> fe_face_values(*mapping,
                                   *fe,
                                   face_quadrature_formula,
                                   update_values | update_quadrature_points |
                                     update_JxW_values);

  const unsigned int dofs_per_cell = fe->n_dofs_per_cell();
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
      cell_matrix = 0;
      cell_rhs    = 0;
      for (const unsigned int q_index : fe_values.quadrature_point_indices())
        {
          for (const unsigned int i : fe_values.dof_indices())
            for (const unsigned int j : fe_values.dof_indices())
              cell_matrix(i, j) +=
                (coefficient.value(
                   fe_values.quadrature_point(q_index)) * // a(x_q)
                 fe_values.shape_grad(i, q_index) *       // grad phi_i(x_q)
                 fe_values.shape_grad(j, q_index) *       // grad phi_j(x_q)
                 fe_values.JxW(q_index));                 // dx
          for (const unsigned int i : fe_values.dof_indices())
            cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q)
                            forcing_term.value(
                              fe_values.quadrature_point(q_index)) * // f(x_q)
                            fe_values.JxW(q_index));                 // dx
        }

      if (cell->at_boundary())
        //  for(const auto face: cell->face_indices())
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
          if (neumann_ids.find(cell->face(f)->boundary_id()) !=
              neumann_ids.end())
            {
              fe_face_values.reinit(cell, f);
              for (const unsigned int q_index :
                   fe_face_values.quadrature_point_indices())
                for (const unsigned int i : fe_face_values.dof_indices())
                  cell_rhs(i) += fe_face_values.shape_value(i, q_index) *
                                 neumann_boundary_condition.value(
                                   fe_face_values.quadrature_point(q_index)) *
                                 fe_face_values.JxW(q_index);
            }

      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(
        cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
    }
}



template <int dim>
void
Poisson<dim>::solve()
{
  if (use_direct_solver == true)
    {
      SparseDirectUMFPACK system_matrix_inverse;
      system_matrix_inverse.initialize(system_matrix);
      system_matrix_inverse.vmult(solution, system_rhs);
    }
  else
    {
      SolverCG<Vector<double>> solver(solver_control);
#ifdef DEAL_II_WITH_TRILINOS
      TrilinosWrappers::PreconditionAMG amg;
      amg.initialize(system_matrix);
      solver.solve(system_matrix, solution, system_rhs, amg);
#else
      solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
#endif
    }
  constraints.distribute(solution);
}



template <int dim>
void
Poisson<dim>::estimate()
{
  if (estimator_type == "exact")
    {
      error_per_cell = 0;
      QGauss<dim> quad(fe->degree + 1);
      VectorTools::integrate_difference(*mapping,
                                        dof_handler,
                                        solution,
                                        exact_solution,
                                        error_per_cell,
                                        quad,
                                        VectorTools::H1_seminorm);
    }
  else if (estimator_type == "kelly")
    {
      std::map<types::boundary_id, const Function<dim> *> neumann;
      for (const auto id : neumann_ids)
        neumann[id] = &neumann_boundary_condition;

      QGauss<dim - 1> face_quad(fe->degree + 1);
      KellyErrorEstimator<dim>::estimate(*mapping,
                                         dof_handler,
                                         face_quad,
                                         neumann,
                                         solution,
                                         error_per_cell,
                                         ComponentMask(),
                                         &coefficient);
    }
  else if (estimator_type == "residual")
    {
      // h_T || f+\Delta u_h ||_0,T
      // + \sum over faces
      // 1/2 (h_F)^{1/2} || [n.\nabla u_h] ||_0,F

      QGauss<dim - 1> face_quad(fe->degree + 1);
      QGauss<dim>     quad(fe->degree + 1);

      std::map<types::boundary_id, const Function<dim> *> neumann;
      for (const auto id : neumann_ids)
        neumann[id] = &neumann_boundary_condition;

      KellyErrorEstimator<dim>::estimate(*mapping,
                                         dof_handler,
                                         face_quad,
                                         neumann,
                                         solution,
                                         error_per_cell,
                                         ComponentMask(),
                                         &coefficient);

      FEValues<dim> fe_values(*mapping,
                              *fe,
                              quad,
                              update_hessians | update_JxW_values |
                                update_quadrature_points);

      std::vector<double> local_laplacians(quad.size());


      double residual_L2_norm = 0;

      unsigned int cell_index = 0;
      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          fe_values.reinit(cell);

          fe_values.get_function_laplacians(solution, local_laplacians);
          residual_L2_norm = 0;
          for (const auto q_index : fe_values.quadrature_point_indices())
            {
              const auto arg =
                (local_laplacians[q_index] +
                 forcing_term.value(fe_values.quadrature_point(q_index)));
              residual_L2_norm += arg * arg * fe_values.JxW(q_index);
            }
          error_per_cell[cell_index] +=
            cell->diameter() * std::sqrt(residual_L2_norm);

          ++cell_index;
        }
    }
  else
    {
      AssertThrow(false, ExcNotImplemented());
    }
  auto global_estimator = error_per_cell.l2_norm();
  error_table.add_extra_column("estimator", [global_estimator]() {
    return global_estimator;
  });
  error_table.error_from_exact(*mapping, dof_handler, solution, exact_solution);
}



template <int dim>
void
Poisson<dim>::mark()
{
  if (marking_strategy == "global")
    {
      for (const auto &cell : triangulation.active_cell_iterators())
        cell->set_refine_flag();
    }
  else if (marking_strategy == "fixed_fraction")
    {
      GridRefinement::refine_and_coarsen_fixed_fraction(
        triangulation,
        error_per_cell,
        coarsening_and_refinement_factors.second,
        coarsening_and_refinement_factors.first);
    }
  else if (marking_strategy == "fixed_number")
    {
      GridRefinement::refine_and_coarsen_fixed_number(
        triangulation,
        error_per_cell,
        coarsening_and_refinement_factors.second,
        coarsening_and_refinement_factors.first);
    }
  else
    {
      Assert(false, ExcInternalError());
    }
}

template <int dim>
void
Poisson<dim>::output_results(const unsigned cycle) const
{
  DataOut<dim>          data_out;
  DataOutBase::VtkFlags flags;
  flags.write_higher_order_cells = true;
  data_out.set_flags(flags);
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");
  auto interpolated_exact = solution;
  VectorTools::interpolate(*mapping,
                           dof_handler,
                           exact_solution,
                           interpolated_exact);
  data_out.add_data_vector(interpolated_exact, "exact");
  data_out.add_data_vector(error_per_cell, "estimator");
  data_out.build_patches(*mapping,
                         std::max(mapping_degree, fe_degree),
                         DataOut<dim>::curved_inner_cells);
  std::string   fname = output_filename + "_" + std::to_string(cycle) + ".vtu";
  std::ofstream output(fname);
  data_out.write_vtu(output);
}



template <int dim>
void
Poisson<dim>::print_system_info()
{
  std::cout << "Number of cores  : " << MultithreadInfo::n_cores() << std::endl
            << "Number of threads: " << MultithreadInfo::n_threads()
            << std::endl;
}



template <int dim>
void
Poisson<dim>::run()
{
  print_system_info();
  make_grid();
  for (unsigned int cycle = 0; cycle < n_refinement_cycles; ++cycle)
    {
      setup_system();
      assemble_system();
      solve();
      estimate();
      output_results(cycle);
      if (cycle < n_refinement_cycles - 1)
        {
          mark();
          refine_grid();
        }
    }
  error_table.output_table(std::cout);
}

template class Poisson<1>;
template class Poisson<2>;
template class Poisson<3>;