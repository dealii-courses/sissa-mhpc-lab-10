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
#include "base_block_problem.h"

#include <deal.II/base/work_stream.h>

#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_tools.h>

using namespace dealii;

template <int dim>
BaseBlockProblem<dim>::BaseBlockProblem(
  const std::vector<std::string> component_names,
  const std::string &            problem_name)
  : BaseProblem<dim>(component_names.size(), problem_name)
  , component_names(component_names)
{}



template <int dim>
void
BaseBlockProblem<dim>::setup_system()
{
  TimerOutput::Scope timer_section(this->timer, "setup_system");
  if (!this->fe)
    {
      this->fe = FETools::get_fe_by_name<dim>(this->fe_name);
      this->mapping =
        std::make_unique<MappingQGeneric<dim>>(this->mapping_degree);
      const auto vars = dim == 1 ? "x" : dim == 2 ? "x,y" : "x,y,z";
      this->forcing_term.initialize(vars,
                                    this->forcing_term_expression,
                                    this->constants);
      this->exact_solution.initialize(vars,
                                      this->exact_solution_expression,
                                      this->constants);

      this->dirichlet_boundary_condition.initialize(
        vars, this->dirichlet_boundary_conditions_expression, this->constants);

      this->neumann_boundary_condition.initialize(
        vars, this->neumann_boundary_conditions_expression, this->constants);
    }

  this->dof_handler.distribute_dofs(*this->fe);

  // renumber dofs in a blockwise manner.
  std::vector<unsigned int> blocks(this->n_components);
  unsigned int              i = 0;
  Assert(component_names.size() > 0, ExcInternalError());
  AssertDimension(this->n_components, component_names.size());
  blocks[0] = i;
  for (unsigned int j = 1; j < this->n_components; ++j)
    {
      if (component_names[j] == component_names[j - 1])
        blocks[j] = i;
      else
        blocks[j] = ++i;
    }
  DoFRenumbering::component_wise(this->dof_handler, blocks);

  dofs_per_block = DoFTools::count_dofs_per_fe_block(this->dof_handler, blocks);

  locally_owned_dofs =
    this->dof_handler.locally_owned_dofs().split_by_block(dofs_per_block);


  IndexSet non_blocked_locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(this->dof_handler,
                                          non_blocked_locally_relevant_dofs);
  locally_relevant_dofs =
    non_blocked_locally_relevant_dofs.split_by_block(dofs_per_block);

  this->pcout << "Number of degrees of freedom: " << this->dof_handler.n_dofs()
              << std::endl;


  this->constraints.clear();
  this->constraints.reinit(non_blocked_locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(this->dof_handler, this->constraints);

  for (const auto &id : this->dirichlet_ids)
    VectorTools::interpolate_boundary_values(*this->mapping,
                                             this->dof_handler,
                                             id,
                                             this->dirichlet_boundary_condition,
                                             this->constraints);
  this->constraints.close();


  TrilinosWrappers::BlockSparsityPattern dsp(locally_owned_dofs,
                                             locally_owned_dofs,
                                             locally_relevant_dofs,
                                             this->mpi_communicator);

  DoFTools::make_sparsity_pattern(this->dof_handler,
                                  dsp,
                                  this->constraints,
                                  false);
  // SparsityTools::distribute_sparsity_pattern(dsp,
  //                                            locally_owned_dofs,
  //                                            mpi_communicator,
  //                                            locally_relevant_dofs);


  system_block_matrix.reinit(dsp);

  block_solution.reinit(locally_owned_dofs, this->mpi_communicator);
  system_block_rhs.reinit(locally_owned_dofs, this->mpi_communicator);

  locally_relevant_block_solution.reinit(locally_owned_dofs,
                                         locally_relevant_dofs,
                                         this->mpi_communicator);

  this->error_per_cell.reinit(this->triangulation.n_active_cells());

  // Now call anything that may be needed
  this->setup_system_call_back();
}



template <int dim>
void
BaseBlockProblem<dim>::assemble_system_one_cell(
  const typename DoFHandler<dim>::active_cell_iterator &,
  ScratchData &,
  CopyData &)
{
  Assert(false, ExcPureFunctionCalled());
}



template <int dim>
void
BaseBlockProblem<dim>::copy_one_cell(const CopyData &copy)
{
  this->constraints.distribute_local_to_global(copy.matrices[0],
                                               copy.vectors[0],
                                               copy.local_dof_indices[0],
                                               system_block_matrix,
                                               system_block_rhs);
}



template <int dim>
void
BaseBlockProblem<dim>::assemble_system()
{
  TimerOutput::Scope timer_section(this->timer, "assemble_system");
  QGauss<dim>        quadrature_formula(this->fe->degree + 1);
  QGauss<dim - 1>    face_quadrature_formula(this->fe->degree + 1);

  ScratchData scratch(*this->mapping,
                      *this->fe,
                      quadrature_formula,
                      update_values | update_gradients |
                        update_quadrature_points | update_JxW_values,
                      face_quadrature_formula,
                      update_values | update_quadrature_points |
                        update_JxW_values);

  CopyData copy(this->fe->n_dofs_per_cell());

  auto worker = [&](const auto &cell, auto &scratch, auto &copy) {
    assemble_system_one_cell(cell, scratch, copy);
  };

  auto copier = [&](const auto &copy) { copy_one_cell(copy); };

  using CellFilter =
    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>;

  WorkStream::run(CellFilter(IteratorFilters::LocallyOwnedCell(),
                             this->dof_handler.begin_active()),
                  CellFilter(IteratorFilters::LocallyOwnedCell(),
                             this->dof_handler.end()),
                  worker,
                  copier,
                  scratch,
                  copy);

  system_block_matrix.compress(VectorOperation::add);
  system_block_rhs.compress(VectorOperation::add);
}



template <int dim>
void
BaseBlockProblem<dim>::solve()
{
  TimerOutput::Scope timer_section(this->timer, "solve");
  // SolverCG<LA::MPI::BlockVector> solver(solver_control);
  // LA::MPI::PreconditionAMG  amg;
  // amg.initialize(system_matrix);
  // solver.solve(system_matrix, solution, system_rhs, amg);
  // constraints.distribute(solution);
  // locally_relevant_solution = solution;
}



template class BaseBlockProblem<1>;
template class BaseBlockProblem<2>;
template class BaseBlockProblem<3>;