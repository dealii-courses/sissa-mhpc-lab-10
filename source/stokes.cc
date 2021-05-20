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
#include "stokes.h"

#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/linear_operator_tools.h>
#include <deal.II/lac/solver_gmres.h>

using namespace dealii;

namespace
{
  std::vector<std::string>
  get_component_names(int dim)
  {
    std::vector<std::string> names(dim + 1, "u");
    names[dim] = "p";
    return names;
  }
} // namespace

template <int dim>
Stokes<dim>::Stokes()
  : BaseBlockProblem<dim>(get_component_names(dim),
                          "Stokes<" + std::to_string(dim) + ">")
  , velocity(0)
  , pressure(dim)
{
  // Output the vector result.
  this->add_data_vector.connect([&](auto &data_out) {
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      interpretation(dim + 1,
                     DataComponentInterpretation::component_is_part_of_vector);

    interpretation[dim] = DataComponentInterpretation::component_is_scalar;

    data_out.add_data_vector(this->locally_relevant_block_solution,
                             this->component_names,
                             DataOut<dim>::type_dof_data,
                             interpretation);
  });
}



template <int dim>
void
Stokes<dim>::assemble_system_one_cell(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  ScratchData &                                         scratch,
  CopyData &                                            copy)
{
  auto &cell_matrix = copy.matrices[0];
  auto &cell_rhs    = copy.vectors[0];

  cell->get_dof_indices(copy.local_dof_indices[0]);

  const auto &fe_values = scratch.reinit(cell);
  cell_matrix           = 0;
  cell_rhs              = 0;

  for (const unsigned int q_index : fe_values.quadrature_point_indices())
    {
      for (const unsigned int i : fe_values.dof_indices())
        {
          const auto eps_v = fe_values[velocity].symmetric_gradient(
            i, q_index); // SymmetricTensor<2,dim>
          const auto div_v =
            fe_values[velocity].divergence(i, q_index);         // double
          const auto q = fe_values[pressure].value(i, q_index); // double

          for (const unsigned int j : fe_values.dof_indices())
            {
              const auto eps_u = fe_values[velocity].symmetric_gradient(
                j, q_index); // SymmetricTensor<2,dim>
              const auto div_u =
                fe_values[velocity].divergence(j, q_index);         // double
              const auto p = fe_values[pressure].value(j, q_index); // double

              cell_matrix(i, j) +=
                (scalar_product(eps_v, eps_u) - p * div_v - q * div_u) *
                fe_values.JxW(q_index); // dx
            }
          for (const unsigned int i : fe_values.dof_indices())
            {
              const auto comp_i = this->fe->system_to_component_index(i).first;
              cell_rhs(i) +=
                (fe_values.shape_value(i, q_index) * // phi_i(x_q)
                 this->forcing_term.value(fe_values.quadrature_point(q_index),
                                          comp_i) * // f(x_q)
                 fe_values.JxW(q_index));           // dx
            }
        }
    }

  if (cell->at_boundary())
    //  for(const auto face: cell->face_indices())
    for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
      if (this->neumann_ids.find(cell->face(f)->boundary_id()) !=
          this->neumann_ids.end())
        {
          auto &fe_face_values = scratch.reinit(cell, f);
          for (const unsigned int q_index :
               fe_face_values.quadrature_point_indices())
            for (const unsigned int i : fe_face_values.dof_indices())
              {
                const auto comp_i =
                  this->fe->system_to_component_index(i).first;
                cell_rhs(i) +=
                  fe_face_values.shape_value(i, q_index) *
                  this->neumann_boundary_condition.value(
                    fe_face_values.quadrature_point(q_index), comp_i) *
                  fe_face_values.JxW(q_index);
              }
        }
}


template <int dim>
void
Stokes<dim>::solve()
{
  TimerOutput::Scope                timer_section(this->timer, "solve");
  SolverGMRES<LA::MPI::BlockVector> solver(this->solver_control);


  LA::MPI::PreconditionAMG amg;
  amg.initialize(this->system_block_matrix.block(0, 0));


  solver.solve(this->system_block_matrix,
               this->block_solution,
               this->system_block_rhs,
               PreconditionIdentity());
  this->constraints.distribute(this->block_solution);
  this->locally_relevant_block_solution = this->block_solution;
}



template class Stokes<1>;
template class Stokes<2>;
template class Stokes<3>;