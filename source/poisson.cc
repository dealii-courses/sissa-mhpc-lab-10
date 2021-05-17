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

using namespace dealii;

template <int dim>
Poisson<dim>::Poisson()
  : BaseProblem<dim>(1, "Poisson<" + std::to_string(dim) + ">")
{
  this->add_parameter("Coefficient expression", coefficient_expression);

  // Make sure we initialize the coefficients.
  this->setup_system_call_back.connect([&]() {
    const auto vars = dim == 1 ? "x" : dim == 2 ? "x,y" : "x,y,z";
    coefficient.initialize(vars, coefficient_expression, this->constants);
  });

  // Output the scalar result.
  this->add_data_vector.connect([&](auto &data_out) {
    data_out.add_data_vector(this->locally_relevant_solution, "solution");
  });
}



template <int dim>
void
Poisson<dim>::assemble_system_one_cell(
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
        for (const unsigned int j : fe_values.dof_indices())
          cell_matrix(i, j) +=
            (coefficient.value(fe_values.quadrature_point(q_index)) * // a(x_q)
             fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
             fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
             fe_values.JxW(q_index));           // dx
      for (const unsigned int i : fe_values.dof_indices())
        cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q)
                        this->forcing_term.value(
                          fe_values.quadrature_point(q_index)) * // f(x_q)
                        fe_values.JxW(q_index));                 // dx
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
              cell_rhs(i) += fe_face_values.shape_value(i, q_index) *
                             this->neumann_boundary_condition.value(
                               fe_face_values.quadrature_point(q_index)) *
                             fe_face_values.JxW(q_index);
        }
}

template class Poisson<2>;
template class Poisson<3>;