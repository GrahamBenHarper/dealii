// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/std_cxx14/memory.h>

#include <deal.II/fe/fe_tensor_product_polynomials.h>

#include <deal.II/lac/vector.h>

#include <algorithm>
#include <sstream>


DEAL_II_NAMESPACE_OPEN


template <int dim, int spacedim>
FE_TensorProductPolynomials<dim, spacedim>::FE_TensorProductPolynomials(
  const unsigned int order,
  const unsigned int n_face_support_points)
  : FE_Poly<TensorProductPolynomials<dim>, dim>(
      TensorProductPolynomials<dim>(),
      FiniteElementData<dim>(this->get_dpo_vector(),
                             1,
                             2,
                             FiniteElementData<dim>::L2),
      std::vector<bool>(4, false), // restriction not implemented
      std::vector<ComponentMask>(4, std::vector<bool>(1, true)))
  , order(order)
  , n_face_support_points(n_face_support_points)
{
  Assert(dim == 2, ExcNotImplemented());
  Assert(order == 0, ExcNotImplemented());
  this->initialize_support_points();
}



template <int dim, int spacedim>
std::vector<unsigned int>
FE_TensorProductPolynomials<dim, spacedim>::get_dpo_vector()
{
  std::vector<unsigned int> dpo(dim + 1, 0);
  dpo[dim - 1] = 1;

  return dpo;
}



template <int dim>
std::string
FE_TensorProductPolynomials<dim, spacedim>::get_name() const
{
  std::ostringstream namebuf;
  namebuf << "FE_TensorProductPolynomials"
          << "<" << dim << ">"
          << "(" << this->order << ", " << this->n_face_support_points << ")";
  return namebuf.str();
}



template <int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, spacedim>>
FE_TensorProductPolynomials<dim>::clone() const
{
  return std_cxx14::make_unique<FE_TensorProductPolynomials<dim, spacedim>>(
    this->order, this->n_face_support_points);
}



template <int dim, int spacedim>
void
FE_TensorProductPolynomials<dim, spacedim>::initialize_support_points()
{
  Assert(dim == 2, ExcNotImplemented());
  dealii::QGauss<dim - 1> face_quadrature(this->n_face_support_points);
  this->weights = face_quadrature.get_weights();
  this->generalized_support_points.resize(4 * face_quadrature.size());
  for (unsigned int q = 0; q < face_quadrature.size(); ++q)
    {
      this->generalized_support_points[0 * face_quadrature.size() + q] =
        dealii::Point<dim>(0, 1 - face_quadrature.point(q)(0));
      this->generalized_support_points[1 * face_quadrature.size() + q] =
        dealii::Point<dim>(1, 1 - face_quadrature.point(q)(0));
      this->generalized_support_points[2 * face_quadrature.size() + q] =
        dealii::Point<dim>(face_quadrature.point(q)(0), 0);
      this->generalized_support_points[3 * face_quadrature.size() + q] =
        dealii::Point<dim>(face_quadrature.point(q)(0), 1);
    }
}



template <int dim, int spacedim>
void
FE_TensorProductPolynomials<dim, spacedim>::convert_generalized_support_point_values_to_dof_values(
  const std::vector<Vector<double>> &support_point_values,
  std::vector<double> &              nodal_values) const
{
  AssertDimension(support_point_values.size(),
                  this->generalized_support_points.size());
  AssertDimension(nodal_values.size(), this->dofs_per_cell);

  const unsigned int q_points_per_face = this->weights.size();
  std::fill(nodal_values.begin(), nodal_values.end(), 0.0);

  std::vector<Vector<double>>::const_iterator value =
    support_point_values.begin();
  for (unsigned int face = 0; face < dealii::GeometryInfo<dim>::faces_per_cell;
       ++face)
    {
      for (unsigned int q = 0; q < q_points_per_face; ++q)
        {
          nodal_values[face] += (*value)[0] * this->weights[q];
          ++value;
        }
    }
}




template <>
void
FE_TensorProductPolynomials<1, 2>::fill_fe_values(
  const Triangulation<1, 2>::cell_iterator &,
  const CellSimilarity::Similarity       cell_similarity,
  const Quadrature<1> &                  quadrature,
  const Mapping<1, 2> &                  mapping,
  const Mapping<1, 2>::InternalDataBase &mapping_internal,
  const dealii::internal::FEValuesImplementation::MappingRelatedData<1, 2>
    &                                          mapping_data,
  const FiniteElement<1, 2>::InternalDataBase &fe_internal,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<1, 2>
    &output_data) const
{
  // convert data object to internal data for this class. fails with an
  // exception if that is not possible
  Assert(dynamic_cast<const InternalData *>(&fe_internal) != nullptr,
         ExcInternalError());
  const InternalData &fe_data =
    static_cast<const InternalData &>(fe_internal); // NOLINT

  // transform gradients and higher derivatives. there is nothing to do
  // for values since we already emplaced them into output_data when
  // we were in get_data()
  if (fe_data.update_each & update_gradients &&
      cell_similarity != CellSimilarity::translation)
    for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
      mapping.transform(make_array_view(fe_data.shape_gradients, k),
                        mapping_covariant,
                        mapping_internal,
                        make_array_view(output_data.shape_gradients, k));

  if (fe_data.update_each & update_hessians &&
      cell_similarity != CellSimilarity::translation)
    {
      for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
        mapping.transform(make_array_view(fe_data.shape_hessians, k),
                          mapping_covariant_gradient,
                          mapping_internal,
                          make_array_view(output_data.shape_hessians, k));

      for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
        for (unsigned int i = 0; i < quadrature.size(); ++i)
          for (unsigned int j = 0; j < 2; ++j)
            output_data.shape_hessians[k][i] -=
              mapping_data.jacobian_pushed_forward_grads[i][j] *
              output_data.shape_gradients[k][i][j];
    }

  if (fe_data.update_each & update_3rd_derivatives &&
      cell_similarity != CellSimilarity::translation)
    {
      for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
        mapping.transform(make_array_view(fe_data.shape_3rd_derivatives, k),
                          mapping_covariant_hessian,
                          mapping_internal,
                          make_array_view(output_data.shape_3rd_derivatives,
                                          k));

      for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
        correct_third_derivatives(output_data,
                                  mapping_data,
                                  quadrature.size(),
                                  k);
    }
}



template <>
void
FE_TensorProductPolynomials<2, 3>::fill_fe_values(
  const Triangulation<2, 3>::cell_iterator &,
  const CellSimilarity::Similarity       cell_similarity,
  const Quadrature<2> &                  quadrature,
  const Mapping<2, 3> &                  mapping,
  const Mapping<2, 3>::InternalDataBase &mapping_internal,
  const dealii::internal::FEValuesImplementation::MappingRelatedData<2, 3>
    &                                          mapping_data,
  const FiniteElement<2, 3>::InternalDataBase &fe_internal,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<2, 3>
    &output_data) const
{
  // convert data object to internal data for this class. fails with an
  // exception if that is not possible
  Assert(dynamic_cast<const InternalData *>(&fe_internal) != nullptr,
         ExcInternalError());
  const InternalData &fe_data =
    static_cast<const InternalData &>(fe_internal); // NOLINT

  // transform gradients and higher derivatives. there is nothing to do
  // for values since we already emplaced them into output_data when
  // we were in get_data()
  if (fe_data.update_each & update_gradients &&
      cell_similarity != CellSimilarity::translation)
    for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
      mapping.transform(make_array_view(fe_data.shape_gradients, k),
                        mapping_covariant,
                        mapping_internal,
                        make_array_view(output_data.shape_gradients, k));

  if (fe_data.update_each & update_hessians &&
      cell_similarity != CellSimilarity::translation)
    {
      for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
        mapping.transform(make_array_view(fe_data.shape_hessians, k),
                          mapping_covariant_gradient,
                          mapping_internal,
                          make_array_view(output_data.shape_hessians, k));

      for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
        for (unsigned int i = 0; i < quadrature.size(); ++i)
          for (unsigned int j = 0; j < 3; ++j)
            output_data.shape_hessians[k][i] -=
              mapping_data.jacobian_pushed_forward_grads[i][j] *
              output_data.shape_gradients[k][i][j];
    }

  if (fe_data.update_each & update_3rd_derivatives &&
      cell_similarity != CellSimilarity::translation)
    {
      for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
        mapping.transform(make_array_view(fe_data.shape_3rd_derivatives, k),
                          mapping_covariant_hessian,
                          mapping_internal,
                          make_array_view(output_data.shape_3rd_derivatives,
                                          k));

      for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
        correct_third_derivatives(output_data,
                                  mapping_data,
                                  quadrature.size(),
                                  k);
    }
}



// explicit instantiations
#include "fe_tensor_product_polynomials.inst"

DEAL_II_NAMESPACE_CLOSE
