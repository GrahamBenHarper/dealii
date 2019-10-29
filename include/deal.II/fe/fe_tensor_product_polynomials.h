// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2019 by the deal.II authors
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


#ifndef dealii_fe_tensor_product_polynomials_h
#define dealii_fe_tensor_product_polynomials_h

#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/fe_base.h>
#include <deal.II/fe/fe_poly.h>

#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN


/**
 * Implementation of the TensorProductPolynomials element. 
 * ADD MORE
 *
 * @ingroup fe
 * @author Graham Harper
 * @date 2019
 */
template <int dim>
class FE_TensorProductPolynomials : public FE_Poly<TensorProductPolynomials<dim>, dim>
{
public:
  /**
   * ADD MORE
   */
  FE_TensorProductPolynomials(const unsigned int order                 = 0,
                    const unsigned int n_face_support_points = 2);

  virtual std::string
  get_name() const override;

  virtual std::unique_ptr<FiniteElement<dim, dim>>
  clone() const override;

  // documentation inherited from the base class
  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double> &              nodal_values) const override;

private:
  /**
   * Order of this element.
   */
  const unsigned int order;

  /**
   * The number of quadrature points used on each face to evaluate node
   * functionals during interpolation.
   */
  const unsigned int n_face_support_points;

  /**
   * The weights used on the faces to evaluate node functionals.
   */
  std::vector<double> weights;

  /**
   * Compute generalized support points and their weights.
   */
  void
  initialize_support_points();
  /**
   * Return information about degrees of freedom per object as needed during
   * construction.
   */
  std::vector<unsigned int>
  get_dpo_vector();
};


DEAL_II_NAMESPACE_CLOSE

#endif
