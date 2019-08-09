// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2019 by the deal.II authors
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

//NOTE: The required calls for FE_PolyTensor seem to be
// poly_space.compute(), poly_space.n(), and poly_space.degree()
// poly_space.n() will be filled by the function compute_n_pols(k).

#ifndef dealii_polynomials_base_h
#define dealii_polynomials_base_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomial_space.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * This class provides a framework for the finite element polynomial
 * classes for use with finite element classes that are derived from
 * @p FE_Poly and @p FE_PolyTensor. An object of this type is stored
 * as a member variable in the finite element class. 
 *
 * <h3>Deriving classes</h3>
 *
 * Any derived class must provide a the most basic properties for shape 
 * functions evaluated on the reference cell. 
 *
 * @ingroup Polynomials
 * @author Graham Harper
 * @date 2019
 */
template <int dim>
class PolynomialsBase
{
  /**
   * Constructor. 
   */
  PolynomialsBase(const unsigned int k);
  
  /**
    * Move constructor.
    */
  PolynomialsBase(PolynomialsBase<dim> &&) = default; // NOLINT

  /**
   * Copy constructor.
   */
  PolynomialsBase(const PolynomialsBase<dim> &) = default;

  /**
   * Virtual destructor. Makes sure that pointers to this class are deleted
   * properly.
   */
  virtual ~PolynomialsBase() = default;
  
  /**
   * Compute the value and the derivatives of the polynomials at 
   * @p unit_point.
   *
   * The size of the vectors must either be zero or equal <tt>n()</tt>.  In
   * the first case, the function will not compute these values.
   *
   * If you need values or derivatives of all polynomials then use this
   * function, rather than using any of the <tt>compute_value</tt>, 
   * <tt>compute_grad</tt> or <tt>compute_grad_grad</tt> functions, see below,
   * in a loop over all tensor product polynomials.
   */
  virtual void
  compute(const Point<dim> &           unit_point,
          std::vector<Tensor<1, dim>> &values,
          std::vector<Tensor<2, dim>> &grads,
          std::vector<Tensor<3, dim>> &grad_grads,
          std::vector<Tensor<4, dim>> &third_derivatives,
          std::vector<Tensor<5, dim>> &fourth_derivatives) const = 0;

  /**
   * Return the number of polynomials.
   */
  virtual unsigned int
  n() const;

  /**
   * Return the degree of the space.
   */
  virtual unsigned int
  degree() const;

  /**
   * A sort of virtual copy constructor, this function returns a copy of
   * the polynomial space object. Derived classes need to override the function
   * here in this base class and return an object of the same type as the
   * derived class.
   *
   * Some places in the library, for
   * example the constructors of FE_Poly and FE_PolyTensor need to make copies
   * of polynomial spaces without knowing their exact type. They do so through
   * this function.
   */
  virtual std::unique_ptr<PolynomialsBase<dim>>
  clone() const = 0;
  
  /**
   * Return the name of the space.
   */
  virtual std::string
  get_name() const = 0;

  /**
   * Return the number of polynomials in the space. This must be overridden by
   * the derived class.
   */
  virtual unsigned int
  compute_n_pols(const unsigned int k) const;

private:
  /**
   * The degree of this object as given to the constructor.
   */
  const unsigned int my_degree;

  /**
   * Number of polynomials.
   */
  const unsigned int n_pols;
};



template <int dim>
inline unsigned int
PolynomialsBase<dim>::n() const
{
  return n_pols;
}



template <int dim>
inline unsigned int
PolynomialsBase<dim>::degree() const
{
  return my_degree;
}



DEAL_II_NAMESPACE_CLOSE

#endif
