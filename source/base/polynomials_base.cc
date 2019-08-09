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


#include <deal.II/base/polynomials_base.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/thread_management.h>

#include <iomanip>
#include <iostream>

DEAL_II_NAMESPACE_OPEN

template <int dim>
PolynomialsBase<dim>::PolynomialsBase(const unsigned int k)
  : my_degree(k)
  , n_pols(compute_n_pols(k))
{
  // nothing to do here for now  
}


template <int dim>
unsigned int
PolynomialsBase<dim>::compute_n_pols(unsigned int /*k*/) const
{
  Assert(false, ExcNotImplemented());
  return 0;
}



template class PolynomialsBase<1>;
template class PolynomialsBase<2>;
template class PolynomialsBase<3>;

DEAL_II_NAMESPACE_CLOSE
