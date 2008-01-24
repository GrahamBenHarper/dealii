//----------------------------  full_tensor_08.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  full_tensor_08.cc  ---------------------------

// a stress test using repeated multiplication of tensors

#include "../tests.h"
#include <base/tensor.h>
#include <base/logstream.h>
#include <fstream>
#include <iomanip>

  

template <int dim>
void test ()
{
  const double lambda = 1.5,
	       mu     = 1.7;
  Tensor<2,dim> t;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      t[i][j] = (1. + (i+lambda)*(mu+13));
	  
  Tensor<2,dim> a;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      a[i][j] = (1. + (i+1)*(j+1));

				   // stress test the whole thing many
				   // times. normalize in each step to
				   // make sure the result remains
				   // representable in floating point
				   // arithmetic. essentially, this
				   // invokes the power method to
				   // compute the largest eigenvector
				   // (eigentensor in this case)
  for (unsigned int i=0; i<1000000; ++i)
    {
      a = t*a;
      a /= a.norm();
    }

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j) 
      deallog << i << ' ' << j << ' ' << a[i][j] << std::endl;
}

  


int main ()
{
  std::ofstream logfile("full_tensor_08/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<2> ();
  test<3> ();
  
  deallog << "OK" << std::endl;
}
