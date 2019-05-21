#include <iostream>
#include <astro/astro.h>
#include <astro/utilities/besselInt.h>
#include <astro/utilities/utilities.h>

/*
#include "../src/astro.h"
#include "../src/utilities/besselInt.h"
#include "../src/utilities/utilities.h" 
*/

double kernel (double x, void * p)
{ return x; }

int main ()
{
  gsl_function gsl_fct;
  gsl_fct.params = NULL;
  gsl_fct.function = &kernel;
  astro::besselInt bessel_integral (&gsl_fct, 1, 2);
  
  const int n = 256;
  const double k_min = 0.1, k_max = 10.0;
  for (int i = 0; i < n; i++)
  {
    double k = astro::x_logarithmic (i, n, k_min, k_max);
    std::cout
      << k << "  "
      << bessel_integral (1.0, 10.0, k) << std::endl;
  }
  
  return 0;
}
