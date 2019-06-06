#include <stdexcept>
#include <limits>

#include <astro/utilities/osc_integrator.h>
#include <astro/utilities/integrator.h>

namespace astro
{
  osc_integrator::osc_integrator
    (Integrand & kernel_, astro::osc_function_type osc_type_, int nu_,
     double lower_bound_, double upper_bound_, astro::spacing_type bin_type_,
     int n_bins_):
  kernel (kernel_), osc_type (osc_type_), bessel_order (nu_),
  lower_bound (lower_bound_), upper_bound (upper_bound_),
  bin_type (bin_type_), n_bins (n_bins_)
  {
    if (!(bin_type == astro::LINEAR_SPACING) && 
        !(bin_type == astro::LOG_SPACING))
    {
      std::string errorMsg ("In astro::osc_integrator: invalid binning type");
      throw std::runtime_error (errorMsg);
    }

    if ((bin_type == astro::LOG_SPACING) &&
        ((lower_bound < 0.0) || (upper_bound < 0.0)))
    {
      std::string errorMsg
        ("In astro::osc_integrator: integral bounds must be positive reals");
      throw std::runtime_error (errorMsg);
    }
    
    n_abscissae = 8;

    if ((bin_type == astro::LOG_SPACING) && (lower_bound == 0.0))
      lower_bound = 1.0e-10;
    bins = astro::fill (bin_type, n_bins, lower_bound, upper_bound);
  }
  
  double osc_integrator::CQUAD_kernel (double k_, double x_)
  {
    double result = 0.0;
    
    if (bin_type == astro::LINEAR_SPACING)
      result = LC_oscillator (k_*x_, osc_type, bessel_order) * kernel (k_,x_);
    
    if (bin_type == astro::LOG_SPACING)
    {
      double y = exp (x_);
      result = y * LC_oscillator (k_*y, osc_type, bessel_order) *
        kernel (k_,y);
    }
    return result;
  }
  
  double osc_integrator::operator ()(double k_, bool warning_)
  {
    double max_periods = 2.0;
    double result = 0.0;

    // #pragma omp parallel for reduction (+:result)
    for (size_t i = 0; i < bins.size ()-1; i++)
    {  
      double res = 0.0;
      if (k_*(bins.at (i+1)-bins.at (i)) < max_periods * 2.0*M_PI)
      {
        std::function<double (double)> kernel = std::bind
          (&osc_integrator::CQUAD_kernel, this, k_, std::placeholders::_1);
        astro::integrator gsl_integrator (kernel);

        double lower_bin_bound = bins.at (i);
        double upper_bin_bound = bins.at (i+1);

        if (bin_type == astro::LOG_SPACING)
        {
          double eps = std::numeric_limits<double>::epsilon ();
          lower_bin_bound = log (lower_bin_bound)*(1.0+eps);
          upper_bin_bound = log (upper_bin_bound)*(1.0-eps);
        }
        
        res = gsl_integrator (lower_bin_bound, upper_bin_bound, true);
      }
      else
      {
        astro::LC_integrator integrate
          (kernel, osc_type, bessel_order, bins.at (i), bins.at (i+1),
           n_abscissae);
        res = integrate (k_, warning_);
      }

      result += res;
    }
    return result;
  }  
}
