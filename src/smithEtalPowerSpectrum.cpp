#include <astro/cosmology/smithEtalPowerSpectrum.h>
#include <astro/cosmology/gaussFilter.h>

const double astro::smithEtalPowerSpectrum::eps = 0.01;

astro::smithEtalPowerSpectrum::smithEtalPowerSpectrum
  (cosmologyBase * cosmological_model_in,
   double sigma_8_in, scaleFilter * sigma_filter_in, double index):
powerSpectrum (cosmological_model_in, sigma_8_in, sigma_filter_in, index)
{
  spectral_parameters.ao = -1.0;
  cosmic_structures = new cosmicStructures (cosmological_model);
  bardeen_ps = new bardeenPowerSpectrum
    (cosmological_model_in, sigma_8_in, sigma_filter_in, index);
}

astro::smithEtalPowerSpectrum::~smithEtalPowerSpectrum ()
{
  delete cosmic_structures;
  delete bardeen_ps;
}

double astro::smithEtalPowerSpectrum::operator () (double k, double a)
{
  if (fabs (a-spectral_parameters.ao) >= eps)
    setConstants (a);
  double y = k*spectral_parameters.rs;
  double delta = 0.5/M_PI/M_PI*gsl_pow_3 (k)*bardeen_ps->operator () (k,a);
  double dq = delta*
    (pow (1.0+delta,spectral_parameters.ac[1])/
     (1.0+spectral_parameters.ac[0]*delta))*exp (-0.25*y*(1.0+0.5*y));
  double dp = spectral_parameters.ac[5]*pow (y,3*spectral_parameters.om[1])/
    (1.0+spectral_parameters.ac[6]*pow (y,spectral_parameters.om[2])+
    pow (spectral_parameters.ac[7]*spectral_parameters.om[3]*y,
         3.0-spectral_parameters.ac[2]));
  double dh = dp/
    (1.0+spectral_parameters.ac[3]/y+spectral_parameters.ac[4]/y/y);
  double dn = dq+dh;
  return 2.0*M_PI*M_PI*dn/gsl_pow_3 (k);
}

void astro::smithEtalPowerSpectrum::setConstants (double a)
{
  spectral_parameters.ao = a;
  double growth_factor = a*
    cosmic_structures->growthFactor (a)/cosmic_structures->growthFactor (1.0);
  spectral_parameters.om[0] = cosmological_model->omegaMatter (a);
  spectral_parameters.om[1] = pow (spectral_parameters.om[0],-0.0307);
  spectral_parameters.om[2] = pow (spectral_parameters.om[0],-0.0585);
  spectral_parameters.om[3] = pow (spectral_parameters.om[0], 0.0743);

  scaleFilter * save_filter = bardeen_ps->getFilter ();
  gaussFilter gauss_filter;
  bardeen_ps->setFilter (&gauss_filter);  
  spectral_parameters.rs = bardeen_ps->nonlinearRadius (1.0/growth_factor);
  double n_eff = bardeen_ps->effectiveIndex (spectral_parameters.rs);
  double c_eff = bardeen_ps->effectiveCurvature (spectral_parameters.rs);
  bardeen_ps->setFilter (save_filter);
  
  double n_eff_2 = n_eff*n_eff;
  spectral_parameters.ac[0] = 1.3884+0.3700*n_eff-0.1452*n_eff_2;
  spectral_parameters.ac[1] = 0.8291+0.9854*n_eff+0.3401*n_eff_2;
  spectral_parameters.ac[2] = 0.8649+0.2989*n_eff+0.1631*c_eff;
  spectral_parameters.ac[3] = -3.5442+0.1908*n_eff;
  spectral_parameters.ac[4] = 0.9589+1.2857*n_eff;
  spectral_parameters.ac[5] = 1.4861+1.8369*n_eff+1.6762*n_eff_2+
    0.7940*n_eff*n_eff_2+0.1670*n_eff_2*n_eff_2-0.6206*c_eff;
  spectral_parameters.ac[6] = 0.9463+0.9466*n_eff+0.3084*n_eff_2-0.9400*c_eff;
  spectral_parameters.ac[7] = -0.2807+0.6669*n_eff+0.3214*n_eff_2-0.0793*c_eff;
  for (int i = 3; i < 8; i++)
    spectral_parameters.ac[i] = exp (spectral_parameters.ac[i]/M_LOG10E);
}
