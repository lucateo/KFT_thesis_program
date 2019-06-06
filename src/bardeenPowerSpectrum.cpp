#include <iostream>
#include <astro/cosmology/bardeenPowerSpectrum.h>

const double astro::bardeenPowerSpectrum::spectral_parameters[6] =
  {2.34, 3.89, 16.1, 5.46, 6.71, 10.2};

astro::bardeenPowerSpectrum::bardeenPowerSpectrum
  (cosmologyBase * cosmological_model_in,
   double sigma_8_in, scaleFilter * sigma_filter_in, double index):
powerSpectrum (cosmological_model_in, sigma_8_in, sigma_filter_in, index)
{
  const double h=cosmological_model->hubbleFunction (1.0);
  const double omegaM=cosmological_model->omegaMatter (1.0);
  const double omegaB=cosmological_model->omegaBaryon (1.0);
  shape_parameter = omegaM*h*exp(-omegaB*(1.0+sqrt(2*h)/omegaM));
  amplitude = sigma_8*sigma_8/sigma2 (8.0);
  cosmic_structures = new cosmicStructures (cosmological_model);
}

astro::bardeenPowerSpectrum::~bardeenPowerSpectrum ()
{ delete cosmic_structures; }

double astro::bardeenPowerSpectrum::transfer_function (const double k)
{
  double q = k/shape_parameter;
  double aq = spectral_parameters[0]*q;
  double ts = 1.0+
    spectral_parameters[1]*q+
    gsl_pow_2 (spectral_parameters[2]*q)+
    gsl_pow_3 (spectral_parameters[3]*q)+
    gsl_pow_4 (spectral_parameters[4]*q);
  return log (1.0+aq)/aq/sqrt (sqrt (ts));
}

double astro::bardeenPowerSpectrum::operator ()
  (const double k, const double a)
{
  double t = transfer_function (k);
  double growth_factor = (a < 1.0) ?
    a*cosmic_structures->growthFactor (a)/
    cosmic_structures->growthFactor (1.0) : 1.0;
  t *= growth_factor;
  return amplitude*pow(k,spectral_index)*t*t;
}

double astro::bardeenPowerSpectrum::getShapeParameter () const
{ return shape_parameter; }
