#include "../include/Bartelmann.h"

testPowerSpectrum::testPowerSpectrum
    (astro::cosmologyBase * cosmological_model_in, const double sigma8_in,
     astro::scaleFilter * filter_in):
  powerSpectrum (cosmological_model_in, sigma8_in, filter_in, 1.0)
  {
    cosmic_structures = new astro::cosmicStructures (cosmological_model);
    powerSpectrum::amplitude = sigma_8*sigma_8/sigma2 (8.0);
  }

  testPowerSpectrum::~testPowerSpectrum ()
  { delete cosmic_structures; }

  double testPowerSpectrum::operator () (const double k, const double a)
  {
    double D_plus = cosmic_structures->Dplus (a);
    double kappa = k/k0;
    return D_plus*D_plus*amplitude*k/gsl_pow_2 (1.0+kappa*kappa);
  }


