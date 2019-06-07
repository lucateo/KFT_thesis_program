#include "../include/KFT2.h"

int main ()
{
  astro::cosmologyBase cosmological_model;
  cosmological_model.setDarkUniverse ();
  astro::tophatFilter filter;

  int a = 1;
  int number_initial = 0;
  double gaussNorm = 1;
  testPowerSpectrum power_spectrum (&cosmological_model, 8.0, &filter);
  power_spectrum.setInitialCondition(number_initial);
  power_spectrum.setGaussNorm(gaussNorm);
  power_spectrum.trial(a);
 
  return 0;
}
