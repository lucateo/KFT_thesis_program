#include "../include/KFT2.h"

int main ()
{
  astro::cosmologyBase cosmological_model;
  cosmological_model.setDarkUniverse ();
  astro::tophatFilter filter;
  
  int number_initial = 0;
  double gaussNorm = 0.1;
  
  testPowerSpectrum power_spectrum (&cosmological_model, 8.0, &filter);
  double a = 0.1;
  
  power_spectrum.setInitialCondition(number_initial);
  power_spectrum.setGaussNorm(gaussNorm);
  power_spectrum.trial(a);
  power_spectrum.writeTest(a);
 
  return 0;
}
