#include "../include/KFT2.h"

int main ()
{
  astro::cosmologyBase cosmological_model;
  cosmological_model.setDarkUniverse ();
  astro::tophatFilter filter;
  
  int number_initial = 9;
  double k0_gauss = 0.1;
  double sigma_gauss = 1.0;
  
  testPowerSpectrum power_spectrum (&cosmological_model, 8.0, &filter, 
    number_initial, k0_gauss, sigma_gauss);
  double a = 1.0;
  
  power_spectrum.trial(a);
  power_spectrum.writeTest(a);
 
  return 0;
}
