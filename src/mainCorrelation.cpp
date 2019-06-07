#include "../include/KFT2.h"
#include <unistd.h>

/**
 * automatizing correlation computation 
   */
int main ()
{
  astro::cosmologyBase cosmological_model;
  cosmological_model.setDarkUniverse ();
  astro::tophatFilter filter;

  int number_a = 10;
  int number_initial = 100;
  int i = 21; 
  
  while ( i< number_initial)
  {
        // for (int j=number_a -1; j< number_a; j++)
        // {
            testPowerSpectrum power_spectrum (&cosmological_model, 8.0, &filter);
            power_spectrum.setInitialCondition(i);
            KFT::kftCosmology C (&cosmological_model, &power_spectrum);
            // double a = astro::x_logarithmic (double(j), double(number_a), 0.01,1.0);
            power_spectrum.writeAllSpectrum(&C,1.0);
            std::cout << "Stop computing, arrived: n = "<<  i<< " a = " << 1  << std::endl;
            usleep(10000000); // sleeps 10 seconds
            std::cout << "Start computing" << std::endl;
        // }
        i=i+2;
    }
  
  // for gaussian initial condition 
  for (int i=0; i< 30; i++)
        {
            testPowerSpectrum power_spectrum (&cosmological_model, 8.0, &filter);
            power_spectrum.setInitialCondition(0);
            double gaussNorm = astro::x_logarithmic (double(i), 30.0, 0.001,100.0);
            power_spectrum.setGaussNorm(gaussNorm);
            KFT::kftCosmology C (&cosmological_model, &power_spectrum);
            power_spectrum.writeAllGaussian(&C,1.0);
            std::cout << "Stop computing, arrived: normalization = "<< gaussNorm << " a = " << 1  << std::endl;
            usleep(10000000); // sleeps 10 seconds
            std::cout << "Start computing" << std::endl;
        }
  return 0;
}
