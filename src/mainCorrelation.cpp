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

  // Parameters to look at
  int number_initial = 100;
  double a = 0.01;
  int i = 21;
  double gaussNorm_fixed = 10.0;
  // Determines the particular program to run
  int determine_program = 3;

  // Loop for dark spectrum condition
  if (determine_program == 0)
  {
    while ( i< number_initial)
    {
        // for (int j=number_a -1; j< number_a; j++)
        // {
            testPowerSpectrum power_spectrum (&cosmological_model, 8.0, &filter);
            power_spectrum.setInitialCondition(4);
            KFT::kftCosmology C (&cosmological_model, &power_spectrum);
            // double a = astro::x_logarithmic (double(j), double(number_a), 0.01,1.0);
            power_spectrum.writeAllSpectrum(&C,a);
            std::cout << "Stop computing, arrived: n = "<<  i<< " a = " << a  << std::endl;
            usleep(10000000); // sleeps 10 seconds
            std::cout << "Start computing" << std::endl;
        // }
        i=i+2;
    }
  }
  // Loop for gaussian initial condition
  if (determine_program == 1)
  {
    for (int i=0; i< 30; i++)
        {
            testPowerSpectrum power_spectrum (&cosmological_model, 8.0, &filter);
            power_spectrum.setInitialCondition(0);
            double gaussNorm = astro::x_logarithmic (double(i), 20.0, 10.0,1000.0);
            power_spectrum.setGaussNorm(gaussNorm);
            KFT::kftCosmology C (&cosmological_model, &power_spectrum);
            power_spectrum.writeAllGaussian(&C,a);
            std::cout << "Stop computing, arrived: normalization = "<< gaussNorm << " a = " << a  << std::endl;
            usleep(10000000); // sleeps 10 seconds
            std::cout << "Start computing" << std::endl;
        }
  }

  // Dark spectrum for fixed condition  
  if (determine_program == 2)
  {
        testPowerSpectrum power_spectrum (&cosmological_model, 8.0, &filter);
        power_spectrum.setInitialCondition(4);
        KFT::kftCosmology C (&cosmological_model, &power_spectrum);
        power_spectrum.writeAllSpectrum(&C,a);
  }

  // Gaussian initial for fixed condition 
  if (determine_program == 3)
  {
        testPowerSpectrum power_spectrum (&cosmological_model, 8.0, &filter);
        power_spectrum.setInitialCondition(0);
        power_spectrum.setGaussNorm(gaussNorm_fixed);
        KFT::kftCosmology C (&cosmological_model, &power_spectrum);
        power_spectrum.writeAllGaussian(&C,a);
  }
    return 0;
}
