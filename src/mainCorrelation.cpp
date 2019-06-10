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
  double a = 3;
  int i_initial = 3;
  double gaussNorm_fixed = 0.01;
  int initial_fixed = 4;
  // Determines the particular program to run
  // 0 = loop dark matter, 1 = loop gaussian,
  // 2 = fixed dark matter, 3 = fixed gaussian
  int determine_program = 0;

  // Loop for dark spectrum condition
  if (determine_program == 0)
  {
    while ( i_initial < 11)
    {
      for (int j= 0; j< 10; j++)
      {
        testPowerSpectrum power_spectrum (&cosmological_model, 8.0, &filter);
        power_spectrum.setInitialCondition(i_initial);
        KFT::kftCosmology C (&cosmological_model, &power_spectrum);
        double a_loop = astro::x_logarithmic (double(j), 10.0, 3.0,5.0);
        power_spectrum.writeAllSpectrum(&C,a_loop);
        
        std::cout << "Stop computing, arrived: n = "<<  i_initial<< " a = " 
          << a_loop  << std::endl;
        
        usleep(10000000); // sleeps 10 seconds
        std::cout << "Start computing" << std::endl;
      }
      i_initial=i_initial+1;
    }
  }
  // Loop for gaussian initial condition
  if (determine_program == 1)
  {
    for (int i=1; i< 10; i++)
      {
        for (int j = 0; j< 9; j++)
        {
          double a = astro::x_logarithmic (double(j), 10.0, 0.005,0.1);
          testPowerSpectrum power_spectrum (&cosmological_model, 8.0, &filter);
          power_spectrum.setInitialCondition(0);
          double gaussNorm = astro::x_logarithmic (double(i), 10.0, 0.001,2.0);
          power_spectrum.setGaussNorm(gaussNorm);
          KFT::kftCosmology C (&cosmological_model, &power_spectrum);
          power_spectrum.writeAllGaussian(&C,a);
          std::cout << "Stop computing, arrived: normalization = "<< gaussNorm 
            << " a = " << a  << std::endl;
          usleep(10000000); // sleeps 10 seconds
          std::cout << "Start computing" << std::endl;
        }
      }
  }

  // Dark spectrum for fixed condition  
  if (determine_program == 2)
  {
    testPowerSpectrum power_spectrum (&cosmological_model, 8.0, &filter);
    power_spectrum.setInitialCondition(initial_fixed);
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
