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
  double a = 1.0;
  // Initial condition for the loop
  int i_initial = 15;
  double k0_fixed = 10.0;
  double sigma_fixed = 0.01;
  int initial_fixed = 18;
  // Determines the particular program to run
  // 0 = loop dark matter, 1 = loop gaussian,
  // 2 = fixed dark matter, 3 = fixed gaussian, 4 = gaussian loop sigma only
  // 5 = gaussian loop sigma and a only
  int determine_program = 2;
  
  if (determine_program == 0)
  {
    while ( i_initial < 23)
    {
      for (int j= 0; j< 5; j++)
      {
        double a_loop = astro::x_logarithmic (double(j), 5.0, 0.1,1.0);
        testPowerSpectrum power_spectrum (&cosmological_model, 8.0, &filter,
          i_initial, k0_fixed, sigma_fixed);
        // power_spectrum.setInitialCondition(i_initial);
        KFT::kftCosmology C (&cosmological_model, &power_spectrum);
        power_spectrum.writeAllSpectrum(&C,a_loop);
        
        std::cout << "Stop computing, arrived: n = " <<  i_initial << " a = " 
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
    for (int i=0; i< 2; i++)
      {
        for (int j = 0; j< 10; j++)
        {
          for (int m = 0; m < 10; m++)
          {
            double a_loop = astro::x_logarithmic (double(i), 2.0, 0.005,0.1);
            double k0_loop = astro::x_logarithmic (double(j), 10.0, 0.005,0.1);
            double sigma_loop = astro::x_logarithmic (double(m), 10.0, 0.005,0.1);
            testPowerSpectrum power_spectrum (&cosmological_model, 8.0, &filter, 0,
              k0_loop, sigma_loop);
            KFT::kftCosmology C (&cosmological_model, &power_spectrum);
            power_spectrum.writeAllGaussian(&C,a_loop);
            std::cout << "Stop computing, arrived: k0= "<< k0_loop << " sigma "
              << sigma_loop << " a = " << a_loop  << std::endl;
            usleep(10000000); // sleeps 10 seconds
            std::cout << "Start computing" << std::endl;
          }
        }
      }
  }

  // Dark spectrum for fixed condition  
  if (determine_program == 2)
  {
    testPowerSpectrum power_spectrum (&cosmological_model, 8.0, &filter, 
      initial_fixed, k0_fixed, sigma_fixed);
    KFT::kftCosmology C (&cosmological_model, &power_spectrum);
    power_spectrum.writeAllSpectrum(&C,a);
  }

  // Gaussian initial for fixed condition 
  if (determine_program == 3)
  {
    testPowerSpectrum power_spectrum (&cosmological_model, 8.0, &filter, 0,
      k0_fixed, sigma_fixed);
    KFT::kftCosmology C (&cosmological_model, &power_spectrum);
    power_spectrum.writeAllGaussian(&C,a);
  }

  if (determine_program == 4)
  {
    for (int i = 0; i<10; i++)
    {
      double sigma_loop = astro::x_logarithmic (double(i), 10.0, 0.1,10.0);
      testPowerSpectrum power_spectrum (&cosmological_model, 8.0, &filter, 0,
        k0_fixed, sigma_loop);
      KFT::kftCosmology C (&cosmological_model, &power_spectrum);
      power_spectrum.writeAllGaussian(&C,a);
      std::cout << "Stop computing, arrived: k0= "<< k0_fixed << " sigma "
        << sigma_loop << " a = " << a  << std::endl;
      usleep(10000000); // sleeps 10 seconds
      std::cout << "Start computing" << std::endl;
    }
  }
  
  if (determine_program == 5)
  {
    for (int i = 0; i<8; i++)
    {
      for (int j=0; j< 5; j++)
      {
        double sigma_loop = astro::x_logarithmic (double(i), 8.0, 5.0,10.0);
        double a_loop = astro::x_logarithmic (double(j), 5.0, 0.01,1.0);
        testPowerSpectrum power_spectrum (&cosmological_model, 8.0, &filter, 0,
          k0_fixed, sigma_loop);
        KFT::kftCosmology C (&cosmological_model, &power_spectrum);
        power_spectrum.writeAllGaussian(&C,a_loop);
        std::cout << "Stop computing, arrived: k0= "<< k0_fixed << " sigma "
          << sigma_loop << " a = " << a_loop  << std::endl;
        usleep(10000000); // sleeps 10 seconds
        std::cout << "Start computing" << std::endl;
      }
    }
  }

    return 0;
}

