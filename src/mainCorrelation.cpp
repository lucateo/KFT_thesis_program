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
  double a = 1;
  // Initial condition for the loop
  int i_initial = 310;
  double k0_fixed = 10.0;
  double sigma_fixed = 0.215443 ;
  int initial_fixed = 4;
  // Determines the particular program to run
  // 0 = loop dark matter, 1 = loop gaussian,
  // 2 = fixed dark matter, 3 = fixed gaussian, 4 = gaussian loop sigma only
  // 5 = gaussian loop sigma and a only, 6 = trial higher order
  // 7 = control spectrum function, 8 = Trial function for higher orders 
  // 9 = curly P_ij trial function
  int determine_program = 9;

  if (determine_program == 0)
  {
    while ( i_initial < 400)
    {
      for (int j= 4; j< 5; j++)
      {
        double a_loop = astro::x_logarithmic (double(j), 5.0, 0.01,1.0);
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
      i_initial=i_initial+20;
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
    for (int i = 0; i<10; i++)
    {
      for (int j=0; j< 3; j++)
      {
        double sigma_loop = astro::x_logarithmic (double(i), 10.0, 8.0,30.0);
        double a_loop = astro::x_logarithmic (double(j), 3.0, 0.01,1.0);
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

  if (determine_program == 6)
  {
    testPowerSpectrum power_spectrum (&cosmological_model, 8.0, &filter, 
        initial_fixed, k0_fixed, sigma_fixed);
    KFT::kftCosmology C (&cosmological_model, &power_spectrum);
    double ratio = 2.0;
    double a_local = 1.0;
    int determine = 1;
    int i_initial = 0;
    power_spectrum.writeBiSpectrumFixExpDamping(&C , a_local, ratio,determine,i_initial);
  }

  if (determine_program == 7)
  {
    testPowerSpectrum power_spectrum (&cosmological_model, 8.0, &filter,
      420, k0_fixed, sigma_fixed);
    KFT::kftCosmology C (&cosmological_model, &power_spectrum);
    double a_local = 1.0;
    power_spectrum.writeAllSpectrumControl(&C,a_local);
  }
  
  if (determine_program == 8)
  {
    testPowerSpectrum power_spectrum (&cosmological_model, 8.0, &filter,
      initial_fixed, k0_fixed, sigma_fixed);
    KFT::kftCosmology C (&cosmological_model, &power_spectrum);
    std::ostringstream os_ps, os_cf ;
    os_ps << "data/ps_table/ps_tableHigherOrder_a_" << a << "_n_initial_"
      << initial_fixed << ".d";

    os_cf << "data/cf_table/cf_tableHigherOrder_a_"<< a << "_n_initial_"
      << initial_fixed << ".d";

    std::string ps_table = os_ps.str();
    std::string cf_table = os_cf.str();
    powerSpectraModified P (&C);
    P.initCorrelation (ps_table, cf_table);
    double a_local = 0.01;
    double k = 1.0;
    P.Trial(a_local,k);
  }
  
  if (determine_program == 9)
  {
    testPowerSpectrum power_spectrum (&cosmological_model, 8.0, &filter,
      4, k0_fixed, sigma_fixed);
    KFT::kftCosmology C (&cosmological_model, &power_spectrum);
    double a_local = 1.0;
    double k_prime = 1.0;
    double mu = 0.5;
    int determine = 2;
    int i_initial =75;
    power_spectrum.writeAllHigherOrder(&C,a_local, k_prime,mu, determine, i_initial);
  }
  return 0;
}

