#include "../include/Bartelmann.h"

/**
 * Main program, essentially simply putting together your power spectrum
 * with implemented KFT methods. The program accepts two command-line options:
 * -t {type}, with type = 0, 1, 2
 * -n {n}, with n specifying the number of output entries in the power
 * spectrum table
 * The types mean:
 *   0: write the input power spectrum to file "testSpectrum.txt"
 *   1: create and tabulate the momentum-correlation functions
 *   2: calculate power spectra from these input tables
 */
int main (int argc, char * argv[])
{
  astro::clArguments ca (argc, argv);
  const int type = ca.get ("-t", 0);
  const int n = ca.get ("-n", 128);
  
  /**
   * Set internal variables
   */
  double a_initial = 1.0e-3;
  const double k_min = 1.0e-3, k_max = 100.0;
  const double q_min = 1.0e-5, q_max = 1.0e5;
  std::string ps_table = "ps_table.d";
  std::string cf_table = "cf_table.d";
  
  /**
   * Initialise cosmological model, power spectrum, and KFT cosmology class
   */
  astro::cosmologyBase cosmological_model;
  cosmological_model.setDarkUniverse ();
  astro::tophatFilter filter;
  testPowerSpectrum power_spectrum (&cosmological_model, 8.0, &filter);
  KFT::kftCosmology C (&cosmological_model, &power_spectrum);
  
  /**
   * Choose calculation type
   */
  switch (type)
  {
    case (0):
    {
      /**
       * Tabulate input power spectrum to file
       */
      astro::functionWriter write ("testSpectrum.txt");
      write.push_back ([&] (double k) { return power_spectrum (k, 1.0); });
      write (k_min, k_max, n, astro::LOG_SPACING);
      break;
    }
    case (1):
    {
      /**
       * Calculate and tabulate momentum-correlation functions
       */
      KFT::iniCorrTable corr_table
        (C.get_power_spectrum (), a_initial, q_min, q_max);
      corr_table.print_tables (ps_table, cf_table);
      break;
    }
    case (2):
    {
      /**
       * Calculate and tabulate power spectra
       */
      KFT::powerSpectra P (&C);
      P.initCorrelation (ps_table, cf_table);
      
      double a_final = 1.0;
      astro::functionWriter write ("powerSpectra.txt");
      write.push_back ([&] (double k) { return P.meanF (k, a_final); });
      write.push_back ([&] (double k) { return P.linearP (k, a_final); });
      write.push_back ([&] (double k) { return P.curlyP (k, a_final); });
      write.push_back ([&] (double k) { return P.BornApproxP (k, a_final); });
      write.add_header ("# Different types of cosmic density-fluctuation");
      write.add_header ("# power spectra as functions of wave number");
      write.add_header ("# column 1: wave number k in h/Mpc");
      write.add_header ("# column 2: mean-field interaction term");
      write.add_header ("# column 3: linearly evolved power spectrum");
      write.add_header ("# column 4: curly P power spectrum from KFT");
      write.add_header ("# column 5: mean-field non-linear KFT power spectrum");
      write (k_min, k_max, n, astro::LOG_SPACING);
      break;
    }
    default:
    {
      /**
       * Complain if you get here
       */
      throw std::invalid_argument ("unknown calculation type chosen");
      break;
    }
  }
  
  return 0;
}
