#include <iostream>
#include <sstream>
#include <fstream>

#include "astro/utilities/utilities.h"
#include "astro/io/clArguments.h"

#include "../src/cosmology/kftCosmology.h"
#include "../src/cosmology/interactionTerm.h"
#include "../src/cosmology/powerSpectra.h"

int main (int argc, char * argv[])
{
  // Simple class meant for easy handling of command-line arguments. 
  astro::clArguments ca (argc, argv);
  
  //get(const std::string  	opt, const int def) Returns the value of an 
  //option as an integer number, or the default def if the option was not set
  // -n option????
  const int n = ca.get ("-n", 128);
  const double k_min = 1.0e-3, k_max = 100.0;
  
  //Class bundling cosmological functions needed for KFT
  KFT::kftCosmology C;
 
  // Class bundling methods for calculating different power spectra. 
  // Constructor initializing the class with a pointer to a KFT cosmology 
  KFT::powerSpectra P (&C);
  
  //Initializes an interpolated correlation with the tables contained in the files named by the two strings given. 
  // ps = power spectrum, cf = correlation function
  P.initCorrelation ("ps_table.d", "cf_table.d");
  
  std::ofstream output ("powerSpectra.txt");
  output
    << "# Table of different types of cosmic density-fluctuation" << std::endl
    << "# power spectra as functions of wave number k in h/Mpc" << std::endl
    << "# column 1: wave number k in h/Mpc" << std::endl
    << "# column 2: mean interaction term in Born's approximation" << std::endl
    << "# column 3: linearly evolved power spectrum" << std::endl
    << "# column 4: numerically computed power spectrum" << std::endl
    << "# column 5: curly P power spectrum from KFT" << std::endl
    << "# column 6: non-linear KFT power spectrum in Born's approximation"
    << std::endl << std::endl;
  
  for (int i = 0; i < n; i++)
  {
    double k = astro::x_logarithmic (i, n, k_min, k_max);
    std::stringstream line;
    line << k
      << "  " << P.meanF (k, 1.0)
      << "  " << P.linearP (k, 1.0)
      << "  " << P.numericalP (k, 1.0)
      << "  " << P.curlyP (k, 1.0)
      << "  " << P.BornApproxP (k, 1.0);
    if (ca.set ("-v"))
      std::cout << line.str () << std::endl;
    output << line.str () << std::endl;
  }
  
  output.close ();
  
  return 0;
}
