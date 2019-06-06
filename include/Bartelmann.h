#include <astro/cosmology/cosmologyBase.h>
#include <astro/cosmology/cosmicStructures.h>
#include <astro/cosmology/powerSpectrum.h>
#include <astro/cosmology/tophatFilter.h>
#include <astro/io/clArguments.h>
#include <astro/io/functionWriter.h>
#include <astro/utilities/functionTable.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <functional>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

#include <astro/utilities/osc_integrator.h>
#include <astro/utilities/integrator.h>
#include <astro/utilities/utilities.h>
#include <astro/utilities/LC_integrator.h>
#include <astro/utilities/LC_2d_integrator.h>

#include <KFT/particleDynamics/newZeldovichParticleDynamics.h>

#include <KFT/cosmology/kftCosmology.h>
#include <KFT/cosmology/powerSpectra.h>
#include <KFT/iniCorr/iniCorrTable.h>

/**
 * Specify your power spectrum here. The class must be derived from the
 * abstract powerSpectrum class of libastro for the libKFT to accept it.
 * It should suffice to modify the operator () function.
 */
class testPowerSpectrum: public astro::powerSpectrum
{
protected:
  const double k0 = 0.025;
  astro::cosmicStructures * cosmic_structures;
public:
  testPowerSpectrum
    (astro::cosmologyBase * cosmological_model_in, const double sigma8_in,
     astro::scaleFilter * filter_in);

  ~testPowerSpectrum ();

  double operator () (const double k, const double a);
  };

