#include <astro/cosmology/cosmologyBase.h>
#include <astro/cosmology/cosmicStructures.h>
#include <astro/cosmology/powerSpectrum.h>
#include <astro/cosmology/tophatFilter.h>
#include <astro/io/clArguments.h>
#include <astro/io/functionWriter.h>
#include <astro/utilities/functionTable.h>
#include <astro/utilities/osc_integrator.h>
#include <astro/utilities/integrator.h>
#include <astro/utilities/utilities.h>
#include <astro/utilities/LC_integrator.h>
#include <astro/utilities/LC_2d_integrator.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <functional>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>


#include <KFT/particleDynamics/newZeldovichParticleDynamics.h>
#include <KFT/cosmology/kftCosmology.h>
#include <KFT/cosmology/powerSpectra.h>
#include <KFT/iniCorr/iniCorrTable.h>

/**
 * \class testPowerSpectrum
 * 
 * \brief Class providing power spectra evolution tools used in my thesis 
 *
 * Specify your power spectrum here. The class must be derived from the
 * abstract powerSpectrum class of libastro for the libKFT to accept it.
 * It should suffice to modify the operator () function.
 *
 * \author Matthias Bartelmann, modified by Luca Teodori
 */
class testPowerSpectrum: public astro::powerSpectrum
{
  protected:
  /// Radiation-matter equality scale  
  const double k0 = 0.025;
  astro::cosmicStructures * cosmic_structures;
  
  public:
  /** 
   * The constructor of the class, it initialize the power spectrum with the
   * cosmological model
   * \param sigma8_in the normalization of the power spectrum (usually it is 0.8)
   * \param filter_in The filter to compute \f$ \sigma_8 \f$
   * \param cosmological_model_in pointer to a cosmology model (usually a a standard \f$ \Lambda \f$CDM one
   */
  testPowerSpectrum
    (astro::cosmologyBase * cosmological_model_in, const double sigma8_in,
     astro::scaleFilter * filter_in );

  ~testPowerSpectrum ();

  // variables
  int m_initial_condition = 0; ///< determines the initial condition to use
  double m_sigma_gauss = 0.001; ///< determines sigma for gauss initial condition
  double m_gauss_k0 = 0.001; ///< determines  \f$ k_0 \f$ for gauss initial condition
  double m_gaussNorm = 1; ///< determines \f$ k_0^2/(2 \sigma^2) \f$ factor in gaussian initial condition
  double a_initial = 1.0e-3; ///< deterines the initial time
  const double k_min = 1.0e-3, k_max = 100.0;
  const double q_min = 1.0e-5, q_max = 1.0e5;
  const int n_bins = 128; ///< number of division in writing functions

  /// It overloads the initial, linearly evolved power spectrum (based on the initial
  /// number m_initial_condition
  double operator () (const double k, const double a);
   /**
     * Sets the number for the initial condition, for m_number_in =0 it is the gaussian shape
     * \f[ P_{\delta}^i(k) = A \exp(-\frac{(k - k_0)^2}{2\sigma^2})  \f]
     * whereas for m_number_in \f$ \neq 0 \f$, it is related to
     * \f[ P_{\delta}^i(k) = Ak \qty( 1 + \qty(\frac{k}{k_0} )^2   )^{-n/2}  \f]
     * in particular, for half integer m_number_in \f$ \le 10 \f$, the parameter is indeed the
     * \f$ n \f$ index, whereas if m_number_in > 0, then n = m_number_in/10
     * (see implementation of operator () for details)
     */
  void setInitialCondition(int n) {m_initial_condition = n;};
  /// It sets sigma_gauss in gaussian initial condition case to the power spectrum
  void setSigmaGauss(double sigma_gauss) {m_sigma_gauss = sigma_gauss;};
  /// It sets k0 in gaussian initial condition case to the power spectrum
  void setGauss_k0(double k0_gauss) {m_gauss_k0 = k0_gauss;};
  /// It sets in gaussian initial condition case to the power spectrum
  void setGaussNorm(double gaussNorm) {m_gaussNorm = gaussNorm;};
  /// Trial function, prints the initial power spectrum
  void trial (double a);
  /// Trial function, it writes to file the initial power spectrum
  void writeTest (double a);
  /// write in file correlation tables
  void writeCorrelations(KFT::kftCosmology * C);
  /// write in file the various quantities for the spectra (it needs 
  /// writeCorrelations to be run before for the relevant initial condition)
  void writeSpectrum (KFT::kftCosmology * C, double a);
  /// write both correlation tables and text file with quantities relevant for
  /// power spectrum
  void writeAllSpectrum(KFT::kftCosmology * C, double a);
  /// The version for Gaussian initial condition of writeAllSpectrum
  void writeAllGaussian(KFT::kftCosmology * C, double a);
};

