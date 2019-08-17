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
     astro::scaleFilter * filter_in, int initial_condition, double k0_gauss,
     double sigma_gauss);
  /// Destructor of the class
  ~testPowerSpectrum ();

  // variables
   /**
     * Sets the number for the initial condition, for m_initial_condition=0 it is the gaussian shape
     * \f[ P_{\delta}^i(k) = A \exp(-\frac{(k - k_0)^2}{2\sigma^2})  \f]
     * whereas for m_number_in \f$ \neq 0 \f$, it is related to
     * \f[ P_{\delta}^i(k) = Ak \qty( 1 + \qty(\frac{k}{k_0} )^2   )^{-n/2}  \f]
     * in particular, for m_number_in \f$ \le 100 \f$, the parameter is indeed the
     * \f$ n \f$ index, whereas if m_number_in > 100, then n = m_number_in/100
     * (see implementation of operator () for details)
     */
  int m_initial_condition = 4; 
  double m_sigma_gauss = 1; ///< determines sigma for gauss initial condition
  double m_gauss_k0 = 0.01; ///< determines  \f$ k_0 \f$ for gauss initial condition
  double m_gaussNorm; ///< determines \f$ k_0^2/(2 \sigma^2) \f$ factor in gaussian initial condition
  double a_initial = 1.0e-3; ///< deterines the initial time
  const double k_min = 1.0e-3, k_max = 1.0e3; //< determines plotting \f$ k \f$ range
  const double q_min = 1.0e-5, q_max = 1.0e5; //< determines \f$ q \f$ integration range
  const int n_bins = 128; ///< number of division in writing functions

  /// It overloads the initial, linearly evolved power spectrum (based on the initial
  /// number m_initial_condition
  double operator () (const double k, const double a);
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
  /// Power spectrum write in text with control at each step
  void writeAllSpectrumControl(KFT::kftCosmology * C, double a);
  /// The version for Gaussian initial condition of writeAllSpectrum
  void writeAllGaussian(KFT::kftCosmology * C, double a);
  /// For higher order curly P
  void writeAllHigherOrder(KFT::kftCosmology * C, double a, double k_prime, 
      double mu, int determine, int i_initial);
  /// Bispectrum with all the columns with various terms
  /// \param determine if determine = 0, computes correlation table and prints,
  /// if determine = 1 it does not compute correlation table and prints the header,
  /// if determine = 2 it does not compute correlation table and it does not print
  /// the header
  /// \param i_initial from where it starts the loop for plotting, if
  /// determine = 0 or 1, put this to 0
  void writeBiSpectrumFull(KFT::kftCosmology * C, double a, double ratio,
      int determine, int i_initial);
  /// Bispectrum with removed exp damping from single curly P factors
  /// \param determine if determine = 0, computes correlation table and prints,
  /// if determine = 1 it does not compute correlation table and prints the header,
  /// if determine = 2 it does not compute correlation table and it does not print
  /// the header
  /// \param i_initial from where it starts the loop for plotting, if
  /// determine = 0 or 1, put this to 0
  void writeBiSpectrumFixExpDamping(KFT::kftCosmology * C, double a,
    double ratio, int determine, int i_initial);
};

/**
 * \class powerSpectraModified
 * 
 * \brief Class providing power spectra extra tools 
 * 
 * class with extra tools for higher order correlators computations
 *
 * \author Luca Teodori
 */
class powerSpectraModified: public KFT::powerSpectra
{
  public:
  /// Constructor  
  powerSpectraModified(KFT::kftCosmology *cosmology_in);
  /// It computes the \f$ \mathcal{P}_{ij} \f$, with:
  /// \param mu_ij it is \f$ \mu_{ij} = \vec{L_{p_i}} \cdot \vec{L_{p_j}} \f$
  /// \param mu_i it is \f$ \mu_{i} = \vec{L_{p_i}} \cdot \vec{K_{ij}} \f$
  /// \param mu_j it is \f$ \mu_{ij} = \vec{L_{p_j}} \cdot \vec{K_{ij}} \f$
  double curlyP_ij (double a, double mu_ij, double mu_i, double mu_j, 
      double K_ij, double L_i, double L_j);

  void Trial(double a, double k);
  /// curly P_ij without exponential damping factor (it is not working!)
  double curlyP_ijNoDamping (double a, double mu_ij, double mu_i,
    double mu_j, double K_ij, double L_i, double L_j);
  /// To try to correct for bispectrum damping, to put in single \f$ \mathcal{P} \f$
  /// factors
  double DampingBispectrum(double a, double k_factor);
  /// try to compute \f$ \mathcal{P} \mathcal{P}\mathcal{P}\f$ factor
  double PPP_term (double k1, double a, double k2, double mu);
};
