#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <functional>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>


#include <astro/utilities/osc_integrator.h>
#include <astro/cosmology/cosmologyBase.h>
#include <astro/cosmology/cosmicStructures.h>
#include <astro/io/functionWriter.h>
#include <astro/utilities/integrator.h>
#include <astro/utilities/utilities.h>
#include <astro/utilities/LC_integrator.h>
#include <astro/utilities/LC_2d_integrator.h>

#include <KFT/particleDynamics/newZeldovichParticleDynamics.h>

/**
 * \ingroup PowerSpectrumEvolution
 *
 * \class Powermine
 *
 * \brief Class providing power spectra evolution tools used in my thesis
 *
 * Class providing all the functions and tools I used to compute the non linear
 * evolution of power spectra using the KFT approach. It is based on the òibastro
 * and libKFT.
 *
 * \author Luca Teodori
*/
class Powermine
{
    public:
    /// Pointer to the cosmological model of interest, as a default we will use
    /// a standard \f$ \Lambda \f$CDM model.
    astro::cosmologyBase * cosmological_model;
    /// Pointer to the library cosmicStructures, that provides utilities like
    /// growth factors and such; it is initialized with the cosmological_model;
    astro::cosmicStructures * cosmic_structures;
    /// New Zel'dovich propagator utilities.
    KFT::newZeldovichParticleDynamics * propagator;
    /**
     * Constructor for Power spectrum, it contains the initial power spectrum
     * and the cosmological model we are using.
     * \param number_in tells the initial condition, see p_initial;
     * \param cos_model is the cosmological model we are considering
    */
    Powermine (int number_in, astro::cosmologyBase * cos_model);
    ~Powermine(); ///< Destructor

    /**
     * Initial power spectrum, for m_number_in =0 it is the gaussian shape
     * \f[ P_{\delta}^i(k) = A \exp(-\frac{(k - k_0)^2}{2\sigma^2})  \f]
     * whereas for m_number_in \f$ \neq 0 \f$, it is related to
     * \f[ P_{\delta}^i(k) = Ak \qty( 1 + \qty(\frac{k}{k_0} )^2   )^{-n/2}  \f]
     * in particular, for half integer m_number_in \f$ \le 0 \f$, the parameter is indeed the
     * \f$ n \f$ index, whereas if m_number_in > 0, then \f$ n = m_number_in/10\f$
     * (see its implementation for details).
     */
    double p_initial (double k);
    double sigma_1(); ///< \f$ \sigma_1^2 \f$ factor
    double QFactor (double a, double k);///< The exponential damping \f$ Q_\mathrm{D} \f$ factor
    double B_1 (double q); ///< Computes the integral in \f$ j_2(q) \f$
    double B_2 (double q);  ///< Computes the integral in \f$ j_1(q) \f$
    double Integral2DLevin(double k);///< the 2D Levin integral for CurlyP
    double CurlyP(double a, double k); ///< The complete \f$ \bar{\mathcal{P}} \f$
        
    double integral_y (double a, double k); ///< Integral arising in computation of \f$ S_\mathrm{I} \f$
    double gradV (double a, double k); ///< The Born averaged gradient of the interacion potential
    double integrand_SI(double a, double k); ///< The integrand of \f$ S_\mathrm{I} \f$
    double S_I (double a, double k);///< The Born averaged interacting action
    double fullPowerSpectrum (double a, double k); ///< The complete non linearly evolved power spectrum

    // Plotting functions
    void writeB1B2(); ///< Function for writing data in text file, used then for plotting/analysis
    void writeQ(double a);///< Function for writing data in text file, used then for plotting/analysis
    void writeSfunctions(double a);///< Function for writing data in text file, used then for plotting/analysis
    void writeAll(double a);///< Function for writing data in text file, used then for plotting/analysis
    void printTestSI(double a);///< Function for print in terminal
    void printTestB1B2();///< Function for print in terminal
    void printTestintegralLevin(double a);///< Function for print in terminal
    void printTestFullP(double a);

    // Parameters
    double k_0=10;
    double A = 10; ///< Normalization of the power spectrum, to be determined by \f$ \sigma_8 \f$
    double nbins = 32; ///< Parameter for astro::osc_integrator utilities
    double k_max = k_0*(1.0e5); ///< The maximum \f$ k \f$ used in integration
    double k_min = 1.0e-5;///< The minimum \f$ k \f$ used in integration
    double q_min = 1.0e-5;///< The minimum \f$ q \f$ used in integration
    double q_max = 1.0e5;///< The maximum \f$ q \f$ used in integration
    double sigma_gauss = 1; ///< sigma for Gauss initial condition
    int m_number_in;///< It determines the initial condition
    double a_min = 0.001;///< The scale factor corresponding to the initial time considered
};


/**
 * J function, defined as
 * \f[ J(y) :=  \int_{-1}^{1} \dd{\mu} \frac{1 - \mu y}{1 + y^2 - 2\mu y} = 1 + \frac{1 - y^2}{2y}
 * \ln \frac{1+y}{\abs{1-y}}   \f]
 */
double J (double y);



