#include <iostream>
#include <iomanip>
#include <fstream>

#include <astro/cosmology/cosmologyBase.h>
#include <astro/cosmology/cosmicStructures.h>
#include <astro/io/functionWriter.h>
#include <astro/utilities/utilities.h>

#include <KFT/particleDynamics/newZeldovichParticleDynamics.h>



class Powermine
{
    public:
    Powermine ();    
    Powermine (int number_in, astro::cosmologyBase cos_model)
    {
        m_number_in=number_in;
        // Cosmology setting
        astro::cosmicStructures cos_struct (&cos_model); 
        KFT::newZeldovichParticleDynamics prop (&cos_model, a_min); 
    };
    ~Powermine();

    /** Initial power spectrum */
    double p_initial (double k);
    double sigma_1();

   
    // Parameters
    const double k_0=10;
    const double A = 100; // normalization
    const double nbins = 32; // seems good also with 32
    const double k_max = 1000.0;
    const double k_min = 0.001;
    const double q_min = 0.001;
    const double q_max = 10000.0;
    const double sigma_gauss = 1; // sigma for Gauss initial condition
    int m_number_in = 4;

    // Cosmology setting
    const double a_min = 0.001;
    //  cos_model.setDarkUniverse ();
    astro::cosmicStructures cos_struct (&cos_model);
    KFT::newZeldovichParticleDynamics prop (&cos_model, a_min);

};

double QFactor (Powermine spectrum, double a, double k);

double B_1 (Powermine spectrum, double q); ///< Computes the integral in \f$ j_2(q) \f$
double B_2 (Powermine spectrum, double q);  ///< Computes the integral in \f$ j_1(q) \f$
double integral_mu (Powermine spectrum, double q, double k); ///< Computes the integral in \f$ \mu \f$
double integral_q (Powermine spectrum, double k);  ///< Computes the integral in \f$ q \f$


/// Prefactor of \f$ \mathcal{P} \f$
double preFactorCurlyP(Powermine spectrum, double k, double a);
/// Prefactor of \f$ \grad{V} \f$
double preFactorGradV(Powermine spectrum, double a, double k);

/// Function for writing in a text file q, B_1 and B_2
void writeCurlyP1(Powermine spectrum );
/// Function for writing in a text file q, integral in \f$ \mu \f$
void writeCurlyP2(Powermine spectrum, double k);






