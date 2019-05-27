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
    Powermine(int number_in);

    double p_initial (double k);



    private:
    // Parameters
    const double k_0=10;
    const double A = 100; // normalization
    const double nbins = 32; // seems good also with 32
    const double k_max = 1000.0;
    const double k_min = 0.001;
    const double q_min = 0.001;
    const double q_max = 10000.0;
    const double sigma_gauss = 1; // sigma for Gauss initial condition

    // Cosmology setting
    const double a_min = 0.001;
    astro::cosmologyBase cos_model;
    //  cos_model.setDarkUniverse ();
    astro::cosmicStructures cos_struct (&cos_model);
    KFT::newZeldovichParticleDynamics propagator (&cos_model, a_min);




};
