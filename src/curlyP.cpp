#include <iostream>
#include <iomanip>
#include <fstream>
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

#include <KFT/particleDynamics/newZeldovichParticleDynamics.h>
// GSL_POSINF infinity representation
// TODO: mu integral and then q integral
// ISSUES: LINEAR SPACE seems to give nan,
// CHECK: routines indeed compute what they should
// POSSIBLE THINGS TO DO: put everything in classes (and Doxygen documentation)

// Parameters
const double k_0=10;
const double A = 100; // normalization
const double nbins = 32; // seems good also with 32
const double k_max = 10000.0;
const double k_min = 0.001;
const double q_min = 0.001;
const double q_max = 10000.0;
const double sigma_gauss = 1; // sigma for Gauss initial condition

// Cosmology setting
const double a_min = 0.001;
astro::cosmologyBase cos_model;
//  cos_model.setDarkUniverse ();
astro::cosmicStructures cos_struct (&cos_model);
KFT::newZeldovichParticleDynamics prop (&cos_model, a_min);


// Initial power spectrum choice, if it is greater than 10, it means that
// I am choosing number_in/10
const int number_in = 4 ;

double sigma_1; // factor \sigma_1^2
// I think you have to set it in main function
void determine_sigma()
{
    if (number_in == 0)
        // Gaussian case
        sigma_1 = M_SQRT2 * M_SQRTPI * sigma_gauss;
    else if (number_in <=10)
        // Power case
        sigma_1 = (A*k_0*k_0)/(2*M_PI * M_PI *(number_in - 2));
    else // Non integer power
         sigma_1 = (A*k_0*k_0)/(2*M_PI * M_PI *((number_in/10) - 2));
};


// Initial power spectrum
double p_initial (double k)
{
    switch (number_in) {
        // Gaussian case
        case 0 : return A * exp(-gsl_pow_2(k - k_0)/(2*sigma_gauss)) ;
        // Power case, small integer powers
        case 4 : return A*k*(1/gsl_pow_2(1 + gsl_pow_2(k/k_0)));
        case 6 : return A*k*(1/gsl_pow_3(1 + gsl_pow_2(k/k_0)));
        case 8 : return A*k*(1/gsl_pow_4(1 + gsl_pow_2(k/k_0)));
        case 10 : return A*k*(1/gsl_pow_5(1 + gsl_pow_2(k/k_0)));
        // Power case, half integer powers
        case 3 : return A*k*(1/gsl_pow_3( sqrt(1 + gsl_pow_2(k/k_0) ) ) );
        case 5 : return A*k*(1/gsl_pow_5( sqrt(1 + gsl_pow_2(k/k_0) ) ) );
        case 7 : return A*k*(1/gsl_pow_7( sqrt(1 + gsl_pow_2(k/k_0) ) ) );
        case 9 : return A*k*(1/gsl_pow_9( sqrt(1 + gsl_pow_2(k/k_0) ) ) );
        // non integer powers, i.e. consider number_in/10
        default : return A*k*(1/exp((number_in/10)*log(1 + gsl_pow_2(k/k_0))));
    }
};

double B_1 (double q)
{
    std::function<double (double, double)> p= [](double q, double k)
        {return (1/(2*M_PI*M_PI))* p_initial(k);};
    astro::osc_integrator::Integrand p1 =  p;
    astro::osc_integrator Integral	(p1, astro::OSC_SPH_BESSEL,
        2,k_min,k_max, astro::LOG_SPACING, nbins );
    return Integral(q);
};

double B_2 (double q)
{
    std::function<double (double, double)> p = [](double q, double k)
        {return -(1/(2*M_PI*M_PI*k*q))* p_initial(k);};
    astro::osc_integrator::Integrand p1 = p;
    astro::osc_integrator Integral	(p1, astro::OSC_SPH_BESSEL,
        1,k_min,k_max, astro::LOG_SPACING, nbins );
    return Integral(q);
};

double integral_mu (double q, double k)
{
    std::function<double (double, double)> p = [q](double qk, double mu)
        {return 2 * exp(mu * mu * B_1(q)) - 1 ;};
    astro::osc_integrator::Integrand p1 = p;
    astro::osc_integrator Integral	(p1, astro::OSC_COSINE,
        0,0,1, astro::LOG_SPACING, nbins );
    return Integral(q);
};

double integral_q (double k)
{
    std::function<double (double)> p = [k](double q)
        {return q*q * exp(B_2(q)) * integral_mu(q,k) ;};
    astro::integrator kernel (p);
    return kernel(q_min, q_max);
};


double preFactorCurlyP(double k, double a)
{
    double g_qp = prop.g_qp(a);
    return 2* M_PI *exp(-g_qp *g_qp*k*k);
};



// Plotting functions

void writeCurlyP1()
{
    int n = 1000;
    std::ofstream openFile;
    openFile.open("curlyP.txt");
    for (int i=0; i< n; i++)
    {
        double q = astro::x_logarithmic (i, n, q_min, q_max);
        openFile << std::setw(15) << q << "\t" << std::setw(15) << B_1(q)  << "\t" <<
            std::setw(15) << B_2(q) << std::endl;
    }
	openFile.close();
};

void writeCurlyP2(double k)
{
    int n = 100;
    std::ofstream openFile;
    openFile.open("curlyP2.txt");
    for (int i=0; i< n; i++)
    {
        double q = astro::x_logarithmic (i, n, q_min, q_max);

        openFile << std::setw(15) << q << "\t" << std::setw(15)/*  << k  << "\t" <<
            std::setw(15) */ << integral_mu(q,k) << std::endl;
    }
	openFile.close();
};


int main ()
{
    determine_sigma();
    std::cout << preFactorCurlyP(10,0.1) << std::endl;






    return 0;
}
