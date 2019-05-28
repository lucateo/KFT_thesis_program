#include <iostream>
#include <iomanip>
#include <fstream>

#include <astro/cosmology/cosmologyBase.h>
#include <astro/cosmology/cosmicStructures.h>
#include <astro/io/functionWriter.h>
#include <astro/utilities/utilities.h>

#include <KFT/particleDynamics/newZeldovichParticleDynamics.h>

#include "../include/powerSpectraTeodori.h"



Powermine::Powermine(int number_in,  astro::cosmologyBase cos_model ) 
{
    m_number_in=number_in;
        // Cosmology setting
    astro::cosmicStructures cos_struct (&cos_model); 
    KFT::newZeldovichParticleDynamics prop (&cos_model, a_min); 

}
double Powermine::p_initial (double k)
{
     switch (m_number_in) {
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
        default : return A*k*(1/exp((m_number_in/10)*log(1 + gsl_pow_2(k/k_0))));
    }
}
double Powermine::sigma_1 ()
{
    if (m_number_in == 0)
        // Gaussian case
        return M_SQRT2 * M_SQRTPI * sigma_gauss;
    else if (m_number_in <=10)
        // Power case
        return (A*k_0*k_0)/(2*M_PI * M_PI *(m_number_in - 2));
    else // Non integer power
         return (A*k_0*k_0)/(2*M_PI * M_PI *((m_number_in/10) - 2));
}

// Then for the real Q factor I have to multiply this by sigma_1
double QFactor (Powermine spectrum, double a, double k)
{

    double propagator = spectrum.prop.g_qp(a);
    double Q = propagator * propagator * 0.33*k*k*spectrum.sigma_1();
    return Q;
}


double B_1 (Powermine spectrum, double q)
{
    std::function<double (double, double)> p= [spectrum](double q, double k)
        {return (1/(2*M_PI*M_PI))* spectrum.p_initial(k);};
    astro::osc_integrator::Integrand p1 =  p;
    astro::osc_integrator Integral	(p1, astro::OSC_SPH_BESSEL,
        2,spectrum.k_min,spectrum.k_max, astro::LOG_SPACING, spectrum.nbins );
    return Integral(q);
};

double B_2 (Powermine spectrum, double q)
{
    std::function<double (double, double)> p = [spectrum](double q, double k)
        {return -(1/(2*M_PI*M_PI*k*q))* spectrum.p_initial(k);};
    astro::osc_integrator::Integrand p1 = p;
    astro::osc_integrator Integral	(p1, astro::OSC_SPH_BESSEL,
        1,spectrum.k_min,spectrum.k_max, astro::LOG_SPACING, spectrum.nbins );
    return Integral(q);
}

double integral_mu (Powermine spectrum, double q, double k)
{
    std::function<double (double, double)> p = [spectrum, q](double qk, double mu)
        {return 2 * exp(mu * mu * B_1(spectrum, q)) - 1 ;};
    astro::osc_integrator::Integrand p1 = p;
    astro::osc_integrator Integral	(p1, astro::OSC_COSINE,
        0,0,1, astro::LOG_SPACING, spectrum.nbins );
    return Integral(q);
}

double integral_q (Powermine spectrum, double k)
{
    std::function<double (double)> p = [spectrum, k](double q)
        {return q*q * exp(B_2(spectrum,q)) * integral_mu(spectrum,q,k) ;};
    astro::integrator kernel (p);
    return kernel(spectrum.q_min, spectrum.q_max);
}


double preFactorCurlyP(Powermine spectrum, double k, double a)
{
    double g_qp = spectrum.prop.g_qp(a);
    return 2* M_PI *exp(-g_qp *g_qp*k*k);
}

double preFactorGradV(Powermine spectrum, double a, double k)
{
    double D_plus = spectrum.cos_struct.Dplus(a);
    double g = spectrum.prop.g(a)    ;
    return 3*a*k*k*D_plus * D_plus /(8*M_PI*M_PI*g*g);
}


// Plotting functions

void writeCurlyP1(Powermine spectrum)
{
    int n = 1000;
    std::ofstream openFile;
    openFile.open("curlyP.txt");
    for (int i=0; i< n; i++)
    {
        double q = astro::x_logarithmic (i, n, spectrum.q_min, spectrum.q_max);
        openFile << std::setw(15) << q << "\t" << std::setw(15) << B_1(spectrum,q)  << "\t" <<
            std::setw(15) << B_2(spectrum, q) << std::endl;
    }
	openFile.close();
}

void writeCurlyP2(Powermine spectrum, double k)
{
    int n = 100;
    std::ofstream openFile;
    openFile.open("curlyP2.txt");
    for (int i=0; i< n; i++)
    {
        double q = astro::x_logarithmic (i, n, spectrum.q_min, spectrum.q_max);

        openFile << std::setw(15) << q << "\t" << std::setw(15)/*  << k  << "\t" <<
            std::setw(15) */ << integral_mu(spectrum,q,k) << std::endl;
    }
	openFile.close();
}



