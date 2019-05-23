#include <iostream>
#include <iomanip>
#include <fstream>
#include <functional>
#include <math.h>
//#include <cerf.h>
#ifdef __cplusplus
extern "C"
{
#endif
// C header here
#include <complex.h>

#ifdef __cplusplus
}
#endif

#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf.h>


#include <astro/utilities/osc_integrator.h>
#include <astro/cosmology/cosmologyBase.h>
#include <astro/cosmology/cosmicStructures.h>
#include <astro/io/functionWriter.h>
#include <KFT/particleDynamics/newZeldovichParticleDynamics.h>
#include <astro/utilities/utilities.h>

// GSL_POSINF infinity representation

// imaginary unit
gsl_complex iu = gsl_complex_rect(0, 1);

std::function<gsl_complex (double, double)> erfi_func = [](double q, double B_1)
{
    // negative imaginary unit
    gsl_complex iu_neg = gsl_complex_rect(0, -1);
    double prefactor_real =   M_SQRTPI /(2*sqrt(B_1)) * gsl_sf_exp(q*q/(4*B_1));
    double _Complex prova = I;
   // double _Complex erfi_part = cerfi(prova);
    //gsl_complex erfi_part =

    //gsl_complex prefactor = gsl_complex_mul (iu_neg , M_SQRTPI) ;
    return iu_neg;
}; // apparently you have to put ; before main function

    extern "C" {
      double _Complex cerf (double _Complex z);
    }

int main ()
{
    double k_0=1;
    double nbins = 64;
    double upperBound =  1000.0* k_0;//GSL_POSINF; it gives nan, seems to diverge!!!!!!

    // The part involving Bessel2
    // function < what returns (what has as input)> name = [](input){definition}
    std::function<double (double, double)> p_initial = [k_0](double q, double k)
        {return k* (1/gsl_pow_4(1 + gsl_pow_2(k/k_0)));};
    astro::osc_integrator::Integrand p1 = p_initial;
    astro::osc_integrator Integral1	(p1, astro::OSC_SPH_BESSEL,
        2,0.0,upperBound, astro::LOG_SPACING, nbins );
    std::cout << Integral1(10) << std::endl;

    // The part involving Bessel1
    std::function<double (double, double)> p_initial2 = [k_0](double q, double k)
        {return -(1/q)* (1/gsl_pow_4(1 + gsl_pow_2(k/k_0)));};

    astro::osc_integrator::Integrand p2 = p_initial2;
    astro::osc_integrator Integral2	(p2, astro::OSC_SPH_BESSEL,
        1,0.0,upperBound, astro::LOG_SPACING, nbins );
    std::cout << Integral2(10) << std::endl;
    //double _Complex prova(0,1);

    double _Complex erfi_part = cerf(5);
    //std::cout << prova << std::endl;
    std::cout << erfi_part << std::endl;


	return 0;
}
