#include <iostream>
#include <iomanip>
#include <fstream>
#include <functional>

#include <gsl/gsl_math.h>
#include <astro/utilities/osc_integrator.h>
#include <astro/cosmology/cosmologyBase.h>
#include <astro/cosmology/cosmicStructures.h>
#include <astro/io/functionWriter.h>
#include <KFT/particleDynamics/newZeldovichParticleDynamics.h>
#include <astro/utilities/utilities.h>


int main ()
{
    double k_0=1;
    double nbins = 64;
    double upperBound = 1000.0* k_0;

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


	return 0;
}
