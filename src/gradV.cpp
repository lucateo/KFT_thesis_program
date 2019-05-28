#include <iostream>
#include <iomanip>
#include <fstream>

#include <astro/cosmology/cosmologyBase.h>
#include <astro/cosmology/cosmicStructures.h>
#include <astro/io/functionWriter.h>
#include <astro/utilities/utilities.h>

#include <KFT/particleDynamics/newZeldovichParticleDynamics.h>

const double a_min = 0.001;
astro::cosmologyBase cos_model;
//  cos_model.setDarkUniverse ();
astro::cosmicStructures cos_struct (&cos_model);
KFT::newZeldovichParticleDynamics prop (&cos_model, a_min);



double preFactorGradV(double a, double k)
{
    double D_plus = cos_struct.Dplus(a);
    double g = prop.g(a)    ;
    return 3*a*k*k*D_plus * D_plus /(8*M_PI*M_PI*g*g);
};

double integral_y (double k, double a)
{
    
    std::function<double (double)> p = [k,a](double y)
        {return q*q * exp(B_2(q)) * integral_mu(q,k) ;};
    astro::integrator kernel (p);
    return kernel(q_min, q_max);
};




int main ()
{
	double prefact = preFactorGradV(0.5,10);
    std::cout << prefact << std::endl;




	return 0;
}
