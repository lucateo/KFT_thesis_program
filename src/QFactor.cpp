#include <iostream>
#include <iomanip>
#include <fstream>

#include <astro/cosmology/cosmologyBase.h>
#include <astro/cosmology/cosmicStructures.h>
#include <astro/io/functionWriter.h>
#include <KFT/particleDynamics/newZeldovichParticleDynamics.h>
#include <astro/utilities/utilities.h>

//#include "../src/particleDynamics/newZeldovichParticleDynamics.h"

int main ()
{
	const double a_min = 0.001;
	astro::cosmologyBase cos_model;
	//  cos_model.setDarkUniverse ();
	astro::cosmicStructures cos_struct (&cos_model);
	KFT::newZeldovichParticleDynamics prop (&cos_model, a_min);

    	// number of steps and parameters
	const int n = 256;
	const double a_max = 1;
	const double k_min = 0.1;
	const double k_max = 100;

        //  writing Q stuff
       std::ofstream Qfile;
        Qfile.open("Q.txt");
        for (int i=0; i< n; i++)
        {
            double k = astro::x_logarithmic (i, n, k_min, k_max);
            for(int j =0; j<n; j++)
            {
                double a = astro::x_logarithmic (j, n, a_min, a_max);
                double propagator = prop.g_qp(a);
                double Q = propagator * propagator * 0.33*k*k;
                Qfile << std::setw(10) << Q  << "\t";
            }
            Qfile << std::endl;
        }
	Qfile.close();




	return 0;
}
