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
	
        std::ofstream outfile;
        outfile.open("g_qp.txt"); //std::ios_base::app if you want to append it

        outfile  <<  "Computing g_qp(a1,a2), horizontal is a2, vertical is a1 " << std::endl;
	// number of steps and parameters
	const int n = 256;
	const double a_max = 1;
	const double k_min = 0.1;
	const double k_max = 100;
	
        // propagator g_qp (a,a')
        // put the first row with the values of a'
        for(int j =0; j<n; j++)
        {
            double a_temp2 = astro::x_logarithmic (j, n, a_min, a_max);
            outfile << std::setw(10) << a_temp2  << "\t";
        }
        outfile << std::endl;
        for (int i=0; i< n; i++)
        {
            double a_temp1 = astro::x_logarithmic (i, n, a_min, a_max);
            outfile << std::setw(10) << a_temp1 << "\t";
            for(int j =0; j<n; j++)
            {
                double a_temp2 = astro::x_logarithmic (j, n, a_min, a_max);
                double g_qptemp = prop.g_qp (a_temp1, a_temp2);
                outfile << std::setw(10) << g_qptemp  << "\t";
            }
            outfile << std::endl;
        }
	outfile.close();
	
	return 0;
}
