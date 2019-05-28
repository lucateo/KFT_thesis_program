#include <iostream>
#include <iomanip>
#include <fstream>

#include <astro/cosmology/cosmologyBase.h>
#include <astro/cosmology/cosmicStructures.h>
#include <astro/io/functionWriter.h>
#include <astro/utilities/utilities.h>

#include <KFT/particleDynamics/newZeldovichParticleDynamics.h>

#include "../include/powerSpectraTeodori.h"

int main ()
{
    astro::cosmologyBase cos_model;
    Powermine spectrum (4, cos_model);
    double prefact = preFactorGradV(spectrum, 0.5,10);
    std::cout << prefact << std::endl;




	return 0;
}
