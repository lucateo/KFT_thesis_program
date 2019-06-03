#include "../include/powerSpectraTeodori.h"
// Test program 
int main ()
{
    // Initialize the power spectrum
    astro::cosmologyBase cos_model;
    Powermine spectrum (0, &cos_model);

    // Parameters used in trial functions
    double a = 0.01;
    double k = 10;
    double q = 1;
    double mu = 0.5;

    // Various test, uncomment the one you want
    //spectrum.printTestSI(a);
    //spectrum.printTestB1B2();
    spectrum.printTestintegralLevin(a);
    //spectrum.printTestFullP(a);

    // write functions (write in a txt file)
   /*  spectrum.writeIntegralmu(k); */
    // spectrum.writeQ(a);
    //spectrum.writeB1B2();
    /* spectrum.writeSfunctions(a); */
    // spectrum.writeAll(a);
   
    return 0;
}
