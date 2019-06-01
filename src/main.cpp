#include "../include/powerSpectraTeodori.h"
// Test program, the functions with comments are too slow right now
int main ()
{
    astro::cosmologyBase cos_model;
    Powermine spectrum (4, &cos_model);

    double a = 0.8;
    double k = 10;
    double q = 1;
    double mu = 0.5;

    // write functions
   /*  spectrum.writeIntegralmu(k); */
    // spectrum.writeQ(a);
    //spectrum.writeB1B2();
    /* spectrum.writeSfunctions(a); */
    spectrum.writeAll(a);
    std::cout << "Test 2d Levin:    "  << "\t" << spectrum.Integral2DLevin(k) << std::endl;
   /*  std::cout << "Test B_1 function:"  << "\t" << spectrum.B_1(q) << std::endl; */
    /* std::cout << "Test B_2 function:"  << "\t" << spectrum.B_2(q) << std::endl; */



    // double Qtest = spectrum.QFactor(a,k);
    // double p_initialTest = spectrum.p_initial(k);
    // double sigma_test = spectrum.sigma_1();
    // double B1Test = spectrum.B_1(q);
    // double B2Test = spectrum.B_2(q);
    // double muTest = spectrum.integral_mu(q,k);
    // // double qTest = spectrum.integral_q(k);
    // // double curlyPtest = spectrum.CurlyP(a,k) ;
    //
    // double qfirstTest = spectrum.integralq_first(k,mu);
    // // double musecondTest = spectrum.integralmu_second(k);
    //
    // double IntegralYTest = spectrum.integral_y(a,k) ;
    // double gradVTest  = spectrum.gradV(a,k);
    // double Integrand_SITest = spectrum.integrand_SI(a,k) ;
    // double SITest = spectrum.S_I(a,k) ;
    // std::cout << "Test Q factor:    "  << "\t" << Qtest << std::endl;
    // std::cout << "Test p initial:   " << "\t" << p_initialTest << std::endl;
    // std::cout << "Test sigma_1:     "  << "\t" << sigma_test << std::endl;
    // std::cout << "Test B_1 function:"  << "\t" << B1Test << std::endl;
    // std::cout << "Test B_2 function:"  << "\t" << B2Test << std::endl;
    // std::cout << "Test integral mu: "  << "\t" << muTest << std::endl;
    // // std::cout << "Test integral q:"  << "\t" << qTest << std::endl;
    // // std::cout << "Test curlyP :     "  << "\t" << curlyPtest  << std::endl;
    //
    // std::cout << "Test integral q first:"  << "\t" << qfirstTest << std::endl;
    // // std::cout << "Test integral mu sec:"  << "\t" << musecondTest << std::endl;
    // std::cout << "Test integral y:  "  << "\t" << IntegralYTest  << std::endl;
    /* std::cout << "Test gradient V:  "  << "\t" << gradVTest << std::endl; */
   /*  std::cout << "Test integrand S_I:"  << "\t" << Integrand_SITest  << std::endl; */
    /* std::cout << "Test S_I:          "  << "\t" << SITest << std::endl; */

	return 0;
}
