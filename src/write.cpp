#include "../include/powerSpectraTeodori.h"

// Plotting functions
void Powermine::writea1a2()
{
    int n = 1000;
    std::ofstream openFile;
    openFile.open("data/curlyP.txt");
    for (int i=0; i< n; i++)
    {
        double q = astro::x_logarithmic (i, n, 0.0001, 1000.0);
        openFile << std::setw(15) << q << "\t" << std::setw(15) << a_1(q)  << "\t" <<
            std::setw(15) << a_2(q) << std::endl;
    }
	openFile.close();
}

void Powermine::writeQ(double a)
{
    int n = 30;
    std::ostringstream os;
    os << "data/Q_a_" << a << "_n_initial_"<< m_number_in << "_A_" << A << "_k_0_" << k_0 << ".txt";

    std::ofstream openFile;
    openFile.open(os.str());
    for (int i=0; i< n; i++)
    {
        double k = astro::x_logarithmic (i, n, k_min, k_max);

        openFile << std::setw(15) << k << "\t" << std::setw(15)  << QFactor(a,k) << std::endl;
    }
	openFile.close();
}

void Powermine::writeSfunctions(double a)
{
    int n = 30;
    std::ostringstream os;
    os << "data/Sfunction_a_" << a << "_n_initial_"<< m_number_in << "_A_" << A << "_k_0_" << k_0 << ".txt";

    std::ofstream openFile;
    openFile.open(os.str());
    for (int i=0; i< n; i++)
    {
        double k = astro::x_logarithmic (i, n, k_min, k_max);

        openFile << std::setw(15) << k << "\t" << std::setw(15)  << integral_y(a,k) << "\t" <<
            std::setw(15)  << gradV(a,k) <<  "\t" << std::setw(15)  << integrand_SI(a,k)
            << "\t" << std::setw(15) << S_I(a,k)  << std::endl;
    }
	openFile.close();
}

void Powermine::writeAll(double a)
{
    int n = 30;
    std::ostringstream os;
    os << "data/All_a_" << a << "_n_initial_"<< m_number_in << "_A_" << A << "_k_0_" << k_0 << ".txt";

    std::ofstream openFile;
    openFile.open(os.str());

    // Line with the parameters
    openFile << std::setw(15) << a << "\t" << std::setw(15)  << m_number_in << "\t" <<
            std::setw(15)  << A <<  "\t" << std::setw(15)  << k_0 << std::endl;
    openFile << std::endl;

    for (int i=0; i< n; i++)
    {
        double k = astro::x_logarithmic (i, n, 0.001,1000.0);

        openFile << std::setw(15) << k << "\t" << std::setw(15)  << QFactor(a,k) << "\t" <<
            std::setw(15)  << integrand_SI(a,k) <<  /* "\t" << std::setw(15)  << Integral2DLevin(a,k) << */ std::endl;
    }
	openFile.close();
}

void Powermine::printTestSI(double a)
{
    std::cout << std::setw(15) << "k" << "\t" << std::setw(15)  << "integral y" << "\t" <<
            std::setw(15)  << "grad V" <<  "\t" << std::setw(15)  << "Integrand S_I"
            << "\t" << std::setw(15) << "S_I"  << std::endl;
    int n = 30;
    for (int i=0; i< n; i++)
    {
        double k = astro::x_logarithmic (i, n, 10.8, 100.0);

        std::cout << std::setw(15) << k << "\t" << std::setw(15)  << integral_y(a,k)  << "\t" << 
             std::setw(15)  << gradV(a,k) <<  "\t" << std::setw(15)  << integrand_SI(a,k)
             << "\t" << std::setw(15) << S_I(a,k)   << std::endl;
    }
}

void Powermine::printTesta1a2()
{
    std::cout << std::setw(15) << "q" << "\t" << std::setw(15)  << "B_1" << "\t" <<
            std::setw(15)  << "B_2"<<  std::endl;
    int n = 30;
    for (int i=0; i< n; i++)
    {
        double q = astro::x_logarithmic (i, n, 0.0001, 1000.0);

        std::cout << std::setw(15) << q << "\t" << std::setw(15)  << a_1(q)  << "\t" <<
             std::setw(15)  << a_2(q)  << std::endl;
    }
}

void Powermine::printTestintegralLevin(double a)
{
    std::cout << std::setw(15) << "k" << "\t" << std::setw(15)  << "integral" <<  std::endl;
    int n = 30;
    for (int i=0; i< n; i++)
    {
        double k = astro::x_logarithmic (i, n, 0.0001, 1000.0);

        std::cout << std::setw(15) << k << "\t" << std::setw(15)  << Integral2DLevin(a,k) << std::endl;
    }
}

void Powermine::printTestFullP(double a)
{
    std::cout << std::setw(15) << "k" << "\t" << std::setw(15)  << "Full P" <<  std::endl;
    int n = 30;
    for (int i=0; i< n; i++)
    {
        double k = astro::x_logarithmic (i, n, 0.001, 1000.0);

        std::cout << std::setw(15) << k << "\t" << std::setw(15)  << fullPowerSpectrum(a,k) << std::endl;
    }
}




// Dead end functions (keep it just for reference)

// double Powermine::integral_mu (double q, double k)
// {
//     std::function<double (double, double)> p = [k, this](double qk, double mu)
//         {return 2 * exp(mu * mu * B_1(qk/k)) - 1 ;};
//     astro::osc_integrator::Integrand p1 = p;
//     astro::osc_integrator Integral	(p1, astro::OSC_COSINE,
//         0,0,1, astro::LOG_SPACING, nbins );
//     return Integral(q);
// }
//
// double Powermine::integral_q (double k)
// {
//     std::function<double (double)> p = [k, this](double q)
//         {return q*q * exp(B_2(q)) * integral_mu(q,k) ;};
//     astro::integrator kernel (p);
//     return kernel(q_min, q_max);
// }
//
// double Powermine::integralq_first(double k, double mu)
// {
//     std::function<double (double, double)> p = [k, this](double kmu , double q)
//         {return q*q * exp(B_2(q)) *(exp( (kmu*kmu/(k*k))*B_1(q) ) - 1);};
//     astro::osc_integrator::Integrand p1 = p;
//     astro::osc_integrator Integral	(p1, astro::OSC_COSINE,
//         0,q_min,q_max, astro::LOG_SPACING, nbins );
//     return Integral(mu);
// }
//
// double Powermine::integralmu_second (double k)
// {
//     std::function<double (double)> p = [k, this](double mu)
//         {return 2*integralq_first(k,mu) ;};
//     astro::integrator kernel (p);
//     double mu_min = 0;
//     double mu_max = 1;
//     return kernel(mu_min, mu_max);
// }
//
//

