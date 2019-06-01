#include "../include/powerSpectraTeodori.h"


Powermine::Powermine(int number_in,  astro::cosmologyBase * cos_model ) :
    m_number_in (number_in), cosmological_model (cos_model)
{
        // Cosmology setting
    cosmic_structures = new astro::cosmicStructures (cosmological_model);
    propagator = new KFT::newZeldovichParticleDynamics (cosmological_model, a_min);

}

Powermine::~Powermine()
{
    delete cosmic_structures;
    delete propagator;
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

double Powermine::QFactor (double a, double k)
{
    double prop = propagator->g_qp(a);
    double Q = prop * prop * 0.33*k*k*sigma_1();
    return Q;
}


double Powermine::B_1 (double q)
{
    std::function<double (double, double)> p= [this](double q, double k)
        {return (1/(2*M_PI*M_PI))* p_initial(k);};
    astro::osc_integrator::Integrand p1 =  p;
    astro::osc_integrator Integral	(p1, astro::OSC_SPH_BESSEL,
        2,k_min,k_max, astro::LOG_SPACING, nbins );
    return Integral(q);
};

double Powermine::B_2 (double q)
{
    std::function<double (double, double)> p = [this](double q, double k)
        {return -(1/(2*M_PI*M_PI*k*q))* p_initial(k);};
    astro::osc_integrator::Integrand p1 = p;
    astro::osc_integrator Integral	(p1, astro::OSC_SPH_BESSEL,
        1,k_min,k_max, astro::LOG_SPACING, nbins );
    return Integral(q);
}

double Powermine::Integral2DLevin(double k)
{
    std::function<double (double, double, double)>
    kernel = [this] (double k, double q, double mu)
    {    return q*q * exp(B_2(q)) *(exp(mu * mu * B_1(q)) - 1);};
    astro::LC_2d_integrator
    Test (kernel, astro::OSC_COSINE, 0, -1.0, 1.0, q_min, q_max, 8, 8);
    return Test (k);
}

double Powermine::CurlyP(double a, double k)
{
    double g_qp = propagator->g_qp(a);
    double prefact = 2* M_PI *exp(-g_qp *g_qp*k*k);
    return prefact * Integral2DLevin(k);
}

// Plotting functions
void Powermine::writeB1B2()
{
    int n = 1000;
    std::ofstream openFile;
    openFile.open("data/curlyP.txt");
    for (int i=0; i< n; i++)
    {
        double q = astro::x_logarithmic (i, n, q_min, q_max);
        openFile << std::setw(15) << q << "\t" << std::setw(15) << B_1(q)  << "\t" <<
            std::setw(15) << B_2(q) << std::endl;
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
            std::setw(15)  << integrand_SI(a,k) <<  /* "\t" << std::setw(15)  << CurlyP(a,k) << */ std::endl;
    }
	openFile.close();
}







// GradV functions
double J (double y) {return 1 + (1 - y*y)/(2*y) *log((1 + y)/(1-y));}

double Powermine::integral_y (double a, double k)
{
    double g_qp = propagator->g_qp(a);
    std::function<double (double)> p = [this,k,a,g_qp](double y)
        {return y*y * exp(-sigma_1()*g_qp*g_qp *k*k*y*y/3 ) * p_initial(k*y)*J(y) ;};
    astro::integrator kernel (p);
    return kernel(0.001, 0.999);
}

double Powermine::gradV (double a, double k)
{
    double D_plus = cosmic_structures->Dplus(a);
    double g = propagator->g(a);
    double prefact = 3*a*k*k*D_plus * D_plus /(8*M_PI*M_PI*g*g);
    return prefact*integral_y(a,k);
}


// S_I functions

double Powermine::integrand_SI(double a, double k)
{
    double gdot = propagator->g_dot(a);
    double h =  propagator->h(a);
    double g = propagator->g(a);
    double prefact =gdot * h/(g*g);
    std::function<double (double)> p = [this,k](double a1)
        {return propagator->g(a1) *( propagator->jacobian(a1)/a1) * gradV(a1,k)  ;};
    astro::integrator kernel (p);
    return gradV(a,k) + prefact * kernel(0.0, a);

}

double Powermine::S_I (double a, double k)
{
    std::function<double (double)> p = [this,a,k](double a1)
        {return propagator->g_qp(a,a1) *( propagator->jacobian(a1)/a1) * integrand_SI(a1,k)  ;};
    astro::integrator kernel (p);
    return 2*k*kernel(0,a);
}

double Powermine::fullPowerSpectrum (double a, double k)
{
    double Q_D = QFactor(a,k);
    double curlyP = CurlyP(a,k);
    double SI = S_I(a,k);
    return exp(-Q_D + SI) * curlyP;
}




// Dead functions (keep it just for reference)

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
/*  */

