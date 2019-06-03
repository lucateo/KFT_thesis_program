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

// \sigma_1^2 in the main text
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

// Q_D in the main test
double Powermine::QFactor (double a, double k)
{
    double prop = propagator->g_qp(a);
    double Q = prop * prop * 0.33*k*k*sigma_1();
    return Q;
}

// B_1, B_2 are defined in chapter 4
double Powermine::B_1 (double q)
{
    // Trying to integrate step by step
    double sum= 0.0;
    std::vector<double> bins (nbins+1);
    for (int i = 0; i < nbins+1; i++)
        bins[i] = astro::x_logarithmic (i, nbins+1, k_min, k_max);

    // Integrand function (without j_2)
    std::function<double (double, double)> integrand = [this](double q, double k)
        {return (1/(2*M_PI*M_PI))* p_initial(k);};

    for (int i = 0; i < nbins; i++)
    {
        astro::LC_integrator integrate (integrand, astro::OSC_SPH_BESSEL, 2, bins[i], bins[i+1], 16);
        sum += integrate (q);
    }
    return sum;

    // astro::osc_integrator::Integrand p1 =  p;
    // astro::osc_integrator Integral	(p1, astro::OSC_SPH_BESSEL,
    //     2,k_min,k_max, astro::LOG_SPACING, nbins );
    //
    // return Integral(q);
};

double Powermine::B_2 (double q)
{
    double sum= 0.0;
    std::vector<double> bins (nbins+1);
    for (int i = 0; i < nbins+1; i++)
        bins[i] = astro::x_logarithmic (i, nbins+1, k_min, k_max);

    // Integrand function (without j_1)
    std::function<double (double, double)> integrand = [this](double q, double k)
        {return -(1/(2*M_PI*M_PI*k*q))* p_initial(k);};

    for (int i = 0; i < nbins; i++)
    {
        astro::LC_integrator integrate (integrand, astro::OSC_SPH_BESSEL, 1, bins[i], bins[i+1], 16);
        sum += integrate (q);
    }
    return sum;

    // astro::osc_integrator::Integrand p1 = p;
    // astro::osc_integrator Integral	(p1, astro::OSC_SPH_BESSEL,
    //     1,k_min,k_max, astro::LOG_SPACING, nbins );
    //
    // return Integral(q);
}

// The integral for \mathcal{P} without the prefactor (it is unstable!)
double Powermine::Integral2DLevin(double k)
{
    double sum= 0.0;
    std::vector<double> bins (nbins+1);
    for (int i = 0; i < nbins+1; i++)
        bins[i] = astro::x_logarithmic (i, nbins+1, q_min, q_max);

    // Integrand function (without cosine)
    std::function<double (double, double, double)>
    integrand = [this] (double k, double q, double mu)
    {    return q*q * exp(B_2(q)) *(exp(mu * mu * B_1(q)) - 1);};

    for (int i = 0; i < nbins; i++)
    {
        astro::LC_2d_integrator integrate
            (integrand, astro::OSC_COSINE, 0, -1.0, 1.0, bins[i], bins[i+1], 8, 8);
        sum += integrate (k);
    }
    return sum;
}

// Simply the \mathcal{P} with the prefactor
double Powermine::CurlyP(double a, double k)
{
    double g_qp = propagator->g_qp(a);
    double prefact = 2* M_PI *exp(-g_qp *g_qp*k*k);
    return prefact * Integral2DLevin(k);
}

// J(y) in the main text
double J (double y) {return 1 + (1 - y*y)/(2*y) *log((1 + y)/(1-y));}

// Just the integral in the y variable for gradient V
double Powermine::integral_y (double a, double k)
{
    double g_qp = propagator->g_qp(a);

    std::function<double (double)> p = [this,k,a,g_qp](double y)
        {return y*y * exp(-sigma_1()*g_qp*g_qp *k*k*y*y/3 ) * p_initial(k*y)*J(y) ;};

    astro::integrator kernel (p);
    return kernel(0.0, 1.0);
}

// averaged potential gradient
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
    double prefact =gdot * h/(g*g); // prefactor

    // integrand
    std::function<double (double)> p = [this,k](double a1)
        {return propagator->g(a1) *( propagator->jacobian(a1)/a1) * gradV(a1,k)  ;};

    astro::integrator kernel (p);
    return gradV(a,k) + prefact * kernel(0.0, a);

}

// the full S_I (without the -i factor)
double Powermine::S_I (double a, double k)
{
    std::function<double (double)> p = [this,a,k](double a1)
        {return propagator->g_qp(a,a1) *( propagator->jacobian(a1)/a1) * integrand_SI(a1,k)  ;};

    astro::integrator kernel (p);
    return 2*k*kernel(0,a);
}

// The full power spectrum
double Powermine::fullPowerSpectrum (double a, double k)
{
    double Q_D = QFactor(a,k);
    double curlyP = CurlyP(a,k);
    double SI = S_I(a,k);
    return exp(-Q_D + SI) * curlyP;
}



