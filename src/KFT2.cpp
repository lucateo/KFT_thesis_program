#include "../include/KFT2.h"

testPowerSpectrum::testPowerSpectrum(astro::cosmologyBase * cosmological_model_in,
        const double sigma8_in, astro::scaleFilter * filter_in):
  // The last one is exponent of large scale
  powerSpectrum (cosmological_model_in, sigma8_in, filter_in, 1.0)
  {
    cosmic_structures = new astro::cosmicStructures (cosmological_model);
    powerSpectrum::amplitude = sigma_8*sigma_8/sigma2 (8.0);
}

testPowerSpectrum::~testPowerSpectrum ()
{ delete cosmic_structures; }

double testPowerSpectrum::operator () (const double k, const double a)
{
double D_plus = cosmic_structures->Dplus (a);
double prefactor = D_plus*D_plus*amplitude;
double kappa = k/k0;
switch (m_initial_condition) {
        // Gaussian case, k here is actually k/gaussian_k0, so
        // gaussian_k0 sets the unit of measure
        case 0 : return prefactor * exp(-m_gaussNorm*gsl_pow_2(k -1));
        // Power case, small integer powers
        case 4 : return prefactor*k/gsl_pow_2(1 + gsl_pow_2(kappa));
        case 6 : return prefactor*k/gsl_pow_3(1 + gsl_pow_2(kappa));
        case 8 : return prefactor*k/gsl_pow_4(1 + gsl_pow_2(kappa));
        case 10 : return prefactor*k/gsl_pow_5(1 + gsl_pow_2(kappa));
        // Power case, half integer powers
        case 3 : return prefactor*k/gsl_pow_3( sqrt(1 + gsl_pow_2(kappa) ) ) ;
        case 5 : return prefactor*k/gsl_pow_5( sqrt(1 + gsl_pow_2(kappa) ) ) ;
        case 7 : return prefactor*k/gsl_pow_7( sqrt(1 + gsl_pow_2(kappa) ) ) ;
        case 9 : return prefactor*k/gsl_pow_9( sqrt(1 + gsl_pow_2(kappa) ) ) ;
        // non integer powers, i.e. consider number_in/10
        default : return prefactor*k/(exp((m_initial_condition/10)*log(1 + gsl_pow_2(kappa))));
    }
}

// trial function (to see if everything is ok)
void testPowerSpectrum::trial (double a)
{
    for (int i=0; i< n_bins; i++)
    {
        double k = astro::x_logarithmic (i, n_bins, k_min,k_max);
        std::cout << k << "\t" << operator () (k,a) << std::endl;
    }
}

// write functions
void testPowerSpectrum::writeCorrelations (KFT::kftCosmology * C )
{
    std::ostringstream os_ps, os_cf;
    os_ps << "data/ps_table_a_" << 1.0 << "_n_initial_"<< m_initial_condition << ".txt";
    os_cf << "data/cf_table_a_"<< 1.0 << "_n_initial_"<< m_initial_condition << ".txt";

    std::string ps_table = os_ps.str();
    std::string cf_table = os_cf.str();

    KFT::iniCorrTable corr_table
        (C->get_power_spectrum (), a_initial, q_min, q_max);

    corr_table.print_tables (ps_table, cf_table);
}

void testPowerSpectrum::testPowerSpectrum::writeSpectrum (KFT::kftCosmology * C, double a)
{
    std::ostringstream os, os_ps, os_cf ;
    os_ps << "data/ps_table_a_" << a << "_n_initial_"<< m_initial_condition << ".txt";
    os_cf << "data/cf_table_a_"<< a << "_n_initial_"<< m_initial_condition << ".txt";
    os << "data/powerSpectra_a_" << a << "_n_initial_"<< m_initial_condition << ".txt";
    std::string ps_table = os_ps.str();
    std::string cf_table = os_cf.str();
    std::string power_file = os.str();

    KFT::powerSpectra P (C);
    P.initCorrelation (ps_table, cf_table);

    astro::functionWriter write (power_file );
    write.push_back ([&] (double k) { return P.meanF (k, a); });
    write.push_back ([&] (double k) { return P.linearP (k, a); });
    write.push_back ([&] (double k) { return P.curlyP (k, a); });
    write.push_back ([&] (double k) { return P.BornApproxP (k, a); });
    write.add_header ("# Different types of cosmic density-fluctuation");
    write.add_header ("# power spectra as functions of wave number");
    write.add_header ("# column 1: wave number k in h/Mpc");
    write.add_header ("# column 2: mean-field interaction term");
    write.add_header ("# column 3: linearly evolved power spectrum");
    write.add_header ("# column 4: curly P power spectrum from KFT");
    write.add_header ("# column 5: mean-field non-linear KFT power spectrum");
    write (k_min, k_max, n_bins, astro::LOG_SPACING);
}

void testPowerSpectrum::writeAllSpectrum(KFT::kftCosmology * C, double a)
{
    std::ostringstream os, os_ps, os_cf ;
    os_ps << "data/ps_table_a_" << a << "_n_initial_"<< m_initial_condition << ".txt";
    os_cf << "data/cf_table_a_"<< a << "_n_initial_"<< m_initial_condition << ".txt";
    os << "data/powerSpectra_a_" << a << "_n_initial_"<< m_initial_condition << ".txt";
    std::string ps_table = os_ps.str();
    std::string cf_table = os_cf.str();
    std::string power_file = os.str();

    KFT::iniCorrTable corr_table
        (C->get_power_spectrum (), a_initial, q_min, q_max);
    corr_table.print_tables (ps_table, cf_table);

    KFT::powerSpectra P (C);
    P.initCorrelation (ps_table, cf_table);

    astro::functionWriter write (power_file );
    write.push_back ([&] (double k) { return P.meanF (k, a); });
    write.push_back ([&] (double k) { return P.linearP (k, a); });
    write.push_back ([&] (double k) { return P.curlyP (k, a); });
    write.push_back ([&] (double k) { return P.BornApproxP (k, a); });
    write.add_header ("# Different types of cosmic density-fluctuation");
    write.add_header ("# power spectra as functions of wave number");
    write.add_header ("# column 1: wave number k in h/Mpc");
    write.add_header ("# column 2: mean-field interaction term");
    write.add_header ("# column 3: linearly evolved power spectrum");
    write.add_header ("# column 4: curly P power spectrum from KFT");
    write.add_header ("# column 5: mean-field non-linear KFT power spectrum");
    write (k_min, k_max, n_bins, astro::LOG_SPACING);
}

// For the Gaussian initial conditions
void testPowerSpectrum::writeAllGaussian(KFT::kftCosmology * C, double a)
{
    std::ostringstream os, os_ps, os_cf ;
    os_ps << "data/ps_table_a_" << a << "_normGauss_"<< m_gaussNorm << ".txt";
    os_cf << "data/cf_table_a_"<< a << "_normGauss_"<< m_gaussNorm << ".txt";
    os << "data/GaussSpectra_a_" << a << "_normGauss_"<< m_gaussNorm << ".txt";
    std::string ps_table = os_ps.str();
    std::string cf_table = os_cf.str();
    std::string power_file = os.str();

    KFT::iniCorrTable corr_table
        (C->get_power_spectrum (), a_initial, q_min, q_max);
    corr_table.print_tables (ps_table, cf_table);

    KFT::powerSpectra P (C);
    P.initCorrelation (ps_table, cf_table);

    astro::functionWriter write (power_file );
    write.push_back ([&] (double k) { return P.meanF (k, a); });
    write.push_back ([&] (double k) { return P.linearP (k, a); });
    write.push_back ([&] (double k) { return P.curlyP (k, a); });
    write.push_back ([&] (double k) { return P.BornApproxP (k, a); });
    write.add_header ("# Different types of cosmic density-fluctuation");
    write.add_header ("# power spectra as functions of wave number");
    write.add_header ("# column 1: wave number k in units of k0");
    write.add_header ("# column 2: mean-field interaction term");
    write.add_header ("# column 3: linearly evolved power spectrum");
    write.add_header ("# column 4: curly P power spectrum from KFT");
    write.add_header ("# column 5: mean-field non-linear KFT power spectrum");
    write (k_min, k_max, n_bins, astro::LOG_SPACING);
}



