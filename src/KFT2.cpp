#include "../include/KFT2.h"

testPowerSpectrum::testPowerSpectrum(astro::cosmologyBase * cosmological_model_in,
  const double sigma8_in, astro::scaleFilter * filter_in, int initial_condition,
  double k0_gauss, double sigma_gauss):
  // The last one is exponent of large scale
  powerSpectrum (cosmological_model_in, sigma8_in, filter_in, 1.0)
  {
    m_initial_condition = initial_condition;
    m_gauss_k0 = k0_gauss;
    m_sigma_gauss = sigma_gauss;
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
  // Gaussian case
  if (m_initial_condition == 0 )
    return prefactor * exp(-gsl_pow_2(k -m_gauss_k0)/(2*m_sigma_gauss));
  // even n case
  else if (m_initial_condition == 4 )
    return prefactor*k/gsl_pow_2(1.0 + gsl_pow_2(kappa));
  else if (m_initial_condition == 6 )
    return prefactor*k/gsl_pow_3(1.0 + gsl_pow_2(kappa));
  else if (m_initial_condition == 8 )
    return prefactor*k/gsl_pow_4(1.0 + gsl_pow_2(kappa));
  else if (m_initial_condition == 10 )
    return prefactor*k/gsl_pow_5(1.0 + gsl_pow_2(kappa));
  else if (m_initial_condition == 12 )
    return prefactor*k/gsl_pow_6(1.0 + gsl_pow_2(kappa));
  else if (m_initial_condition == 14 )
    return prefactor*k/gsl_pow_7(1.0 + gsl_pow_2(kappa));
  else if (m_initial_condition == 16 )
    return prefactor*k/gsl_pow_8(1.0 + gsl_pow_2(kappa));
  else if (m_initial_condition == 18 )
    return prefactor*k/gsl_pow_9(1.0 + gsl_pow_2(kappa));
  // odd n case
  else if (m_initial_condition == 3 )
    return prefactor*k/gsl_pow_3(sqrt(1.0 + gsl_pow_2(kappa)));
  else if (m_initial_condition == 5 )
    return prefactor*k/gsl_pow_5(sqrt(1.0 + gsl_pow_2(kappa)));
  else if (m_initial_condition == 7 )
    return prefactor*k/gsl_pow_7(sqrt(1.0 + gsl_pow_2(kappa)));
  else if (m_initial_condition == 9 )
    return prefactor*k/gsl_pow_9(sqrt(1.0 + gsl_pow_2(kappa)));
  else if (m_initial_condition == 11 )
    return prefactor*k/gsl_pow_int(sqrt(1.0 + gsl_pow_2(kappa)), 11);
  else if (m_initial_condition == 13 )
    return prefactor*k/gsl_pow_int(sqrt(1.0 + gsl_pow_2(kappa)), 13);
  else if (m_initial_condition == 15 )
    return prefactor*k/gsl_pow_int(sqrt(1.0 + gsl_pow_2(kappa)), 15);
  else if (m_initial_condition == 17 )
    return prefactor*k/gsl_pow_int(sqrt(1.0 + gsl_pow_2(kappa)), 17);
  else if (m_initial_condition == 19 )
    return prefactor*k/gsl_pow_int(sqrt(1.0 + gsl_pow_2(kappa)), 19);
  // other cases
  else if (m_initial_condition < 100)
    return prefactor*k/(exp(double(double(m_initial_condition)/2.0)
             *log(1.0 + gsl_pow_2(kappa))));
  else
    return prefactor*k/(exp(double(double(m_initial_condition)/200.0)
             *log(1.0 + gsl_pow_2(kappa))));
}

// trial function (to see if everything is ok)
void testPowerSpectrum::trial (double a)
{
  double D_plus = cosmic_structures->Dplus (a);
  std::cout << D_plus << "\t" << amplitude << std::endl;
  for (int i=0; i< n_bins; i++)
  {
    double k = astro::x_logarithmic (i, n_bins, k_min,k_max);
    std::cout << k << "\t" << operator () (k,a) << std::endl;
  }
}

void testPowerSpectrum::writeTest (double a)
{
  std::ofstream openFile;
  openFile.open("data/Test.txt");
  for (int i=0; i< n_bins; i++)
  {
    double k = astro::x_logarithmic (double(i), double(n_bins), k_min, k_max);

    openFile << std::setw(15) << k << "\t" << std::setw(15)
      << operator () (k, a) << std::endl;
  }
	openFile.close();
}

// write functions
void testPowerSpectrum::writeCorrelations (KFT::kftCosmology * C )
{
  std::ostringstream os_ps, os_cf;
  os_ps << "data/ps_table/ps_table_a_" << 1.0 << "_n_initial_"
    << m_initial_condition << ".d";
  os_cf << "data/cf_table/cf_table_a_"<< 1.0 << "_n_initial_"
    << m_initial_condition << ".d";

  std::string ps_table = os_ps.str();
  std::string cf_table = os_cf.str();

  KFT::iniCorrTable corr_table
    (C->get_power_spectrum (), a_initial, q_min, q_max);

  corr_table.print_tables (ps_table, cf_table);
}

void testPowerSpectrum::writeSpectrum (KFT::kftCosmology * C, double a)
{
  std::ostringstream os, os_ps, os_cf ;
  os_ps << "data/ps_table/ps_table_a_" << a << "_n_initial_"
    << m_initial_condition << ".d";

  os_cf << "data/cf_table/cf_table_a_"<< a << "_n_initial_"
    << m_initial_condition << ".d";

  os << "data/powerSpectra_a_" << a << "_n_initial_"<< m_initial_condition
    << ".txt";

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
  if (m_initial_condition < 100)
  {
    os_ps << "data/ps_table/ps_table_a_" << a << "_n_initial_"
      << m_initial_condition << ".d";

    os_cf << "data/cf_table/cf_table_a_"<< a << "_n_initial_"
      << m_initial_condition << ".d";

    os << "data/powerSpectra_a_" << a << "_n_initial_"<< m_initial_condition
      << ".txt";
  }
  else if (m_initial_condition > 100)
  {
    double m_initial_print = double(m_initial_condition)/100.0;
    os_ps << "data/ps_table/ps_table_a_" << a << "_n_initial_"
      << m_initial_print << ".d";

    os_cf << "data/cf_table/cf_table_a_"<< a << "_n_initial_"
      << m_initial_print << ".d";

    os << "data/powerSpectra_a_" << a << "_n_initial_"<< m_initial_print
      << ".txt";
  }
  std::string ps_table = os_ps.str();
  std::string cf_table = os_cf.str();
  std::string power_file = os.str();

  KFT::iniCorrTable corr_table
    (C->get_power_spectrum (), a_initial, q_min, q_max);

  corr_table.print_tables (ps_table, cf_table);

  KFT::powerSpectra P (C);
  P.initCorrelation (ps_table, cf_table);

  astro::functionWriter write (power_file);
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

void testPowerSpectrum::writeAllSpectrumControl(KFT::kftCosmology * C, double a)
{
  std::ostringstream os, os_ps, os_cf ;
  if (m_initial_condition < 100)
  {
    os_ps << "data/ps_table/ps_table_a_" << a << "_n_initial_"
      << m_initial_condition << ".d";

    os_cf << "data/cf_table/cf_table_a_"<< a << "_n_initial_"
      << m_initial_condition << ".d";

    os << "data/powerSpectra_a_" << a << "_n_initial_"<< m_initial_condition
      << ".txt";
  }
  else if (m_initial_condition > 100)
  {
    double m_initial_print = double(m_initial_condition)/100.0;
    os_ps << "data/ps_table/ps_table_a_" << a << "_n_initial_"
      << m_initial_print << ".d";

    os_cf << "data/cf_table/cf_table_a_"<< a << "_n_initial_"
      << m_initial_print << ".d";

    os << "data/powerSpectra_a_" << a << "_n_initial_"<< m_initial_print
      << ".txt";
  }
  std::string ps_table = os_ps.str();
  std::string cf_table = os_cf.str();
  std::string power_file = os.str();

  // 0 = do correlation table, 1=do not do correlation table and start from scratch
  // 2= start from middle
  int determine = 0;
  if (determine == 0)
  {
    KFT::iniCorrTable corr_table
      (C->get_power_spectrum (), a_initial, q_min, q_max);

    corr_table.print_tables (ps_table, cf_table);
  }
  KFT::powerSpectra P (C);
  P.initCorrelation (ps_table, cf_table);
  std::ofstream outfile;
  if (determine == 1 || determine == 0)
  {
    outfile.open(power_file);
    outfile << "k" << "\t" << "Mean Field" << "\t" << "Linearly evolved" << "\t"
      << "curly P" << "\t" << "Full curly P" << std::endl;
  }
  else if (determine==2)
    outfile.open(power_file, std::ios_base::app);

  for (int i = 0; i < n_bins; i++)
  {
    double k = astro::x_logarithmic(i,n_bins, k_min, 100000.0);
    outfile << k << "\t" << P.meanF(k,a) << "\t" << P.linearP(k,a) << "\t"
      << P.curlyP(k,a) << "\t" << P.BornApproxP(k,a) << std::endl;
    std::cout << "i = " << i << std::endl;
  }
  outfile.close();
}

// For the Gaussian initial conditions
void testPowerSpectrum::writeAllGaussian(KFT::kftCosmology * C, double a)
{
  std::ostringstream os, os_ps, os_cf ;
  os_ps << "data/ps_table/ps_table_a_" << a << "_k0_"<< m_gauss_k0 << "_sigma_"
    << m_sigma_gauss << ".d";

  os_cf << "data/cf_table/cf_table_a_"<< a << "_k0_"<< m_gauss_k0 << "_sigma_"
    << m_sigma_gauss << ".d";

  os << "data/GaussSpectrakmax10000_a_" << a << "_k0_"<< m_gauss_k0 << "_sigma_"
    << m_sigma_gauss << ".txt";
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
  write (k_min, 10000.0, n_bins, astro::LOG_SPACING);
}



