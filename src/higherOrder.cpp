#include "../include/KFT2.h"

powerSpectraModified::powerSpectraModified(KFT::kftCosmology * cosmology_in):
  powerSpectra(cosmology_in) {}

double powerSpectraModified::curlyP_ij (double a, double mu_ij, double mu_i,
    double mu_j, double K_ij, double L_i, double L_j)
 {
   KFT::CorrLine::funcArgs arguments;
   arguments.mu_ij = mu_ij;
   arguments.mu_i  =  mu_i;
   arguments.mu_j  = mu_j;
   KFT::MomentumCorrLine corr_line (*corr_interp, arguments, n_bins, q_min, q_lim);
   arguments.K_ij = K_ij;
   double g_qp = cosmology->g_qp (a);
   arguments.L_i  = g_qp*L_i;
   arguments.L_j  = g_qp*L_j;
   return
     gsl_pow_2 (cosmology->growth_factor (a)/cosmology->g_qp (a))*
     corr_line ();
}

double powerSpectraModified::bispectrumFree (double k1, double a,
    double k2, double mu)
{
  double k3_module = sqrt(k1*k1 + k2*k2);
  double factor1 = k1*k1 + k1*k2*mu;
  double factor2 = k2*k2 + k2*k1*mu;

  double P_31_0 = curlyP_ij(a, -1, 1, -1, k1, k1,k1);
  double P_21_0 = curlyP_ij(a, mu, 1, mu, k2, k2,k1);
  double P_32_0 = curlyP_ij(a, -1, 1, -1, k2, k2,k2);

  double P_31_1= curlyP_ij(a, -factor1/(k1*k3_module) , factor1/(k1*k3_module) ,
      -1,k1,k1+k2, k1);
  double P_32_1=curlyP_ij(a,-factor2/(k2*k3_module), -factor2/(k2*k3_module),
      1,k2,k1+k2,k2);

  double P_21_2= curlyP_ij(a, mu, -mu, -1, k1, k2,k1);
  double P_32_2= curlyP_ij(a, -factor2/(k2*k3_module), -1,factor2/(k2*k3_module),
      k1+k2, +k2+k1,k2);
  
  double P_31_3= curlyP_ij(a,-factor1/(k1*k3_module),1,-factor1/(k1*k3_module),
      +k2+k1,+k2+k1,k1);
  double P_21_3= curlyP_ij(a, mu, 1,mu,k2,k2,k1);
  
  // To reproduce Bernardeau Q_N
  double P1 = curlyP(k1,a);
  double P2 = curlyP(k2,a);
  double P3 = curlyP(k3_module,a);
  double Q_N_denominator = P1*P2 + P2*P3 + P1*P3;
  
  double bispectrum = 2*M_PI*(P_31_0 + P_21_0 + P_32_0) +
    P_31_1*P_32_1 + P_21_2*P_32_2 + P_31_3*P_21_3;
  return bispectrum/Q_N_denominator;
}

// I should integrate also over angles, crazy!
double powerSpectraModified::PPP_term (double k1, double a,
    double k2, double mu)
{
  std::function<double (double)> PPP = [this,k1,a,k2,mu](double k_32)
  {
    double k3_module = sqrt(k1*k1 + k2*k2);
    double factor1 = k1*k1 + k1*k2*mu;
    double factor2 = k2*k2 + k2*k1*mu;
    double P_32 = curlyP_ij(a, -factor1/(k2*k3_module), 1,mu,k2,k2,k1);
    double P_21 =1;
    double P_31  = 1;
    return P_32*P_21*P_31/(8*M_PI*M_PI*M_PI); 
  };

  astro::integrator kernel (PPP);
  return kernel(q_min,q_max);
}


void testPowerSpectrum::writeAllHigherOrder(KFT::kftCosmology * C, double a,
    double k2, double mu)
{
  std::ostringstream os, os_ps, os_cf ;
  os_ps << "data/ps_table/ps_tableHigherOrder_a_" << a << "_n_initial_"
    << m_initial_condition << ".d";

  os_cf << "data/cf_table/cf_tableHigherOrder_a_"<< a << "_n_initial_"
    << m_initial_condition << ".d";

  os << "data/powerSpectraHigherOrder_a_" << a << "_n_initial_"<< m_initial_condition
    << ".txt";

  std::string ps_table = os_ps.str();
  std::string cf_table = os_cf.str();
  std::string power_file = os.str();

  KFT::iniCorrTable corr_table
    (C->get_power_spectrum (), a_initial, q_min, q_max);

  corr_table.print_tables (ps_table, cf_table);

  powerSpectraModified P (C);
  P.initCorrelation (ps_table, cf_table);

  // double Fixed_mean = P.meanF(k2, a);

  astro::functionWriter write (power_file);
  // write.push_back ([&] (double k) { return P.meanF (k, a); });
  write.push_back ([&] (double k) { return P.curlyP_ij (a, -mu,
        mu ,-1,-k,-k-k2, k); });
  // write.push_back ([&] (double k) { return Fixed_mean; });
  write.add_header ("# Different types of cosmic density-fluctuation");
  write.add_header ("# power spectra as functions of wave number");
  write.add_header ("# column 1: wave number k in h/Mpc");
  // write.add_header ("# column 2: mean-field interaction term");
  write.add_header ("# column 2: curly P power spectrum");
  // write.add_header ("# column 4: fixed mean term");
  write (k_min, k_max, n_bins, astro::LOG_SPACING);
}

void testPowerSpectrum::writeBiSpectrum(KFT::kftCosmology * C, double a,
    double ratio)
{
  double k = 2.5; // Taken from Bernardeau, ratio = k1/k2
  std::ostringstream os, os_ps, os_cf ;
  if (m_initial_condition ==0)
  {
    os_ps << "data/ps_table/ps_tableHigherOrder_a_" << a << "_sigma_"
      << m_sigma_gauss << "_k0_" << m_gauss_k0 << ".d";

    os_cf << "data/cf_table/cf_tableHigherOrder_a_"<< a << "_sigma_"
      << m_sigma_gauss << "_k0_" << m_gauss_k0 << ".d";

    os << "data/powerSpectraHigherOrder_a_" << a << "_sigma_" << m_sigma_gauss 
      << "_k0_" << m_gauss_k0 <<"_ratio_" << ratio << "_k1_" << k << ".txt";
  }
  else
  {
    os_ps << "data/ps_table/ps_tableHigherOrder_a_" << a << "_n_initial_"
      << m_initial_condition << ".d";

    os_cf << "data/cf_table/cf_tableHigherOrder_a_"<< a << "_n_initial_"
      << m_initial_condition << ".d";

    os << "data/powerSpectraHigherOrder_a_" << a << "_n_initial_"<< m_initial_condition
      << "_ratio_" << ratio << "_k1_" << k << ".txt";
  }
  std::string ps_table = os_ps.str();
  std::string cf_table = os_cf.str();
  std::string power_file = os.str();

  KFT::iniCorrTable corr_table
    (C->get_power_spectrum (), a_initial, q_min, q_max);

  corr_table.print_tables (ps_table, cf_table);

  powerSpectraModified P (C);
  P.initCorrelation (ps_table, cf_table);

  // double Fixed_mean = P.meanF(k_prime, a);

  astro::functionWriter write (power_file);
  // write.push_back ([&] (double k) { return P.meanF (k, a); });
  write.push_back ([&] (double mu) { return P.bispectrumFree (k, a, k/ratio,mu); });
  // write.push_back ([&] (double k) { return Fixed_mean; });
  write.add_header ("# Different types of cosmic density-fluctuation");
  write.add_header ("# power spectra as functions of wave number");
  write.add_header ("# column 1: wave number k in h/Mpc");
  // write.add_header ("# column 2: mean-field interaction term");
  write.add_header ("# column 2: curly P power spectrum");
  // write.add_header ("# column 4: fixed mean term");
  write (-1, 1, n_bins, astro::LINEAR_SPACING);
}

void testPowerSpectrum::writeBiSpectrumFull(KFT::kftCosmology * C, double a,
    double ratio)
{
  double k1 = 2.5; // Taken from Bernardeau, ratio = k1/k2
  std::ostringstream os, os_ps, os_cf ;
  if (m_initial_condition ==0)
  {
    os_ps << "data/ps_table/ps_tableHigherOrder_a_" << a << "_sigma_"
      << m_sigma_gauss << "_k0_" << m_gauss_k0 << ".d";

    os_cf << "data/cf_table/cf_tableHigherOrder_a_"<< a << "_sigma_"
      << m_sigma_gauss << "_k0_" << m_gauss_k0 << ".d";

    os << "data/powerSpectraHigherOrder_a_" << a << "_sigma_" << m_sigma_gauss 
      << "_k0_" << m_gauss_k0 <<"_ratio_" << ratio << "_k1_" << k1 << ".txt";
  }
  else
  {
    os_ps << "data/ps_table/ps_tableHigherOrder_a_" << a << "_n_initial_"
      << m_initial_condition << ".d";

    os_cf << "data/cf_table/cf_tableHigherOrder_a_"<< a << "_n_initial_"
      << m_initial_condition << ".d";

    os << "data/powerSpectraHigherOrder_a_" << a << "_n_initial_"<< m_initial_condition
      << "_ratio_" << ratio << "_k1_" << k1 << ".txt";
  }
  std::string ps_table = os_ps.str();
  std::string cf_table = os_cf.str();
  std::string power_file = os.str();
  // 0 = do correlation table, 1=do not do correlation table and start from scratch
  // 2= start from middle
  int determine = 2;
  if (determine == 0)
  {
    KFT::iniCorrTable corr_table
      (C->get_power_spectrum (), a_initial, q_min, q_max);

    corr_table.print_tables (ps_table, cf_table);
  }
  powerSpectraModified P (C);
  P.initCorrelation (ps_table, cf_table);
  double k2 = double(double(k1)/double(ratio));
  std::ofstream outfile;
  // outfile.open(power_file, std::ios_base::app);
  if (determine == 1 || determine == 0)
  {
    outfile.open(power_file);
    outfile << "mu" << "\t" << "First Factor" << "\t" << "Second Factor" << "\t"
      << "QN denominator" << "\t" << "Bispectrum" << std::endl;
  }
  else if (determine==2)
    outfile.open(power_file, std::ios_base::app);

  for (int i = 79; i < n_bins; i++)
  {
    double mu = astro::x_linear(i,n_bins,-1.0,1.0);
    double k3_module = sqrt(k1*k1 + k2*k2+ 2*mu*k1*k2);
    double factor1 = k1*k1 + k1*k2*mu;
    double factor2 = k2*k2 + k2*k1*mu;

    double P_31_0 = P.curlyP_ij(a, -1, 1, -1, k1, k1,k1);
    double P_21_0 = P.curlyP_ij(a, mu, 1, mu, k2, k2,k1);
    double P_32_0 = P.curlyP_ij(a, -1, 1, -1, k2, k2,k2);

    double P_31_1= P.curlyP_ij(a, -factor1/(k1*k3_module) , factor1/(k1*k3_module) ,
        -1,k1,k1+k2, k1);
    std::cout << "number = " << i << " mu = " << mu << " factor 2 = " << factor2 << " k3 = " << k3_module << std::endl;
    double P_32_1=P.curlyP_ij(a,-factor2/(k2*k3_module), -factor2/(k2*k3_module),
        1,k2,k1+k2,k2);

    double P_21_2= P.curlyP_ij(a, mu, -mu, -1, k1, k2,k1);
    double P_32_2= P.curlyP_ij(a, -factor2/(k2*k3_module), -1,factor2/(k2*k3_module),
        k1+k2, +k2+k1,k2);
    
    double P_31_3= P.curlyP_ij(a,-factor1/(k1*k3_module),1,-factor1/(k1*k3_module),
        +k2+k1,+k2+k1,k1);
    double P_21_3= P.curlyP_ij(a, mu, 1,mu,k2,k2,k1);
    // To reproduce Bernardeau Q_N
    double P1 = P.curlyP(k1,a);
    double P2 = P.curlyP(k2,a);
    double P3 = P.curlyP(k3_module,a);
    double Q_N_denominator = P1*P2 + P2*P3 + P1*P3;
    double first_factors = 2*M_PI*(P_31_0 + P_21_0 + P_32_0);
    double second_factors = P_31_1*P_32_1 + P_21_2*P_32_2 + P_31_3*P_21_3;
    double bispectrum = first_factors + second_factors;

    outfile << mu << "\t" << first_factors << "\t" << second_factors << "\t"
      << Q_N_denominator << "\t" << bispectrum << std::endl;
  }
  outfile.close();
};
