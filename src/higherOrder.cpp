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
  return P_31_1*P_32_1 + P_21_2*P_32_2 + P_31_3*P_21_3;
}


void testPowerSpectrum::writeAllHigherOrder(KFT::kftCosmology * C, double a,
    double k2, double mu)
{
  std::ostringstream os, os_ps, os_cf ;
  os_ps << "data/ps_table/ps_table_a_" << a << "_n_initial_"
    << m_initial_condition << ".d";

  os_cf << "data/cf_table/cf_table_a_"<< a << "_n_initial_"
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

  double Fixed_mean = P.meanF(k2, a);

  astro::functionWriter write (power_file);
  write.push_back ([&] (double k) { return P.meanF (k, a); });
  write.push_back ([&] (double k) { return P.curlyP_ij (a, -mu,
        mu ,-1,-k,-k-k2, k); });
  write.push_back ([&] (double k) { return Fixed_mean; });
  write.add_header ("# Different types of cosmic density-fluctuation");
  write.add_header ("# power spectra as functions of wave number");
  write.add_header ("# column 1: wave number k in h/Mpc");
  write.add_header ("# column 2: mean-field interaction term");
  write.add_header ("# column 3: curly P power spectrum");
  write.add_header ("# column 4: fixed mean term");
  write (k_min, k_max, n_bins, astro::LOG_SPACING);
}

void testPowerSpectrum::writeBiSpectrum(KFT::kftCosmology * C, double a,
    double k_prime, double mu)
{
  std::ostringstream os, os_ps, os_cf ;
  os_ps << "data/ps_table/ps_table_a_" << a << "_n_initial_"
    << m_initial_condition << ".d";

  os_cf << "data/cf_table/cf_table_a_"<< a << "_n_initial_"
    << m_initial_condition << ".d";

  os << "data/powerSpectraHigherOrder_a_" << a << "_n_initial_"<< m_initial_condition
    << "_kprime_" << k_prime << "_mu_" << mu << ".txt";

  std::string ps_table = os_ps.str();
  std::string cf_table = os_cf.str();
  std::string power_file = os.str();

  KFT::iniCorrTable corr_table
    (C->get_power_spectrum (), a_initial, q_min, q_max);

  corr_table.print_tables (ps_table, cf_table);

  powerSpectraModified P (C);
  P.initCorrelation (ps_table, cf_table);

  double Fixed_mean = P.meanF(k_prime, a);

  astro::functionWriter write (power_file);
  write.push_back ([&] (double k) { return P.meanF (k, a); });
  write.push_back ([&] (double k) { return P.bispectrumFree (k, a, k_prime,mu); });
  write.push_back ([&] (double k) { return Fixed_mean; });
  write.add_header ("# Different types of cosmic density-fluctuation");
  write.add_header ("# power spectra as functions of wave number");
  write.add_header ("# column 1: wave number k in h/Mpc");
  write.add_header ("# column 2: mean-field interaction term");
  write.add_header ("# column 3: curly P power spectrum");
  write.add_header ("# column 4: fixed mean term");
  write (k_min, k_max, n_bins, astro::LOG_SPACING);
};
