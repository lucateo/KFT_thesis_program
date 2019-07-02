#include "../include/KFT2.h"

powerSpectraModified::powerSpectraModified(KFT::kftCosmology * cosmology_in):
  powerSpectra(cosmology_in) {}

double powerSpectraModified::curlyP_ij (double k, double a, double l_parallel)
 {
   KFT::CorrLine::funcArgs arguments;
   arguments.mu_ij = l_parallel;
   arguments.mu_i  =  1.0;
   arguments.mu_j  = -1.0;
   KFT::MomentumCorrLine corr_line (*corr_interp, arguments, n_bins, q_min, q_lim);
   arguments.K_ij = k;
   double g_qp = cosmology->g_qp (a);
   arguments.L_i  = g_qp*k;
   arguments.L_j  = arguments.L_i;
   return
     gsl_pow_2 (cosmology->growth_factor (a)/cosmology->g_qp (a))*
     corr_line ();
}

double powerSpectraModified::bispectrumFree (double k, double a,
    double k_prime, double mu)
 {
   double factor1 = -(k_prime * k * mu + k_prime * k_prime )/(k*k + k_prime*k_prime + 2*k*k_prime*mu);
   double factor2 = -(k_prime * k * mu + k* k)/(k*k + k_prime*k_prime + 2*k*k_prime*mu);
   return curlyP_ij(k,a, -k_prime * k * mu / (k*k) )* curlyP_ij(k,a, -k_prime * k * mu / (k_prime*k_prime) )
     + curlyP_ij(k,a, k_prime * k * mu / (k*k) ) * curlyP_ij(k+ k_prime,a, factor1) + 
     curlyP_ij(-k_prime,a, k_prime * k * mu / (k*k) )*curlyP_ij(k+k_prime,a, factor2 );
 }


void testPowerSpectrum::writeAllHigherOrder(KFT::kftCosmology * C, double a,
    double k_prime, double l_parallel)
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

  double Fixed_mean = P.meanF(k_prime, a);

  astro::functionWriter write (power_file);
  write.push_back ([&] (double k) { return P.meanF (k, a); });
  write.push_back ([&] (double k) { return P.curlyP_ij (k, a, l_parallel); });
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
    << ".txt";

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
