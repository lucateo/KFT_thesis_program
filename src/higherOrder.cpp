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

double powerSpectraModified::curlyP_ijNoDamping (double a, double mu_ij, double mu_i,
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
   double damping = corr_interp->get_sigma_p2()*mu_ij*g_qp*g_qp*L_i*L_j;
   return
     gsl_pow_2 (cosmology->growth_factor (a)/cosmology->g_qp (a))*
     exp(-damping)*corr_line ();
}
// In other function, to remove damping change sign to k_factor
double powerSpectraModified::DampingBispectrum(double a, double k_factor)
{
  double g_qp = cosmology->g_qp(a);
  return exp(-corr_interp->get_sigma_p2()*g_qp*g_qp*k_factor);
}

void powerSpectraModified::Trial(double a, double k)
{
  double g_qp = cosmology->g_qp(a);
  double sigma_1 = corr_interp->get_sigma_p2();
  double other = g_qp*g_qp*k*k;
  double all = exp(-sigma_1*other);
  std::cout << "g_qp = " << g_qp << " sigma1 = " << sigma_1 << " other = " <<
    other << " all = " << all << std::endl;
}

// I should integrate also over angles, crazy! It is incomplete
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
    double k2, double mu, int determine, int i_initial)
{
  std::ostringstream os, os_ps, os_cf ;
  os_ps << "data/ps_table/ps_tableHigherOrder_n_initial_"
    << m_initial_condition << ".d";

  os_cf << "data/cf_table/cf_tableHigherOrder_n_initial_"
    << m_initial_condition << ".d";

  os << "data/powerSpectraCurlyPij_a_" << a << "_n_initial_"<< m_initial_condition
    << "_k2_" << k2 << "_mu_"<< mu  << ".txt";

  std::string ps_table = os_ps.str();
  std::string cf_table = os_cf.str();
  std::string power_file = os.str();

  if (determine == 0)
  {
    KFT::iniCorrTable corr_table
      (C->get_power_spectrum (), a_initial, q_min, q_max);

    corr_table.print_tables (ps_table, cf_table);
  }
  powerSpectraModified P (C);
  P.initCorrelation (ps_table, cf_table);
  std::ofstream outfile;
  if (determine == 1 || determine == 0)
  {
    outfile.open(power_file);
    outfile << "k" << "\t" << "P1" << "\t" << "P2" << "\t"
      << "P3" << "\t" << "P4" << std::endl;
  }
  else if (determine==2)
    outfile.open(power_file, std::ios_base::app);

  for (int i = i_initial; i < n_bins; i++)
  {
    double k1 = astro::x_logarithmic(i,n_bins,k_min,k_max);
    double k3_module = sqrt(k1*k1 + k2*k2+ 2*mu*k1*k2);
    double factor1 = k1*k1 + k1*k2*mu;
    double factor2 = k2*k2 + k2*k1*mu;

    double P_31_1= P.curlyP_ij(a, -factor1/(k1*k3_module) , factor1/(k1*k3_module) ,
        -1,k1,k3_module, k1);
    std::cout << "number = " << i << std::endl;
    double P_32_1=P.curlyP_ij(a,-factor2/(k2*k3_module), factor2/(k2*k3_module),
        -1,k2,k3_module,k2);

    double P_21_2= P.curlyP_ij(a, mu, -mu, -1, k1, k2,k1);
    double P_32_2= P.curlyP_ij(a, -factor2/(k2*k3_module), 1,-factor2/(k2*k3_module),
        k3_module, k3_module,k2);

    outfile << k1 << "\t" << P_31_1 << "\t" << P_32_1 << "\t"
      << P_21_2 << "\t" << P_32_2 << std::endl;
  }
  outfile.close();
}

double F2_kernel(double k1,double k2, double mu)
{
  return 0.5 +0.5*mu*(k1/k2 + k2/k1)+ 0.5*mu*mu;
}

void testPowerSpectrum::writeBiSpectrumFull(KFT::kftCosmology * C, double a,
    double ratio, int determine, int i_initial)
{
  double k1 = 2.5; // Taken from Bernardeau, ratio = k1/k2
  std::ostringstream os, os_ps, os_cf ;
  if (m_initial_condition ==0)
  {
    os_ps << "data/ps_table/ps_tableHigherOrder_sigma_"
      << m_sigma_gauss << "_k0_" << m_gauss_k0 << ".d";

    os_cf << "data/cf_table/cf_tableHigherOrder_sigma_"
      << m_sigma_gauss << "_k0_" << m_gauss_k0 << ".d";

    os << "data/powerSpectraHigherOrder_a_"<< a << "_sigma_" << m_sigma_gauss 
      << "_k0_" << m_gauss_k0 <<"_ratio_" << ratio << "_k1_" << k1 << ".txt";
  }
  else
  {
    os_ps << "data/ps_table/ps_tableHigherOrder_n_initial_"
      << m_initial_condition << ".d";

    os_cf << "data/cf_table/cf_tableHigherOrder_n_initial_"
      << m_initial_condition << ".d";

    os << "data/powerSpectraHigherOrder_a_" << a <<"_n_initial_"<< m_initial_condition
      << "_ratio_" << ratio << "_k1_" << k1 << ".txt";
  }
  std::string ps_table = os_ps.str();
  std::string cf_table = os_cf.str();
  std::string power_file = os.str();
  // 0 = do correlation table, 1=do not do correlation table and start from scratch
  // 2= start from middle

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

  for (int i = i_initial; i < n_bins; i++)
  {
    double mu = astro::x_linear(i,n_bins,-1.0,1.0);
    double k3_module = sqrt(k1*k1 + k2*k2+ 2*mu*k1*k2);
    double factor1 = k1*k1 + k1*k2*mu;
    double factor2 = k2*k2 + k2*k1*mu;

    double P_31_1= P.curlyP_ij(a, -factor1/(k1*k3_module) , factor1/(k1*k3_module) ,
        -1,k1,k3_module, k1);
    std::cout << "number = " << i << " mu = " << mu << " factor 2 = " << factor2 << " k3 = " << k3_module << std::endl;
    double P_32_1=P.curlyP_ij(a,-factor2/(k2*k3_module), factor2/(k2*k3_module),
        -1,k2,k3_module,k2);

    double P_21_2= P.curlyP_ij(a, mu, -mu, -1, k1, k2,k1);
    double P_32_2= P.curlyP_ij(a, -factor2/(k2*k3_module), 1,-factor2/(k2*k3_module),
        k3_module, k3_module,k2);
    
    double P_31_3= P.curlyP_ij(a,-factor1/(k1*k3_module),1,-factor1/(k1*k3_module),
        k3_module,k3_module,k1);
    double P_21_3= P.curlyP_ij(a, mu, 1,mu,k2,k2,k1);
    
    double P1 = P.curlyP(k1,a);
    double P2 = P.curlyP(k2,a);
    double P3 = P.curlyP(k3_module,a);
    double Q_N_denominator = P1*P2 + P2*P3 + P1*P3;
    double bispectrum = P_31_1*P_32_1 + P_21_2*P_32_2 + P_31_3*P_21_3;
    
    // For linear Eulerian
    double mu_31 = -(mu*k2*k1 + k1*k1)/(k1*k3_module);
    double mu_32 = -(mu*k2*k1 + k2*k2)/(k2*k3_module);
    double linear_Eulerian = operator()(k1,a)*operator()(k2,a)*F2_kernel(k1,k2,mu)
      +operator()(k2,a)*operator()(k3_module,a)*F2_kernel(k2,k3_module,mu_32) 
      + operator()(k1,a)*operator()(k3_module,a)*F2_kernel(k1,k3_module,mu_31);
    linear_Eulerian = linear_Eulerian/(operator()(k1,a)*operator()(k2,a) +
        operator()(k2,a)*operator()(k3_module,a) + operator()(k1,a)*operator()(k3_module,a));

        outfile << mu << "\t" << linear_Eulerian << "\t" << bispectrum << "\t"
      << Q_N_denominator << "\t" << bispectrum << std::endl;
  }
  outfile.close();
}

void testPowerSpectrum::writeBiSpectrumFixExpDamping(KFT::kftCosmology * C, double a,
    double ratio, int determine, int i_initial)
{
  double k1 = 2.5; // Taken from Bernardeau, ratio = k1/k2
  std::ostringstream os, os_ps, os_cf ;
  if (m_initial_condition ==0)
  {
    os_ps << "data/ps_table/ps_tableHigherOrder_sigma_"
      << m_sigma_gauss << "_k0_" << m_gauss_k0 << ".d";

    os_cf << "data/cf_table/cf_tableHigherOrder_sigma_"
      << m_sigma_gauss << "_k0_" << m_gauss_k0 << ".d";

    os << "data/powerSpectraHigherOrderFixDamping_a_" << a << "_sigma_" << m_sigma_gauss 
      << "_k0_" << m_gauss_k0 <<"_ratio_" << ratio << "_k1_" << k1 << ".txt";
  }
  else
  {
    os_ps << "data/ps_table/ps_tableHigherOrder_n_initial_"
      << m_initial_condition << ".d";

    os_cf << "data/cf_table/cf_tableHigherOrder_n_initial_"
      << m_initial_condition << ".d";

    os << "data/powerSpectraHigherOrderFixDamping_a_" << a << "_n_initial_"<< m_initial_condition
      << "_ratio_" << ratio << "_k1_" << k1 << ".txt";
  }
  std::string ps_table = os_ps.str();
  std::string cf_table = os_cf.str();
  std::string power_file = os.str();
  // 0 = do correlation table, 1=do not do correlation table and start from scratch
  // 2= start from middle
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
  if (determine == 1 || determine == 0)
  {
    outfile.open(power_file);
    outfile << "mu" << "\t" << "Damping corrected" << "\t" << "Example curly P_ij" << "\t"
      << "QN denominator" << "\t" << "Bispectrum" << std::endl;
  }
  else if (determine==2)
    outfile.open(power_file, std::ios_base::app);

  for (int i = i_initial; i < n_bins; i++)
  {
    double mu = astro::x_linear(i,n_bins,-1.0,1.0);
    double k3_module = sqrt(k1*k1 + k2*k2+ 2*mu*k1*k2);
    double factor1 = k1*k1 + k1*k2*mu;
    double factor2 = k2*k2 + k2*k1*mu;

    double P_31_1= P.curlyP_ij(a, -factor1/(k1*k3_module) , factor1/(k1*k3_module) ,
        -1,k1,k3_module, k1);
    std::cout << "number = " << i << " mu = " << mu << " factor 2 = " << factor2 << " k3 = " << k3_module << std::endl;
    double P_32_1=P.curlyP_ij(a,-factor2/(k2*k3_module), factor2/(k2*k3_module),
        -1,k2,k3_module,k2);
    double damping1 = P.DampingBispectrum(a,0);

    double P_21_2= P.curlyP_ij(a, mu, -mu, -1, k1, k2,k1);
    double P_32_2= P.curlyP_ij(a, -factor2/(k2*k3_module), 1,-factor2/(k2*k3_module),
        k3_module, k3_module,k2);
    double damping2 = P.DampingBispectrum(a,k1*k2*mu + k1*k1);
    
    double P_31_3= P.curlyP_ij(a,-factor1/(k1*k3_module),1,-factor1/(k1*k3_module),
        k3_module,k3_module,k1);
    double P_21_3= P.curlyP_ij(a, mu, 1,mu,k2,k2,k1);
    double damping3 = P.DampingBispectrum(a,k1*k2*mu + k2*k2);
    
    // To reproduce Bernardeau Q_N
    double P1 = P.curlyP(k1,a);
    double P2 = P.curlyP(k2,a);
    double P3 = P.curlyP(k3_module,a);
    double Q_N_denominator = P1*P2 + P2*P3 + P1*P3;
    double bispectrum = P_31_1*P_32_1 + P_21_2*P_32_2 + P_31_3*P_21_3;
    double damping_corrected = damping1*P_31_1*P_32_1 + damping2*P_21_2*P_32_2 
      + damping3*P_31_3*P_21_3;

    outfile << mu << "\t" << damping_corrected << "\t" << P_31_3 << "\t"
      << Q_N_denominator << "\t" << bispectrum << std::endl;
  }
  outfile.close();
}
