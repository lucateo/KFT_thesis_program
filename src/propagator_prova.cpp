#include <iostream>

#include <astro/cosmology/cosmologyBase.h>
#include <astro/cosmology/cosmicStructures.h>
#include <astro/io/functionWriter.h>
#include <KFT/particleDynamics/newZeldovichParticleDynamics.h>

//#include "../src/particleDynamics/newZeldovichParticleDynamics.h"

int main ()
{
  const double a_min = 0.001;
  astro::cosmologyBase cos_model;
//  cos_model.setDarkUniverse ();
  astro::cosmicStructures cos_struct (&cos_model);
  KFT::newZeldovichParticleDynamics prop (&cos_model, a_min);
  
  astro::functionWriter write ("propagator.txt");
  write.push_back ([&] (double a) { return prop.g_qp (a); });
  write.push_back ([&] (double a) { return prop.g_pp (a); });
  write.push_back ([&] (double a) { return prop.g (a); });
  write.push_back ([&] (double a) { return prop.g_v (a); });
  write.push_back ([&] (double a) { return prop.g_dot (a); });
  write.push_back ([&] (double a) { return prop.tau (a); });
  write.add_header ("# col. 0: scale factor a");
  write.add_header ("# col. 1: g_qp (a)");
  write.add_header ("# col. 2: g_pp (a)");
  write.add_header ("# col. 3: g (a)");
  write.add_header ("# col. 4: g_v (a)");
  write.add_header ("# col. 5: g_dot (a)");
  write.add_header ("# col. 6: tau (a)");
  write (a_min, 1.0, 512, astro::LOG_SPACING);
  
  return 0;
}
