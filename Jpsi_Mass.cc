#include <math.h>

double Jpsi_Mass(double p_T){

  const double mass_p0 = 5.24092392371867076e-07;
  const double mass_p1 = -7.28663700268195047e-06;
  const double mass_p2 = 1.01285390511199336e-04;
  const double mass_p3 = 3.09659981357079550e+00;
  
  double mass;

  mass = ( mass_p0 * p_T * p_T * p_T + mass_p1 * p_T * p_T + 
	   mass_p2 * p_T + mass_p3);

  return mass;
}
