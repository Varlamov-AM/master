#include <math.h>
double Jpsi_Mass(double p_T, int t){

  const double mass_p0_ALICE0 = 1.24448e-07;
  const double mass_p1_ALICE0 = -3.6652e-07;
  const double mass_p2_ALICE0 = 7.89973e-05;
  const double mass_p3_ALICE0 = 3.09661e+00;

  const double mass_p0_ALICE3_1 = 3.09866e+00;
  const double mass_p1_ALICE3_1 = -2.4389e-01;
  const double mass_p2_ALICE3_1 = 1.14869e+02;
  const double mass_p3_ALICE3_1 = -3.5145e+01;

  const double mass_p0_ALICE3_2 = 3.09900e+00;
  const double mass_p1_ALICE3_2 = -2.75641e-01;
  const double mass_p2_ALICE3_2 = -8.30817e-01;
  const double mass_p3_ALICE3_2 = -6.20492e+00;

  const double mass_p0_ALICE3_3 = 3.09810e+00;
  const double mass_p1_ALICE3_3 = -3.3990e-01;
  const double mass_p2_ALICE3_3 = 1.05611e+02;
  const double mass_p3_ALICE3_3 = -4.2814e+01;

  double mass;

  if(t == 0){
    mass = (mass_p0_ALICE0 * p_T * p_T * p_T + mass_p1_ALICE0 * p_T * 
	    p_T +  mass_p2_ALICE0 * p_T + mass_p3_ALICE0);
    return mass;
  }

  if (t == 1){
    mass = (mass_p0_ALICE3_1 * 
	    (1 - exp(mass_p1_ALICE3_1*(p_T - mass_p2_ALICE3_1) + 
		     mass_p3_ALICE3_1))); 
    return mass;
  }

  if (t == 2){
    mass = (mass_p0_ALICE3_2 * 
	    (1 - exp(mass_p1_ALICE3_2*(p_T - mass_p2_ALICE3_2) + 
		     mass_p3_ALICE3_2)));
    return mass;
  }
  
  if (t == 3){
    mass = (mass_p0_ALICE3_3 * 
	    (1 - exp(mass_p1_ALICE3_3*(p_T - mass_p2_ALICE3_3) + 
		     mass_p3_ALICE3_3)));
    return mass;
  }
  
  return 0.;
}
