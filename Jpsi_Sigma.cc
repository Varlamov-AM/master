#include <math.h> 

double Jpsi_Sigma(double p_T){

  const double sigma_p0 = 1.37049718632674920e-02;
  const double sigma_p1 = -1.61076398059093157e-01;
  const double sigma_p2 = 1.42193475032417238e+00;
  const double sigma_p3 = 2.63952203646109461e-03;
  const double sigma_p4 = 2.07789406363807598e-03;


  double sigma;

  sigma = sigma_p0*exp(sigma_p1*(p_T - sigma_p2)) + sigma_p3*p_T + sigma_p4;

  return sigma;
}
