#include <math.h> 

double Jpsi_Sigma(double p_T, int t){

  const double sigma_p0_ALICE0 =  3.48896e-02;
  const double sigma_p1_ALICE0 = -1.89391e-01;
  const double sigma_p2_ALICE0 = -4.74279e+00;
  const double sigma_p3_ALICE0 =  2.45733e-03;
  const double sigma_p4_ALICE0 =  5.17263e-03;
  
  const double sigma_p0_ALICE3_1 = 3.29870e-02;
  
  const double sigma_p0_ALICE3_2 =  4.65884e-02;
  const double sigma_p1_ALICE3_2 =  1.69952e+02;
  const double sigma_p2_ALICE3_2 = -1.25253e+01;

  const double sigma_p0_ALICE3_3 =  3.06009e-02; 
  const double sigma_p1_ALICE3_3 =  2.87757e+02; 
  const double sigma_p2_ALICE3_3 = -1.68239e+01;

  double sigma;

  if (t == 0){
    
    sigma = (sigma_p0_ALICE0 * 
	     exp(sigma_p1_ALICE0*(p_T - sigma_p2_ALICE0)) + 
	     sigma_p3_ALICE0*p_T + sigma_p4_ALICE0);
    
    return sigma;
  }

  if (t == 1){
    
    sigma =  sigma_p0_ALICE3_1;
    
    return sigma;
  }

  if (t == 2){
    
    sigma = sigma_p0_ALICE3_2/(sqrt(sqrt(p_T + sigma_p1_ALICE3_2)
				    + sigma_p2_ALICE3_2));
    
    return sigma;
  }

  if (t == 3){
    
    sigma = sigma_p0_ALICE3_3/(sqrt(sqrt(sqrt(p_T + sigma_p1_ALICE3_3)
					 + sigma_p2_ALICE3_3)));
    
    return sigma;
  }

  return 0.;
}
