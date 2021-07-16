#include "TLorentzVector.h"
#include "TRandom.h"
#include <math.h>

Double_t smearP(Double_t Ptrue, Double_t Etrue, int smear_type){
  if (smear_type == 0){
    // Generate smeared track 3-momentum from the true 3-momentum
    const Double_t a=0.008, b=0.002; // Momentum resolution of ALICE treaking system
    // const Double_t a=0.00, b=0.000; // for studies ideal resolution
    Double_t sigmaP = Ptrue * sqrt(a*a + b*Ptrue*b*Ptrue);
    Double_t Psmeared = gRandom->Gaus(Ptrue,sigmaP);
    if (Psmeared<0) Psmeared = 0;
    return Psmeared;
  }
  
  if (smear_type == 1){
    // Generate smeared track 3-momentum from the true 3-momentum
    Double_t sigmaP = Ptrue * 0.015;
    Double_t Psmeared = gRandom->Gaus(Ptrue,sigmaP);
    if (Psmeared<0) Psmeared = 0;
    return Psmeared;
  }

  if (smear_type == 2){
    // Generate smeared track 3-momentum from the true 3-momentum
    const Double_t a = 0.018, b = 0.033, c = 0.011;
    Double_t sigmaP = Ptrue * sqrt(a*a/Etrue/Etrue + b*b/Etrue + c*c);
    Double_t Psmeared = gRandom->Gaus(Ptrue,sigmaP);
    if (Psmeared<0) Psmeared = 0;
    return Psmeared;
  }

  if (smear_type == 3){
    // Generate smeared track 3-momentum from the true 3-momentum
    const Double_t a = 0.001, b = 0.01, c = 0.01;
    Double_t sigmaP = Ptrue * sqrt(a*a/Etrue/Etrue + b*b/Etrue + c*c);
    Double_t Psmeared = gRandom->Gaus(Ptrue,sigmaP);
    if (Psmeared<0) Psmeared = 0;
    return Psmeared;
  }
  return 0;
}
