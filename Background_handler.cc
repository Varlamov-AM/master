#include <iostream>
#include <TLorentzVector.h>
#include "Pythia8/Pythia.h"
#include <TH2.h>
#include <string>

using namespace Pythia8;

TLorentzVector resolutionPhoton  (TLorentzVector);
TLorentzVector resolutionElectron(TLorentzVector, int);
bool IsElectronDetectedInCTS(TLorentzVector);
bool IsElectronDetectedInALICE3(TLorentzVector, int);
bool IsPhotonDetectedInPHOS(TLorentzVector);
double Jpsi_Mass(double, int);
double Jpsi_Sigma(double, int);


void Background_handler(Pythia* pythia,
			int* elec, int* posi, int* gamm,
			int elec_num, int posi_num, int gam_num,
			TLorentzVector* elec_data, 
			TLorentzVector* posi_data,
			TLorentzVector* gamma_data,
			TLorentzVector* elec_true_data, 
			TLorentzVector* posi_true_data,
			TLorentzVector* Jpsi_data_ALICE0,
			TLorentzVector* Jpsi_data_ALICE3_1,
			TLorentzVector* Jpsi_data_ALICE3_2,
			TLorentzVector* Jpsi_data_ALICE3_3,
			TH2F* hMassElecPosi_background_ALICE0,
			TH2F* hMassElecPosi_background_ALICE3_1,
			TH2F* hMassElecPosi_background_ALICE3_2,
			TH2F* hMassElecPosi_background_ALICE3_3,
			TH2F* hMassElecPosi_true_background,
			TH2F* hMassGamElecPosi_background_mass_diff_ALICE0, 
			TH2F* hMassGamElecPosi_background_mass_diff_ALICE3_1, 
			TH2F* hMassGamElecPosi_background_mass_diff_ALICE3_2, 
			TH2F* hMassGamElecPosi_background_mass_diff_ALICE3_3,
			double br){

  /* This method is designed to handle background events. As variables, 
     the method receives a pointer to the pythia generator, pointers to 
     arrays of electrons, positrons and photons from the event.

     Then, the following steps are performed:
     
     1. Formation of four-momentum particles.
     2. Smearing four-momentum.
     3. Selection of candidates J/psi.
     4. Forming a combinatorial background, filling background into the 
     corresponding histogram.

  */

  // 4-momentum for particles

  if (elec_num != 0){

    for (int i = 0; i < elec_num; i++){
      TLorentzVector tmp(pythia->event[*(elec + i)].px(), 
			 pythia->event[*(elec + i)].py(), 
			 pythia->event[*(elec + i)].pz(), 
			 pythia->event[*(elec + i)].e());

      *(elec_data + 3 * i + 0) = resolutionElectron(tmp, 0);
      *(elec_data + 3 * i + 1) = resolutionElectron(tmp, 1);
      *(elec_data + 3 * i + 2) = resolutionElectron(tmp, 2);
      *(elec_data + 3 * i + 3) = resolutionElectron(tmp, 2);
      elec_true_data[i] = tmp;
      
    }
  }
  if (posi_num != 0){
    
    for (int i = 0; i < posi_num; i++){
      TLorentzVector tmp(pythia->event[*(posi + i)].px(), 
			 pythia->event[*(posi + i)].py(), 
			 pythia->event[*(posi + i)].pz(), 
			 pythia->event[*(posi + i)].e());
      
      *(posi_data + 3 * i + 0) = resolutionElectron(tmp, 0);
      *(posi_data + 3 * i + 1) = resolutionElectron(tmp, 1);
      *(posi_data + 3 * i + 2) = resolutionElectron(tmp, 2);
      *(posi_data + 3 * i + 3) = resolutionElectron(tmp, 3);
      posi_true_data[i] = tmp;
      
    }
  }
  if (gam_num != 0){
    
    for (int i = 0; i < gam_num; i++){
      TLorentzVector tmp(pythia->event[*(gamm + i)].px(), 
			 pythia->event[*(gamm + i)].py(), 
			 pythia->event[*(gamm + i)].pz(), 
			 pythia->event[*(gamm + i)].e());
      
      gamma_data[i] = resolutionPhoton(tmp);
    }
  }

  if (elec_num != 0 && posi_num != 0){  
    for (int i = 0; i < elec_num; i++){
      for (int j = 0; j < posi_num; j++){
	for (int k = 0; k <= 3; k++){
	  if (k == 0){
	    if (IsElectronDetectedInCTS(elec_true_data[i]) &&
		IsElectronDetectedInCTS(posi_true_data[j])){
	      hMassElecPosi_background_ALICE0->Fill((*(elec_data + 3 * i + k) +
						     *(posi_data + 3 * j + k)).M(),
						    (*(elec_data + 3 * i + k) +
						     *(posi_data + 3 * j + k)).Pt()); 
	    }
	  }
	  if (k == 1){
	    if (IsElectronDetectedInALICE3(elec_true_data[i], 2) &&
		IsElectronDetectedInALICE3(posi_true_data[j], 2)){
	      hMassElecPosi_background_ALICE3_1->Fill((*(elec_data + 3 * i + k) +
						       *(posi_data + 3 * j + k)).M(),
						      (*(elec_data + 3 * i + k) +
						       *(posi_data + 3 * j + k)).Pt()); 
	    }
	  }
	  if (k == 2){
	    if (IsElectronDetectedInALICE3(elec_true_data[i], 2) &&
		IsElectronDetectedInALICE3(posi_true_data[j], 2)){
	      hMassElecPosi_background_ALICE3_2->Fill((*(elec_data + 3 * i + k) +
						       *(posi_data + 3 * j + k)).M(),
						      (*(elec_data + 3 * i + k) +
						       *(posi_data + 3 * j + k)).Pt()); 
	    }
	  }
	  if (k == 3){
	    if (IsElectronDetectedInALICE3(elec_true_data[i], 2) &&
		IsElectronDetectedInALICE3(posi_true_data[j], 2)){
	      hMassElecPosi_background_ALICE3_3->Fill((*(elec_data + 3 * i + k) +
						       *(posi_data + 3 * j + k)).M(),
						      (*(elec_data + 3 * i + k) +
						       *(posi_data + 3 * j + k)).Pt()); 
	    }
	  }
	}
	
	  
	hMassElecPosi_true_background->Fill((elec_true_data[i] + 
					       posi_true_data[j]).M(),
					      (elec_true_data[i] + 
					       posi_true_data[j]).Pt()); 
      }
    }
  }

  
  // find candidates for J/psi

  /* 
     1. Evaluate J/psi mass with defied function Jpsi_Mass(p_T);
     2. Evaluate J/psi sigma with defined function Jpsi_Sigma(p_T);
     3. p_T is taken from (P_{e^+} + P_{e^-}).Pt()
     4. Compare |(P_{e^+} + P_{e^-}).M() - Jpsi_Mass((P_{e^+} + P_{e^-}).Pt())|
     v Jpsi_Sigma((P_{e^+} + P_{e^-}).Pt()). if || < sigma, then forming Jpsi 
     from thats e^{+}e^{-} pair, else continue.
     5. Create arrays of Jpsi TLorentzVectors.

  */

  Int_t Jpsi_num_ALICE0   = 0;
  Int_t Jpsi_num_ALICE3_1 = 0;
  Int_t Jpsi_num_ALICE3_2 = 0;
  Int_t Jpsi_num_ALICE3_3 = 0;
  
  if (elec_num != 0 && posi_num != 0){  
    for (int i = 0; i < elec_num; i++){
      for (int j = 0; j < posi_num; j++){
	double Jpsi_mass[4]; 
	double Jpsi_sigma[4];
	
	for (int l = 0; l < 4; l++){
	  Jpsi_mass[l] = Jpsi_Mass((elec_data[i] + posi_data[j]).Pt(), l);
	  Jpsi_sigma[l] = Jpsi_Sigma((elec_data[i] + posi_data[j]).Pt(), l);
	}
;
	  if (IsElectronDetectedInCTS(elec_true_data[i]) &&
	    IsElectronDetectedInCTS(posi_true_data[j])){
	  if ( fabs((elec_data[i] + posi_data[j]).M() - Jpsi_mass[0]) <= Jpsi_sigma[0]){
	    Jpsi_data_ALICE0[Jpsi_num_ALICE0] = elec_data[i] + posi_data[j];
	    Jpsi_num_ALICE0++;
	  }
	}
	if (IsElectronDetectedInALICE3(elec_true_data[i], 2) &&
	    IsElectronDetectedInALICE3(posi_true_data[j], 2)){
	  if ( fabs((elec_data[i] + posi_data[j]).M() - Jpsi_mass[1]) <= Jpsi_sigma[1]){
	    Jpsi_data_ALICE3_1[Jpsi_num_ALICE3_1] = elec_data[i] + posi_data[j];
	    Jpsi_num_ALICE3_1++;
	  }
	}
	if (IsElectronDetectedInALICE3(elec_true_data[i], 2) &&
	    IsElectronDetectedInALICE3(posi_true_data[j], 2)){
	  if ( fabs((elec_data[i] + posi_data[j]).M() - Jpsi_mass[2]) <= Jpsi_sigma[2]){
	    Jpsi_data_ALICE3_2[Jpsi_num_ALICE3_2] = elec_data[i] + posi_data[j];
	    Jpsi_num_ALICE3_2++;
	  }
	}
	if (IsElectronDetectedInALICE3(elec_true_data[i], 2) &&
	    IsElectronDetectedInALICE3(posi_true_data[j], 2)){
	  if ( fabs((elec_data[i] + posi_data[j]).M() - Jpsi_mass[3]) <= Jpsi_sigma[3]){
	    Jpsi_data_ALICE3_3[Jpsi_num_ALICE3_3] = elec_data[i] + posi_data[j];
	    Jpsi_num_ALICE3_3++;
	  }
	}
      }
    }
  }
  

  if (Jpsi_num_ALICE0 != 0 && gam_num != 0){
    for (int i = 0; i < Jpsi_num_ALICE0; i++){ 
      for (int j = 0; j < gam_num; j++){
	if(IsPhotonDetectedInPHOS(gamma_data[j])){
	  hMassGamElecPosi_background_mass_diff_ALICE0->Fill((Jpsi_data_ALICE0[i] + gamma_data[j]).M() - Jpsi_data_ALICE0[i].M(), (Jpsi_data_ALICE0[i] + gamma_data[j]).Pt(), br);
	}
      }
    }   
  } 

  if (Jpsi_num_ALICE3_1 != 0 && gam_num != 0){
    for (int i = 0; i < Jpsi_num_ALICE3_1; i++){ 
      for (int j = 0; j < gam_num; j++){
	if(IsElectronDetectedInALICE3(gamma_data[j], 1)){
	  hMassGamElecPosi_background_mass_diff_ALICE3_1->Fill((Jpsi_data_ALICE3_1[i] + gamma_data[j]).M() - Jpsi_data_ALICE3_1[i].M(), (Jpsi_data_ALICE3_1[i] + gamma_data[j]).Pt(), br);
	}
      }
    }   
  } 

  if (Jpsi_num_ALICE3_2 != 0 && gam_num != 0){
    for (int i = 0; i < Jpsi_num_ALICE3_2; i++){ 
      for (int j = 0; j < gam_num; j++){
	if(IsElectronDetectedInALICE3(gamma_data[j], 1)){
	  hMassGamElecPosi_background_mass_diff_ALICE3_2->Fill((Jpsi_data_ALICE3_2[i] + gamma_data[j]).M() - Jpsi_data_ALICE3_2[i].M(), (Jpsi_data_ALICE3_2[i] + gamma_data[j]).Pt(), br);
	}
      }
    }   
  }  

  if (Jpsi_num_ALICE3_3 != 0 && gam_num != 0){
    for (int i = 0; i < Jpsi_num_ALICE3_3; i++){ 
      for (int j = 0; j < gam_num; j++){
	if(IsElectronDetectedInALICE3(gamma_data[j], 1)){
	  hMassGamElecPosi_background_mass_diff_ALICE3_3->Fill((Jpsi_data_ALICE3_3[i] + gamma_data[j]).M() - Jpsi_data_ALICE3_3[i].M(), (Jpsi_data_ALICE3_3[i] + gamma_data[j]).Pt(), br);
	}
      }
    }   
  }
}
