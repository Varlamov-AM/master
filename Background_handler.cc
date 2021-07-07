#include <iostream>
#include <TLorentzVector.h>
#include "Pythia8/Pythia.h"
#include <TH2.h>

using namespace Pythia8;

TLorentzVector resolutionPhoton  (TLorentzVector);
TLorentzVector resolutionElectron(TLorentzVector);
bool IsElectronDetectedInCTS(TLorentzVector);

void Background_handler(Pythia* pythia,
			int* elec, int* posi, int* gamm,
			int elec_num, int posi_num, int gam_num,
			TLorentzVector* elec_data, TLorentzVector* posi_data,
			TLorentzVector* gamma_data,
			TLorentzVector* elec_true_data, TLorentzVector* posi_true_data,
			TH2F* hMassElecPosi_background,
			TH2F* hMassElecPosi_true_background){

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

      elec_data[i] = resolutionElectron(tmp);
      elec_true_data[i] = tmp;
      
    }
  }
  if (posi_num != 0){
    
    for (int i = 0; i < posi_num; i++){
      TLorentzVector tmp(pythia->event[*(posi + i)].px(), 
			 pythia->event[*(posi + i)].py(), 
			 pythia->event[*(posi + i)].pz(), 
			 pythia->event[*(posi + i)].e());
      
      posi_data[i] = resolutionElectron(tmp);
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
	if (IsElectronDetectedInCTS(elec_true_data[i]) &&
	    IsElectronDetectedInCTS(posi_true_data[j])){
	  hMassElecPosi_background->Fill((elec_data[i] + posi_data[j]).M(),
					 (elec_data[i] + posi_data[j]).Pt()); 
	  hMassElecPosi_true_background->Fill((elec_true_data[i] + posi_true_data[j]).M(),
					      (elec_true_data[i] + posi_true_data[j]).Pt()); 
	}
      }
    }
  }

  // find candidates for J/psi


}
