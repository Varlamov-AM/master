#include <math.h> 
#include <iostream>


void Slice_and_Fit_Jpsi_inv_mass(){

	TFile *file = new TFile("Pythia_for_Jpsi.root");
	TH2F *hMassElecPosi_background = (TH2F*)file->Get("hMassElecPosi_background");
	TH1D *all = hMassElecPosi_background->ProjectionX("all", 1,50);
	
	const Int_t NBINS = 23;

	double border_arr[NBINS + 1] = {0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,25.,30.,50.};

	TH1D *hMassJpsi_params = new TH1D("hMassJpsi_params", "M_{J/#psi}, GeV/c^{2}", NBINS, border_arr);
	hMassJpsi_params->Sumw2();
	TH1D *hSigmaJpsi_params = new TH1D("hSigmaJpsi_params", "#sigma_{M_{J/#psi}}, GeV/c^{2}", NBINS, border_arr);
	hSigmaJpsi_params->Sumw2();

	TH1D *hMass[23];
	TCanvas *c1 = new TCanvas("c1","c1");
	c1->Divide(5,5);
	for (Int_t iPt=0; iPt<NBINS; iPt++) {
	  Int_t ptMin, ptMax;
	  ptMin = hMassElecPosi_background->GetYaxis()->FindBin(border_arr[iPt]+0.01);
	  ptMax = hMassElecPosi_background->GetYaxis()->FindBin(border_arr[iPt+1]-0.01);
	  TString name = Form("Jpsi_inv_mass_with_p_T_%02d",iPt);
	  TString title = Form("%.0f < p_T < %.0f",border_arr[iPt],border_arr[iPt+1]);
	  hMass[iPt] = hMassElecPosi_background->ProjectionX(name,ptMin,ptMax);
	  hMass[iPt]->SetTitle(title);
	  if (iPt>20) hMass[iPt]->Rebin(2);
	  if (iPt == 22) hMass[iPt]->Rebin(2);
	  c1->cd(iPt+1);
	  hMass[iPt]->Fit("gaus","","",2.8,3.4);
	  Double_t mMean = hMass[iPt]->GetFunction("gaus")->GetParameter(1);
	  Double_t sigma = hMass[iPt]->GetFunction("gaus")->GetParameter(2);
	  Double_t mError = hMass[iPt]->GetFunction("gaus")->GetParError(1);
	  Double_t sError = hMass[iPt]->GetFunction("gaus")->GetParError(2);
	  hMassJpsi_params->SetBinContent(iPt+1,mMean);
	  hMassJpsi_params->SetBinError  (iPt+1,mError);
	  hSigmaJpsi_params->SetBinContent(iPt+1,sigma);
          hSigmaJpsi_params->SetBinError  (iPt+1,sError);
	}


	c1->cd(24);
	hMassJpsi_params->Draw();

	c1->cd(25);
	hSigmaJpsi_params->Draw();

	c1->Print("Results.pdf"); 


	TFile *outFile = new TFile("Jpsi_decays_data_from_pythia.root", "RECREATE");
	
	for(int i = 0; i < NBINS; i++){
	  hMass[i]->Write();
	}
	
	hMassJpsi_params ->Write();
	hSigmaJpsi_params->Write();

	outFile->Close();
	delete outFile;
	return;
}
