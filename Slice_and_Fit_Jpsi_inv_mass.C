#include <math.h> 
#include <iostream>
TF1 *Mass_fit_func = new TF1("Mass_fit_func", "[0]*x*x*x + [1]*x*x + [2]*x +[3]", 0., 50.);
TF1 *Sigma_fit_func = new TF1("Sigma_fit_func", "[0]*exp([1]*(x - [2])) + [3]*x +[3]", 0., 50.);


void Slice_and_Fit_Jpsi_inv_mass(){

	TFile *file = new TFile("Pythia_for_Jpsi_10pow8_ev.root");
	TH2F *hMassElecPosi_background = (TH2F*)file->Get("hMassElecPosi_background");
	TH1D *all = hMassElecPosi_background->ProjectionX("all", 1,50);
	
	const Int_t NBINS = 19;
	const Int_t NNBINS = 38;

	double border_arr[NBINS + 1] = {0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,12.,14.,16.,18.,20.,25.,30.,35.,50.};
	double fit_arr[NNBINS + 1] = {3.05,3.15,3.04,3.16,3.04,3.16,3.03,3.17,3.02,3.18,3.02,3.18, 
				      3.02,3.18,3.01,3.19,3.01,3.19,2.99,3.20,2.97,3.21,2.96,3.23,
				      2.95,3.25,2.93,3.27,2.90,3.30,2.87,3.31,2.80,3.40,2.77,3.41,
				      2.75,3.43};

	TH1D *hMassJpsi_params = new TH1D("hMassJpsi_params", "M_{J/#psi}, GeV/c^{2}", NBINS, border_arr);
	hMassJpsi_params->Sumw2();
	TH1D *hSigmaJpsi_params = new TH1D("hSigmaJpsi_params", "#sigma_{M_{J/#psi}}, GeV/c^{2}", NBINS, border_arr);
	hSigmaJpsi_params->Sumw2();
	

	TF1 *Mass_fit_func = new TF1("Mass_fit_func", "[0]*x*x*x + [1]*x*x + [2]*x +[3]", 0., 50.);
	TF1 *Sigma_fit_func = new TF1("Sigma_fit_func", "[0]*exp([1]*(x - [2])) + [3]*x +[4]", 0., 50.);

	TH1D *hMass[NBINS];
	TCanvas *c1 = new TCanvas("c1","c1");
	c1->Divide(6,4);
	for (Int_t iPt=0; iPt<NBINS; iPt++) {
	  Int_t ptMin, ptMax;
	  ptMin = hMassElecPosi_background->GetYaxis()->FindBin(border_arr[iPt]+0.01);
	  ptMax = hMassElecPosi_background->GetYaxis()->FindBin(border_arr[iPt+1]-0.01);
	  TString name = Form("Jpsi_inv_mass_with_p_T_%02d",iPt);
	  TString title = Form("%.0f < p_T < %.0f",border_arr[iPt],border_arr[iPt+1]);
	  hMass[iPt] = hMassElecPosi_background->ProjectionX(name,ptMin,ptMax);
	  hMass[iPt]->SetTitle(title);
	  //if (iPt > 20) hMass[iPt]->Rebin(2);
	  if (iPt == NBINS - 2) hMass[iPt]->Rebin(2);
	  if (iPt == NBINS - 1) hMass[iPt]->Rebin(3);
	  c1->cd(iPt+1);
	  hMass[iPt]->Fit("gaus","","",fit_arr[2*iPt],fit_arr[2*iPt + 1]);
	  gStyle->SetOptFit(111);
	  Double_t mMean = hMass[iPt]->GetFunction("gaus")->GetParameter(1);
	  Double_t sigma = hMass[iPt]->GetFunction("gaus")->GetParameter(2);
	  Double_t mError = hMass[iPt]->GetFunction("gaus")->GetParError(1);
	  Double_t sError = hMass[iPt]->GetFunction("gaus")->GetParError(2);
	  hMassJpsi_params->SetBinContent(iPt+1,mMean);
	  hMassJpsi_params->SetBinError  (iPt+1,mError);
	  hSigmaJpsi_params->SetBinContent(iPt+1,sigma);
          hSigmaJpsi_params->SetBinError  (iPt+1,sError);
	}


	c1->cd(NBINS + 1);
	hMassJpsi_params->Draw();
	hMassJpsi_params->Fit("Mass_fit_func", "", "", 0., 50.);
	gStyle->SetOptFit(111);


	c1->cd(NBINS + 2);
	hSigmaJpsi_params->Draw();
	hSigmaJpsi_params->Fit("Sigma_fit_func", "", "", 0., 50.);
	gStyle->SetOptFit(111);

	c1->Print("Results_corr_2.pdf"); 


	TFile *outFile = new TFile("Jpsi_decays_data_from_pythia_corr_2.root", "RECREATE");
	
	for(int i = 0; i < NBINS; i++){
	  hMass[i]->Write();
	}
	
	hMassJpsi_params ->Write();
	hSigmaJpsi_params->Write();

	outFile->Close();
	return;
}
