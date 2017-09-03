#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"

void draw() {
		TFile* f1 = new TFile("/moose/Bes3User/lliu/workarea/Dsigma/09mc/ana/root/sigma.root","read");
		TTree *t1 =  (TTree*)f1->Get("FIT");
		TCanvas *c1 = new TCanvas("c1", "sigma mass", 800, 600);
		TCanvas *c2 = new TCanvas("c2", "nbar energy", 800, 600);
		TCanvas *c3 = new TCanvas("c3", "nbar hits", 800, 600);
		TCanvas *c4 = new TCanvas("c4", "nbar second moment", 800, 600);
		TCanvas *c5 = new TCanvas("c5", "nbar cos", 800, 600);
		TH1D *mc09_h1 = new TH1D("mc09_h1", "sigma mass", 100, 1.1, 1.4);
		TH1D *mc09_h2 = new TH1D("mc09_h2", "sigmabar mass", 100, 1.1, 1.4);
		TH1D *mc09_h3 = new TH1D("mc09_h3", "nbar energy", 100, 0.0, 2.0);
		TH1D *mc09_h4 = new TH1D("mc09_h4", "nbar hits", 50, 0, 150);
		TH1D *mc09_h5 = new TH1D("mc09_h5", "nbar second moment", 100, 0, 200);
		TH1D *mc09_h6 = new TH1D("mc09_h6", "nbar cos", 100, -1, 1);



		double chisq, chisqv2, sigma_m, sigmabar_m, nbar_e;
		Double_t nbar_hits, nbar_secmom, nbar_cos;


		t1->SetBranchAddress("Sigmac_m", &sigma_m);
		t1->SetBranchAddress("aSigmac_m", &sigmabar_m);
		t1->SetBranchAddress("nbar_energy", &nbar_e);
		t1->SetBranchAddress("nbar_hits", &nbar_hits);
		t1->SetBranchAddress("nbar_secmom", &nbar_secmom);
		t1->SetBranchAddress("nbar_cos", &nbar_cos);
		for (int i; i<t1->GetEntries(); i++){
				t1->GetEntry(i);
				mc09_h1->Fill(sigma_m);
				mc09_h2->Fill(sigmabar_m);
				mc09_h3->Fill(nbar_e);
				mc09_h4->Fill(nbar_hits);
				mc09_h5->Fill(nbar_secmom);
				mc09_h6->Fill(nbar_cos);
		}
		mc09_h1->SetMarkerStyle(7);
		mc09_h1->SetMarkerSize(0.01);
		mc09_h1->SetMarkerColor(1);
		mc09_h1->GetXaxis()->SetTitle("M_{#Sigma^{-}}(Gev/c^{2})");
		mc09_h1->GetYaxis()->SetTitle("Events/(3.0 MeV/c^{2})");
		mc09_h1->GetYaxis()->SetTitleOffset(1.4);
		mc09_h2->SetMarkerStyle(7);
		mc09_h2->SetMarkerSize(0.01);
		mc09_h2->SetMarkerColor(1);
		mc09_h2->GetXaxis()->SetTitle("M_{#bar{#Sigma}^{-}}(Gev/c^{2})");
		mc09_h2->GetYaxis()->SetTitle("Events/(3.0 MeV/c^{2})");
		mc09_h2->GetYaxis()->SetTitleOffset(1.4);







		gStyle->SetPadTickX(1);
		gStyle->SetPadTickY(1);
		gStyle->SetOptTitle(1); 
		gStyle->SetOptStat(1110); 
		gStyle->SetOptFit(0); 
		gStyle->SetFitFormat("6.3g"); 
		gStyle->SetPalette(1); 
			   
		gStyle->SetTextFont(42);
		gStyle->SetTextSize(0.05);
		gStyle->SetTitleFont(42,"xyz");
		gStyle->SetTitleSize(0.05);
		gStyle->SetLabelFont(42,"xyz");
		gStyle->SetLabelSize(0.05);
		gStyle->SetTitleXOffset(0.8);
		gStyle->SetTitleYOffset(1.1);
		gROOT->ForceStyle();


		c1->cd();
		mc09_h1->Draw("E1");
		mc09_h2->Draw("same");
		c2->cd();
		mc09_h3->Draw();
		c3->cd();
		mc09_h4->Draw();
		c4->cd();
		mc09_h5->Draw();
		c5->cd();
		mc09_h6->Draw();


}

