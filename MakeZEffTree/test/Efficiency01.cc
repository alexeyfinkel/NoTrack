// this is an attempt at getting the NoTrack selection efficiency
//(currently using old sample from Alex's Single Electron skim)

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include "TGraph.h"
#include <TLorentzVector.h>
#include "TRandom3.h"
#include <math.h>
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include <TMath.h>

#include "../src/ZEffTree.h"

struct fitStruct
{	
	TH1* datHist;
	fitStruct(TH1*, TH1*);
	TF1* fitFunct;
	void doFit();
	double getSignalEvt();
	
	private:
	struct fitFunction
	{
		TH1* sigHist;
		
	    fitFunction(TH1* sh)
	    {
	    	sigHist = (TH1*)sh->Clone();
	        sigHist->Scale(1/sigHist->Integral());
	    }
		Double_t fun(Double_t *x,Double_t *par);
		double operator()(double *x, double *par)
		{
			return fun(x, par);
		}
	} fit;
};

fitStruct::fitStruct(TH1* signal, TH1* data) : fit(signal)
{
	datHist = (TH1*)data->Clone();
	//TF1* fitFunct = new TF1("fitFunct", "fun", 0, 150);
}

Double_t fitStruct::fitFunction::fun(Double_t *x,Double_t *par)
{
	//return (TMath::Erfc((par[0]-x[0])/par[1]) * ( par[2]*exp(-x[0]/par[3]) ) + par[4]*sigHist->GetBinContent(sigHist->FindBin(x[0])) );
	return exp(par[0]-x[0]*par[1]) + par[2]*sigHist->GetBinContent(sigHist->FindBin(x[0]));
}

void fitStruct::doFit()
{
	fitFunct = new TF1("fitFunct", fit, 50, 150, 3);
	fitFunct->SetParLimits(0,0,1e6);
	fitFunct->SetParLimits(1,0,100);
	fitFunct->SetParLimits(2,0,1e6);
	fitFunct->SetParameter(0,5);
	fitFunct->SetParameter(1,-0.5);
	fitFunct->SetParameter(2,50);
	//fitFunct->SetParLimits(3,0,100);
	//fitFunct->SetParLimits(4,0,1e6);

	datHist->GetXaxis()->SetRangeUser(50, 150);	
	datHist->Fit("fitFunct","RML");
}

double fitStruct::getSignalEvt()
{
	return fitFunct->GetParameter(2);
}



int Efficiency01(double mean=0.059, double sigma=0.091)
{
	if(mean==0.059 && sigma==0.093)
	{
		std::cout<<"NOTE: Using All-NT mean and sigma for signal smearing!"<<std::endl;
	}
	
	gROOT->SetStyle("Plain");	
	gStyle->SetErrorX(0);
	gStyle->SetOptStat(0);
	
	TCanvas* c1 = new TCanvas("C1", "c1", 800, 600);
	c1->cd();
	
	TLegend* l1 = new TLegend(0.65,0.7,0.9,0.9);
	l1->SetFillColor(10);

	//get an ntuple to smear
	TFile* f1 = new TFile("/local/cms/user/finkel/NoTrack/Ntuple/DE_MC2012B/WP80andNT/makezefftree/DE_MC_2012B_WP80NT_001.root");	
	if(f1 == NULL)
	{
		std::cout<<"Failed to open MC file. Exiting."<<std::endl;
		exit(1);
	}
	//and a data ntuple to fit
	TFile* f2 = new TFile("/local/cms/user/finkel/NoTrack/Ntuple/DE_ReReco_AlexTrgr/makezefftree_AlexTrgr/DE_ReReco_AlexTrgr_002.root");	
	if(f1 == NULL)
	{
		std::cout<<"Failed to open Data file. Exiting."<<std::endl;
		exit(1);
	}
	
	//make relevant hists
	TH1D* smearedM = new TH1D("smearedM","Signal, all-NT smearing",75,0,150);
	TH1D* datMbefore = new TH1D("datMbefore","MZ from Data, no NT cut",75,0,150);
	TH1D* datMafter = new TH1D("datMafter","MZ from Data, with NT cut",75,0,150);
	
	std::cout<<"Ready to fill hists!"<<std::endl;
	
	//open file for the Data ntuple
	ZEffTree *ze2, *ze1;
	ze1 = new ZEffTree(*f1,false);
	ze2 = new ZEffTree(*f2,false);
	const std::string NT = "NTLooseElectronId-EtaDet";
	
	//Fill the data histograms:
	bool run = true;
	while(run) //fill the data hists
	{
			ze2->Entries();

		if( (fabs(ze2->reco.eta[1])>2.5) && (fabs(ze2->reco.eta[1])<3.0) ) 
		{			
			datMbefore->Fill(ze2->reco.mz);
			if( ze2->reco.isSelected(1,NT) )
			{
				datMafter->Fill(ze2->reco.mz);
			}
		}
		
		run = ze2->GetNextEvent();
	}
	
	TLorentzVector elec1, elec2, theZ;
	float pt, eta, phi, E;
	int randSeed=123456; //using constant seed for reproducibility
	TRandom3* rand = new TRandom3(randSeed);
	
	run = true;
	while (run) //fill MC hist and do smearing on the events
	{    	
		ze1->Entries();
		//NOTE! Selecting anything that's WP80 and is from NT region
		if( (fabs(ze1->gen.eta[1])>2.5) && (fabs(ze1->gen.eta[1])<3.0) 
			)//NT not imposed and range not limited for MC, since it does not know about crystal badness
		{
			pt = ze1->gen.pt[0];    	
			eta = ze1->gen.eta[0];
			phi = ze1->gen.phi[0];
			E = pt*cosh(eta);
			elec1.SetPtEtaPhiE(pt,eta,phi,E);
			
			pt = ze1->gen.pt[1]; //believe it or not, HERE is the smearing part:
			pt *= ( 1 + rand->Gaus(mean,sigma) );//note: using PROPORTIONAL correction now    	
			eta = ze1->gen.eta[1];
			phi = ze1->gen.phi[1];    
			E = pt*cosh(eta);	
			elec2.SetPtEtaPhiE(pt,eta,phi,E);
			
			theZ = elec1+elec2;
			
			smearedM->Fill(theZ.M());
		}
		
		run = ze1->GetNextEvent();
	}
	
	std::cout<<"Hists made, starting the fitting process."<<std::endl;
	
	fitStruct before(smearedM,datMbefore);
	before.doFit();
	double sigEvtsBefore = before.getSignalEvt();
	c1->Print("Efficiency_Before.png");
	c1->Clear();
	
	fitStruct after(smearedM,datMafter);
	after.doFit();
	double sigEvtAfter = after.getSignalEvt();
	c1->Print("Efficiency_after.png");
	c1->Close();
	
	std::cout<<"There are "<<sigEvtsBefore<<" signal events before the NT cut and "<<sigEvtAfter<<" events after.\n"<<std::endl;
	std::cout<<"The resulting efficiency is "<<(sigEvtAfter/sigEvtsBefore*100)<<"%."<<std::endl;
	
	return 0;
}

























