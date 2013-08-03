//This is my attempt at making a script to do Kevin-like smearing on Alex's ntuples.

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include "TRandom3.h"
#include <math.h>
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"

#include "../src/ZEffTree.h"

//void smearElectron(TLorentzVector &) ;

int SmearTry01(float mean=0, float sigma=0.01)
{
	if(mean==0 && sigma==0.01)
	{
		std::cout<<"WARNING! Smearing parameters not set. Using default (minimal smearing)!"<<std::endl;
	}
	
	gROOT->SetStyle("Plain");	
	//setTDRStyle();
	gStyle->SetErrorX(0);
	gStyle->SetOptStat(0);
	
	TCanvas* c1 = new TCanvas("C1", "c1", 800, 600);
	c1->cd();

	TLegend* l1 = new TLegend(0.65,0.7,0.9,0.9);
	l1->SetFillColor(10);

	//get an ntuple to smear
	TFile* f1 = new TFile("/local/cms/user/finkel/NoTrack/Ntuple/DE_MC2012B/makezefftree/DE_MC2012B_001.root");	
	if(f1 == NULL)
	{
		std::cout<<"Failed to open MC file. Exiting."<<std::endl;
		exit(1);
	}
	//and a data ntuple for comparison
	TFile* f2 = new TFile("/local/cms/user/finkel/NoTrack/Ntuple/DE_Data2012B/makezefftree/DE_Data2012B_001.root");	
	if(f1 == NULL)
	{
		std::cout<<"Failed to open Data file. Exiting."<<std::endl;
		exit(1);
	}
	
	TFile* file = new TFile("SmearedZ.root","recreate");
	
	//make a tree that accesses the MC ntuple
	ZEffTree *ze1;
    ze1 = new ZEffTree(*f1,false);
    
    //and one for the Data ntuple
    ZEffTree *ze2;
    ze2 = new ZEffTree(*f2,false);
    
    //some hists to make for comparison
    TH1* oldM = new TH1F("oldM","Unsmeared M",75,0,150);
    TH1* newM = new TH1F("newM","Smeared M",75,0,150);
    TH1* datM = new TH1F("datM","M from Data",75,0,150);
    
    int randSeed=123456; //using constant seed for reproducibility
    TRandom3* rand = new TRandom3(randSeed);
    
    TLorentzVector elec1, elec2, theZ;
    float pt, eta, phi, E;
    
    const std::string NT = "NTLooseElectronId-EtaDet";
    bool run = true;
    while (run)
    {    	
    	ze1->Entries();
    	//NOTE! Selecting anything that's WP80 and is from NT
    	if( (fabs(ze1->gen.eta[1])>2.5) && (fabs(ze1->gen.eta[1])<3.0) )//&& (ze1->reco.isSelected(1,NT)) )
    	{
			pt = ze1->gen.pt[0];    	
			eta = ze1->gen.eta[0];
			phi = ze1->gen.phi[0];
			E = pt*cosh(eta);
			elec1.SetPtEtaPhiE(pt,eta,phi,E);
			
			pt = ze1->gen.pt[1];
			pt *= ( 1+ rand->Gaus(mean,sigma) );//note: using PROPORTIONAL correction now    	
			eta = ze1->gen.eta[1];
			phi = ze1->gen.phi[1];    
			E = pt*cosh(eta);	
			elec2.SetPtEtaPhiE(pt,eta,phi,E);
			
			theZ = elec1+elec2;
			
			//std::cout<<"Unsmeared M = "<<ze1->reco.mz<<", \"smeared\" M = "<<theZ.M()<<std::endl;		
			oldM->Fill(ze1->gen.mz);
			newM->Fill(theZ.M());
		}
		
    	run = ze1->GetNextEvent();
    }
    
    run = true;
    while(run)
    {
    	ze2->Entries();
    	
    	//std::cout<<"Ping!"<<std::endl;
    	    	
    	if( (fabs(ze2->reco.eta[1])>2.5) && (fabs(ze2->reco.eta[1])<3.0) )//&& (ze2->reco.isSelected(1,NT)) )
    	{
    		
    		datM->Fill(ze2->reco.mz);
    	}
    	
    	run = ze2->GetNextEvent();
    }
    
    //compute smearing chi-squared
    
    oldM->SetLineWidth(5);
    oldM->SetLineColor(2);
    oldM->Scale( datM->Integral(43,52)/oldM->Integral(43,52) );
    
    newM->Scale( datM->Integral(43,52)/newM->Integral(43,52) );
    
    oldM->GetXaxis()->SetRangeUser(50, 150);
    oldM->GetYaxis()->SetRangeUser(0, 1.2*oldM->GetBinContent(oldM->GetMaximumBin()));
    oldM->Draw("hist");
    
    newM->SetLineWidth(3);
    newM->SetLineColor(4);    
    newM->Draw("same hist");
    
    datM->SetMarkerStyle(20);
    datM->Draw("same E");
    
    l1->AddEntry(oldM, "Unsmeared MC");
    l1->AddEntry(newM, "Smeared MC");
    l1->AddEntry(datM, "MZ from Data");
    l1->Draw();
    
    c1->Update();
    c1->Print("SmearingTryGauss.png");
	
	file->Write();
	
	
    double Chi2=0;
    for(int i=43;i<53;i++) //WARNING: hard-coded bins!
    {
    	Chi2 += pow((newM->GetBinContent(i) - datM->GetBinContent(i)),2)/pow(datM->GetBinError(i),2);
    }    
    Chi2 /= 10;
    
    std::cout<<"Chi squared is "<<Chi2<<" per DOF."<<std::endl;
	
	return 0;
}


//void smearElectron(TLorentzVector &electron)
//{
//}
