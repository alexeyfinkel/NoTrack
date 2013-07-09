//This is my attempt at making a script to do Kevin-like smearing on Alex's ntuples.

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

#include "../src/ZEffTree.h"

//void smearElectron(TLorentzVector &) ;
double doSmearing(double, double, TFile*, TFile*, bool plot=false);//returns chi-squared of smeared mass hist

int SmearTry02(double mean_low=0, double mean_high=0.1, double sigma_low=0, double sigma_high=0.1, int Nstep=10)
{
	if(mean_low==0 && mean_low==0 && sigma_low==0 && sigma_high==0.1 && Nstep==10)
	{
		std::cout<<"WARNING! Smearing parameters not set. Using default!"<<std::endl;
	}
	
	//define some scanning variables; "Nstep" points in range for mean and sigma
	double meanStep = (mean_high-mean_low)/Nstep;
	double sigmaStep = (sigma_high-sigma_low)/Nstep;	
	double theMean = mean_low;
	double theSigma = (sigma_high+sigma_low)/2; //starting with middle-of-range sigma	
	double chi2_new, theChi2;
	
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
	//and make a file to store the results
	//TFile* file = new TFile("SmearedZ.root","recreate");
	
	//Make hists to plot chi-squared
	TH1D* chi2Mean = new TH1D("chi2Mean","Chi2 vs Mean",Nstep+2,mean_low,mean_high+meanStep);
	TH1D* chi2Sigma = new TH1D("chi2Sigma","Chi2 vs Sigma",Nstep+2,sigma_low,sigma_high+sigmaStep);
	
	//this is the smearing function!
	theChi2 = doSmearing(theMean, theSigma, f1, f2);
	chi2Mean->SetBinContent(1,theChi2);
	//scan for best Mean
	for(int i=1;i<(Nstep+1);i++)
	{
		chi2_new = doSmearing(mean_low+i*meanStep, theSigma, f1, f2);
		//meanArr[i] = mean_low+i*meanStep;
		//chi2Mean[i] = chi2_new;
		
		chi2Mean->SetBinContent(i+1,chi2_new);
		if(chi2_new < theChi2) 
		{
			theChi2 = chi2_new;
			theMean = mean_low+i*meanStep;
		}
	}
	//now scan for best sigma	
	for(int i=0;i<(Nstep+1);i++)
	{
		chi2_new = doSmearing(theMean, sigma_low+i*sigmaStep, f1, f2);
		chi2Sigma->SetBinContent(i+1,chi2_new);
		if(chi2_new < theChi2) 
		{
			theChi2 = chi2_new;
			theSigma = sigma_low+i*sigmaStep;
		}
	}
    
    std::cout<<"Best chi squared is "<<theChi2<<" per DOF for mean of "<<theMean<<" and sigma of "<<theSigma<<std::endl;
    
    //closure test
	/*for(int i=0;i<11;i++)
	{
		chi2_new = doSmearing(mean_low+i*meanStep, theSigma, f1, f2, c1);
		if(chi2_new < theChi2) 
		{
			theChi2 = chi2_new;
			theMean = mean_low+i*meanStep;
		}
	}
	
	std::cout<<"Closure Test (mean): New best chi squared is "<<theChi2<<" per DOF for mean of "<<theMean<<" and sigma of "<<theSigma<<std::endl;*/
	
	//fitting function:
	TF1* fun = new TF1("fun","[0]+[1]*pow( x-[2] ,2)",0,1);
	//fit vs. mean
	fun->SetParameter(0,1);
	fun->SetParLimits(0,0,1000);
	fun->SetParameter(1,1);
	//fun->SetParLimits(1,0,1000);
	fun->SetParameter(2,theMean);
	fun->SetParLimits(2,mean_low,mean_high);	
	chi2Mean->Draw();
	chi2Mean->Fit(fun,"","",mean_low,mean_high);
	c1->Update();
	c1->Print("Smearing_Chi2vsMean.png");
	
	theMean = fun->GetParameter(2);
	
	//now vs. sigma
	fun->SetParameter(0,1);
	fun->SetParLimits(0,0,1000);
	fun->SetParameter(1,1);
	fun->SetParLimits(1,1,1e6);
	fun->SetParameter(2,theSigma);
	fun->SetParLimits(2,sigma_low,sigma_high);
	chi2Sigma->Draw();
	chi2Sigma->Fit(fun,"","",sigma_low,sigma_high);
	c1->Update();
	c1->Print("Smearing_Chi2vsSigma.png");
	
	theSigma = fun->GetParameter(2);
	
	//make a graph of the best fit:
	theChi2 = doSmearing(theMean, theSigma, f1, f2, true);
	
	std::cout<<"Best best chi-squared from parabola fitting is "<<theChi2<<" at mean of "<<theMean<<" and sigma of "<<theSigma<<std::endl;
	
	return 0;
}







double doSmearing(double mean, double sigma, TFile* f1, TFile* f2, bool plot)
{
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
		if( (fabs(ze1->gen.eta[1])>2.5) && (fabs(ze1->gen.eta[1])<3.0) && (ze1->reco.isSelected(1,NT)) )
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
		    	
		if( (fabs(ze2->reco.eta[1])>2.5) && (fabs(ze2->reco.eta[1])<3.0) && (ze2->reco.isSelected(1,NT)) )
		{
			
			datM->Fill(ze2->reco.mz);
		}
		
		run = ze2->GetNextEvent();
	}

	//oldM->Scale( datM->Integral(43,52)/oldM->Integral(43,52) );
	newM->Scale( datM->Integral(43,52)/newM->Integral(43,52) );

	double Chi2=0;
	for(int i=43;i<53;i++) //WARNING: hard-coded bins!
	{
		Chi2 += pow((newM->GetBinContent(i) - datM->GetBinContent(i)),2)/pow(datM->GetBinError(i),2);
	}    
	Chi2 /= 10;
	
	if(plot)
	{
		TCanvas* c2 = new TCanvas("C1", "c1", 800, 600);
		c2->cd();
		
		newM->GetXaxis()->SetRangeUser(50, 150);
		newM->GetYaxis()->SetRangeUser(0, 1.3*newM->GetBinContent(newM->GetMaximumBin()));
		newM->SetLineWidth(3);
		newM->SetLineColor(4);
		newM->Draw("hist");

		datM->SetMarkerStyle(20);
		datM->Draw("same E");

		c2->Update();
		c2->Print("SmearingGaussOptimal.png");
		c2->Close();
	}
	
	newM->Delete();
	oldM->Delete();
	datM->Delete();
	
	return Chi2;
}
//void smearElectron(TLorentzVector &electron)
//{
//}
