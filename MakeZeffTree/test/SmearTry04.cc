//This is my attempt at making a script to do Kevin-like smearing on Alex's ntuples.
//Version 4 includes quart subdivision and uses MC reco, like it should have been doing all along!

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
std::vector<double> doSmearing(std::vector<double>, std::vector<double>, TFile*, TFile*, std::vector<TH1*>&, bool plot=false);//returns chi-squared vector of smeared mass hists

int SmearTry04(double mean_low=0, double mean_high=0.1, double sigma_low=0, double sigma_high=0.1, int Nstep=10)
{
	if(mean_low==0 && mean_low==0 && sigma_low==0 && sigma_high==0.1 && Nstep==10)
	{
		std::cout<<"WARNING! Smearing parameters not set. Using default!"<<std::endl;
	}
	
	//define some scanning variables; "Nstep" points in range for mean and sigma
	double meanStep = (mean_high-mean_low)/Nstep;
	double sigmaStep = (sigma_high-sigma_low)/Nstep;	
	std::vector<double> currentMean(11,mean_low), theMean(11,mean_low);
	std::vector<double> theSigma (11,(sigma_high+sigma_low)/2 ), currentSigma(11,(sigma_high+sigma_low)/2); //starting with middle-of-range sigma	
	std::vector<double> chi2_new(11,0), theChi2(11,0);
	
	gROOT->SetStyle("Plain");	
	//setTDRStyle();
	gStyle->SetErrorX(0);
	gStyle->SetOptStat(0);
	
	TCanvas* c1 = new TCanvas("C1", "c1", 1200, 900);
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
	//and a data ntuple for comparison
	TFile* f2 = new TFile("/local/cms/user/finkel/NoTrack/Ntuple/DE_ReReco_AlexTrgr/makezefftree_AlexTrgr/DE_ReReco_AlexTrgr_002.root");	
	if(f1 == NULL)
	{
		std::cout<<"Failed to open Data file. Exiting."<<std::endl;
		exit(1);
	}
	//and make a file to store the results
	TFile* file = new TFile("SmearedZ_01.root","recreate");
	
	std::cout<<"All Files Found! Proceeding to make data hists!"<<std::endl;
	
	//Make the DATA hists (they do not change as different smearings are tried)	
	TH1D* datMAll = new TH1D("datMAll","M from Data, all NT",75,0,150);
	TH1D* datMNTplus = new TH1D("datMNTplus","M from Data, NT+",75,0,150);
	TH1D* datMNTminus = new TH1D("datMNTminus","M from Data, NT-",75,0,150);
	TH1D* datMQ0 = new TH1D("datMQ0","M from Data, NT-, #phi (0,#pi/2)",75,0,150);
	TH1D* datMQ1 = new TH1D("datMQ1","M from Data, NT-, #phi (#pi/2,#pi)",75,0,150);
	TH1D* datMQ2 = new TH1D("datMQ2","M from Data, NT-, #phi (#pi,3#pi/2)",75,0,150);
	TH1D* datMQ3 = new TH1D("datMQ3","M from Data, NT-, #phi (3#pi/2,2#pi)",75,0,150);
	TH1D* datMQ4 = new TH1D("datMQ4","M from Data, NT+, #phi (0,#pi/2)",75,0,150);
	TH1D* datMQ5 = new TH1D("datMQ5","M from Data, NT+, #phi (#pi/2,#pi)",75,0,150);
	TH1D* datMQ6 = new TH1D("datMQ6","M from Data, NT+, #phi (#pi,3#pi/2)",75,0,150);
	TH1D* datMQ7 = new TH1D("datMQ7","M from Data, NT+, #phi (3#pi/2,2#pi)",75,0,150);
	//(Q0-Q3 are NT- quarts; Q4-Q7 are NT+ quarts; Q8 is all of NT-; Q9 is all of NT+; Q10 is all of NT)
	
	std::cout<<"Data Hists Made! Proceeding to make the data hists vector!"<<std::endl;
	
	std::vector<TH1*> datM;
	datM.push_back(datMQ0);
	datM.push_back(datMQ1);
	datM.push_back(datMQ2);
	datM.push_back(datMQ3);
	datM.push_back(datMQ4);
	datM.push_back(datMQ5);
	datM.push_back(datMQ6);
	datM.push_back(datMQ7);
	datM.push_back(datMAll);
	datM.push_back(datMNTplus);
	datM.push_back(datMNTminus);
	
	std::cout<<"Vector Made! Proceeding to fill hists!"<<std::endl;
	
	//open file for the Data ntuple
	ZEffTree *ze2;
	ze2 = new ZEffTree(*f2,false);
	const std::string NT = "NTLooseElectronId-EtaDet";
	
	//Fill the data histograms:
	int run = true;
	while(run) //fill the data hists
	{
			ze2->Entries();

		if( (fabs(ze2->reco.eta[1])>2.5) && (fabs(ze2->reco.eta[1])<3.0) && (ze2->reco.isSelected(1,NT)) //&& (ze2->reco.phi[1]>0) && (ze2->reco.phi[1]<0.182)
			) //Limited eta-phi range for trying to look at specific crystals
		{			
			datM[10]->Fill(ze2->reco.mz);
			
			if( ze2->reco.eta[1]>0 )
			{
				datM[9]->Fill(ze2->reco.mz);
				if(ze2->reco.phi[1]<-1.5708) datM[4]->Fill(ze2->reco.mz);
				else if (ze2->reco.phi[1]<0) datM[5]->Fill(ze2->reco.mz);
				else if (ze2->reco.phi[1]<1.5708) datM[6]->Fill(ze2->reco.mz);
				else if (ze2->reco.phi[1]<3.1416) datM[7]->Fill(ze2->reco.mz);
				
			}
			else
			{
				datM[8]->Fill(ze2->reco.mz);
				if(ze2->reco.phi[1]<-1.5708) datM[0]->Fill(ze2->reco.mz);
				else if (ze2->reco.phi[1]<0) datM[1]->Fill(ze2->reco.mz);
				else if (ze2->reco.phi[1]<1.5708) datM[2]->Fill(ze2->reco.mz);
				else if (ze2->reco.phi[1]<3.1416) datM[3]->Fill(ze2->reco.mz);
			}
		}
		
		run = ze2->GetNextEvent();
	}	
	
	std::cout<<"Data Hists Filled!"<<std::endl;
	
	//test for empty hists!
	for(int j=0; j<11; j++)
	{
		std::cout<<"Hist for Q-"<<j<<" has "<<datM[j]->Integral()<<" entries."<<std::endl;
		if( datM[j]->Integral()==0 )
		{
			
			return 1;
		}
	}
	
	std::cout<<"No Empty Hists Found! Proceeding to make the Chi-Squared hists!"<<std::endl;
	
	//Make hists to plot chi-squared for all the subranges
	TH1D* AllChi2Mean = new TH1D("AllChi2Mean","Chi2 vs Mean All NT",Nstep+2,mean_low,mean_high+meanStep);
	TH1D* AllChi2Sigma = new TH1D("AllChi2Sigma","Chi2 vs Sigma All NT",Nstep+2,sigma_low,sigma_high+sigmaStep);
	TH1D* NTPlusChi2Mean = new TH1D("NTPlusChi2Mean","Chi2 vs Mean NT-",Nstep+2,mean_low,mean_high+meanStep);
	TH1D* NTPlusChi2Sigma = new TH1D("NTPlusChi2Sigma","Chi2 vs Sigma NT-",Nstep+2,sigma_low,sigma_high+sigmaStep);
	TH1D* NTMinusChi2Mean = new TH1D("NTMinusChi2Mean","Chi2 vs Mean NT+",Nstep+2,mean_low,mean_high+meanStep);
	TH1D* NTMinusChi2Sigma = new TH1D("NTMinusChi2Sigma","Chi2 vs Sigma NT+",Nstep+2,sigma_low,sigma_high+sigmaStep);
	TH1D* Q0Chi2Mean = new TH1D("Q0Chi2Mean","Chi2 vs Mean Q0",Nstep+2,mean_low,mean_high+meanStep);
	TH1D* Q0Chi2Sigma = new TH1D("Q0Chi2Sigma","Chi2 vs Sigma Q0",Nstep+2,sigma_low,sigma_high+sigmaStep);
	TH1D* Q1Chi2Mean = new TH1D("Q1Chi2Mean","Chi2 vs Mean Q1",Nstep+2,mean_low,mean_high+meanStep);
	TH1D* Q1Chi2Sigma = new TH1D("Q1Chi2Sigma","Chi2 vs Sigma Q1",Nstep+2,sigma_low,sigma_high+sigmaStep);
	TH1D* Q2Chi2Mean = new TH1D("Q2Chi2Mean","Chi2 vs Mean Q2",Nstep+2,mean_low,mean_high+meanStep);
	TH1D* Q2Chi2Sigma = new TH1D("Q2Chi2Sigma","Chi2 vs Sigma Q2",Nstep+2,sigma_low,sigma_high+sigmaStep);
	TH1D* Q3Chi2Mean = new TH1D("Q3Chi2Mean","Chi2 vs Mean Q3",Nstep+2,mean_low,mean_high+meanStep);
	TH1D* Q3Chi2Sigma = new TH1D("Q3Chi2Sigma","Chi2 vs Sigma Q3",Nstep+2,sigma_low,sigma_high+sigmaStep);
	TH1D* Q4Chi2Mean = new TH1D("Q4Chi2Mean","Chi2 vs Mean Q4",Nstep+2,mean_low,mean_high+meanStep);
	TH1D* Q4Chi2Sigma = new TH1D("Q4Chi2Sigma","Chi2 vs Sigma Q4",Nstep+2,sigma_low,sigma_high+sigmaStep);
	TH1D* Q5Chi2Mean = new TH1D("Q5Chi2Mean","Chi2 vs Mean Q5",Nstep+2,mean_low,mean_high+meanStep);
	TH1D* Q5Chi2Sigma = new TH1D("Q5Chi2Sigma","Chi2 vs Sigma Q5",Nstep+2,sigma_low,sigma_high+sigmaStep);
	TH1D* Q6Chi2Mean = new TH1D("Q6Chi2Mean","Chi2 vs Mean Q6",Nstep+2,mean_low,mean_high+meanStep);
	TH1D* Q6Chi2Sigma = new TH1D("Q6Chi2Sigma","Chi2 vs Sigma Q6",Nstep+2,sigma_low,sigma_high+sigmaStep);
	TH1D* Q7Chi2Mean = new TH1D("Q7Chi2Mean","Chi2 vs Mean Q7",Nstep+2,mean_low,mean_high+meanStep);
	TH1D* Q7Chi2Sigma = new TH1D("Q7Chi2Sigma","Chi2 vs Sigma Q7",Nstep+2,sigma_low,sigma_high+sigmaStep);
	
	std::cout<<"Chi-Squared Hists Made! Proceeding to make the chi2 hists vectors!"<<std::endl;
	
	std::vector<TH1*> Chi2Mean;
	std::vector<TH1*> Chi2Sigma;
	
	Chi2Mean.push_back(Q0Chi2Mean);
	Chi2Mean.push_back(Q1Chi2Mean);
	Chi2Mean.push_back(Q2Chi2Mean);
	Chi2Mean.push_back(Q3Chi2Mean);
	Chi2Mean.push_back(Q4Chi2Mean);
	Chi2Mean.push_back(Q5Chi2Mean);
	Chi2Mean.push_back(Q6Chi2Mean);
	Chi2Mean.push_back(Q7Chi2Mean);
	Chi2Mean.push_back(NTMinusChi2Mean);
	Chi2Mean.push_back(NTPlusChi2Mean);
	Chi2Mean.push_back(AllChi2Mean);
	
	Chi2Sigma.push_back(Q0Chi2Sigma);
	Chi2Sigma.push_back(Q1Chi2Sigma);
	Chi2Sigma.push_back(Q2Chi2Sigma);
	Chi2Sigma.push_back(Q3Chi2Sigma);
	Chi2Sigma.push_back(Q4Chi2Sigma);
	Chi2Sigma.push_back(Q5Chi2Sigma);
	Chi2Sigma.push_back(Q6Chi2Sigma);
	Chi2Sigma.push_back(Q7Chi2Sigma);
	Chi2Sigma.push_back(NTMinusChi2Sigma);
	Chi2Sigma.push_back(NTPlusChi2Sigma);
	Chi2Sigma.push_back(AllChi2Sigma);
	
	std::cout<<"Chi2 Vectors Made! Proceeding to do Smearing!"<<std::endl;
	
	//scan for best Mean for each hist
	theChi2 = doSmearing(theMean, theSigma, f1, file, datM);
	for(int i=1;i<(Nstep+1);i++)
	{
		for(int j=0;j<(11);j++)	{currentMean[j] += meanStep;} //itterate the mean
		
		chi2_new = doSmearing(currentMean, currentSigma, f1, file, datM);	//compute chi-squareds			
		
		for(int j=0;j<(11);j++)
		{
			Chi2Mean[j]->SetBinContent(i+1,chi2_new[j]); //enter the chi-squareds into their respective hists
			if(chi2_new[j] < theChi2[j]) 
			{
				theChi2[j] = chi2_new[j];
				theMean[j] = currentMean[j];
			}
		}
	}
	
	//now scan for best sigma	
	for(int i=0;i<(Nstep+1);i++)
	{
		for(int j=0;j<(11);j++)	{currentSigma[j] = sigma_low+i*sigmaStep;}
		
		chi2_new = doSmearing(theMean, currentSigma, f1, file, datM);
		
		for(int j=0;j<(11);j++)
		{
			Chi2Sigma[j]->SetBinContent(i+1,chi2_new[j]);
			if(chi2_new[j] < theChi2[j]) 
			{
				theChi2[j] = chi2_new[j];
				theSigma[j] = currentSigma[j];
			}
		}
	}
	
	
	//fitting function:
	TF1* fun = new TF1("fun","[0]+[1]*pow( x-[2] ,2)",0,1);
	
	for(int j=0; j<11; j++)
	{
		//TCanvas* c1 = new TCanvas("C1", "c1", 800, 600);
		//c1->cd();
		
		//fit vs. mean
		char fileNameMean[124];
		sprintf(fileNameMean, "SmearingByQuart/Chi2vsMean_Q%d.png",j);
		fun->SetParameter(0,1);
		fun->SetParLimits(0,0,1000);
		fun->SetParameter(1,1);
		//fun->SetParLimits(1,0,1000);
		fun->SetParameter(2,theMean[j]);
		fun->SetParLimits(2,mean_low,mean_high);	
		Chi2Mean[j]->Fit(fun,"","",mean_low,mean_high);
		c1->Update();
		c1->Print(fileNameMean);
		
		theMean[j] = fun->GetParameter(2);
		
		//now vs. sigma
		char fileNameSigma[124];
		sprintf(fileNameSigma, "SmearingByQuart/Chi2vsSigma_Q%d.png",j);
		fun->SetParameter(0,1);
		fun->SetParLimits(0,0,1000);
		fun->SetParameter(1,1);
		fun->SetParLimits(1,1,1e6);
		fun->SetParameter(2,theSigma[j]);
		fun->SetParLimits(2,sigma_low,sigma_high);
		Chi2Sigma[j]->Fit(fun,"","",sigma_low,sigma_high);
		c1->Update();
		c1->Print(fileNameSigma);
		c1->Clear();
		
		theSigma[j] = fun->GetParameter(2);
		
		std::cout<<"Best parabolic chi-squared for Quart-"<<j<<" is "<<theChi2[j]<<" at mean of "<<theMean[j]<<" and sigma of "<<theSigma[j]<<std::endl;
		//make a graph of the best fit:
	}
	theChi2 = doSmearing(theMean, theSigma, f1, file, datM, true);	
	
	//Poster smearing diagram:
	int randSeed=123456; //using constant seed for reproducibility
	TRandom3* rand = new TRandom3(randSeed);

	TLorentzVector elec1, elec2, theZ;
	float pt, eta, phi, E;

	TH1D* oldM = new TH1D("oldM","Unsmeared M",75,0,150);
	TH1D* newM = new TH1D("newM","Smeared M",75,0,150);
	TH1D* theDatM = (TH1D*)datM[9]->Clone();
	
	//make a tree that accesses the MC ntuple
	ZEffTree *ze1;
	ze1 = new ZEffTree(*f1,false);
	run = true;
	while (run) //fill MC hist and do smearing on the events
	{    	
		ze1->Entries();
		//NOTE! Selecting anything that's WP80 and is from NT region
		if( (fabs(ze1->reco.eta[1])>2.5) && (fabs(ze1->reco.eta[1])<3.0) 
			)//NT not imposed and range not limited for MC, since it does not know about crystal badness
		{
			pt = ze1->reco.pt[0];    	
			eta = ze1->reco.eta[0];
			phi = ze1->reco.phi[0];
			E = pt*cosh(eta);
			elec1.SetPtEtaPhiE(pt,eta,phi,E);
			
			pt = ze1->reco.pt[1]; //believe it or not, HERE is the smearing part:
			pt *= ( 1+ rand->Gaus(theMean[10],theSigma[10]) );//note: using PROPORTIONAL correction now    	
			eta = ze1->reco.eta[1];
			phi = ze1->reco.phi[1];    
			E = pt*cosh(eta);	
			elec2.SetPtEtaPhiE(pt,eta,phi,E);
			
			theZ = elec1+elec2;
			
			oldM->Fill(ze1->reco.mz);
			newM->Fill(theZ.M());
		}
		
		run = ze1->GetNextEvent();
	}
	
	oldM->SetTitle(";M_{ee} (GeV);Events (Scaled)");
	oldM->SetLineWidth(3);
	oldM->SetLineColor(2);
	oldM->GetXaxis()->SetRangeUser(50,130);
	oldM->GetYaxis()->SetRangeUser(0, 1.3*oldM->GetBinContent(oldM->GetMaximumBin()) );
	newM->Scale( oldM->GetBinContent(oldM->GetMaximumBin())/newM->GetBinContent(newM->GetMaximumBin()) );
	newM->SetLineWidth(3);
	newM->SetLineColor(4);
	theDatM->Scale(oldM->GetBinContent(oldM->GetMaximumBin())/theDatM->GetBinContent(theDatM->GetMaximumBin()));
	theDatM->SetMarkerStyle(20);
	oldM->Draw("HIST");
	newM->Draw("SAME HIST");
	theDatM->Draw("SAME P");

	l1->AddEntry(oldM,"Unsmeared");	
	l1->AddEntry(newM,"Smeared");
	l1->AddEntry(theDatM,"Data");
	l1->Draw();
	
	c1->Print("SmearingByQuart/SmearPoster.png");
	c1->Clear();
	
	file-> Write();
	
	return 0;
}



std::vector<double> doSmearing(std::vector<double> mean, std::vector<double> sigma, TFile* f1, TFile* file, std::vector<TH1*>& datM, bool plot)
{
	int randSeed=123456; //using constant seed for reproducibility
	TRandom3* rand = new TRandom3(randSeed);

	TLorentzVector elec1, elec2, theZ;
	float pt, eta, phi, E;

	const std::string NT = "NTLooseElectronId-EtaDet";
	std::vector<double> Chi2(11,0);
	
	for(int j=0; j<11; j++)
	{
		//make pre- and post-smearing hists
		TH1D* oldM = new TH1D("oldM","Unsmeared M",75,0,150);
		TH1D* newM = new TH1D("newM","Smeared M",75,0,150);
		
		//make a tree that accesses the MC ntuple
		ZEffTree *ze1;
		ze1 = new ZEffTree(*f1,false);
		bool run = true;
		while (run) //fill MC hist and do smearing on the events
		{    	
			ze1->Entries();
			//NOTE! Selecting anything that's WP80 and is from NT region
			if( (fabs(ze1->reco.eta[1])>2.5) && (fabs(ze1->reco.eta[1])<3.0) 
				)//NT not imposed and range not limited for MC, since it does not know about crystal badness
			{
				pt = ze1->reco.pt[0];    	
				eta = ze1->reco.eta[0];
				phi = ze1->reco.phi[0];
				E = pt*cosh(eta);
				elec1.SetPtEtaPhiE(pt,eta,phi,E);
				
				pt = ze1->reco.pt[1]; //believe it or not, HERE is the smearing part:
				pt *= ( 1+ rand->Gaus(mean[j],sigma[j]) );//note: using PROPORTIONAL correction now    	
				eta = ze1->reco.eta[1];
				phi = ze1->reco.phi[1];    
				E = pt*cosh(eta);	
				elec2.SetPtEtaPhiE(pt,eta,phi,E);
				
				theZ = elec1+elec2;
				
				oldM->Fill(ze1->reco.mz);
				newM->Fill(theZ.M());
			}
			
			run = ze1->GetNextEvent();
		}
		
		int nDat, nMC;

		nDat = datM[j]->Integral();
		nMC = newM->Integral();

		newM->Scale( datM[j]->Integral(43,52)/newM->Integral(43,52) );

		for(int i=43;i<53;i++) //WARNING: hard-coded bins!
		{
			Chi2[j] += pow((newM->GetBinContent(i) - datM[j]->GetBinContent(i)),2)/pow(datM[j]->GetBinError(i),2);
		}    
		Chi2[j] /= 10;

		if(plot)
		{
			TCanvas* c2 = new TCanvas("C1", "c2", 1200, 900);
			c2->cd();
			
			TLegend* l1 = new TLegend(0.65,0.7,0.9,0.9);
			l1->SetFillColor(10);			
			
			char dataLeg[128];
			char mcLeg[128];
			char fileName[128];
			char title[128];
			sprintf(dataLeg,"Data, %d events", nDat);
			sprintf(mcLeg,"MC, %d events", nMC);
			sprintf(fileName,"SmearingByQuart/SmearingGaussOptimal_Q%d.png",j);
			sprintf(title,"Smeared M for Q-%d; #mu=%f, #sigma=%f, #chi^{2}/ndof=%f",j,mean[j],sigma[j],Chi2[j]);
			l1->AddEntry(datM[j], dataLeg);
			l1->AddEntry(newM, mcLeg);
			
			newM->GetXaxis()->SetRangeUser(50, 150);
			newM->GetYaxis()->SetRangeUser(0, 1.3*newM->GetBinContent(newM->GetMaximumBin()));
			newM->SetLineWidth(3);
			newM->SetLineColor(4);
			newM->SetTitle(title);
			newM->Draw("hist");
			l1->Draw();

			datM[j]->SetMarkerStyle(20);
			datM[j]->Draw("same E");
			//datM[j]->Write();

			c2->Update();
			c2->Print(fileName);
			c2->Clear();				
		}
		
		oldM->Delete();
		newM->Delete();
		
		file->Write();
		
	}
	
	return Chi2;
}
