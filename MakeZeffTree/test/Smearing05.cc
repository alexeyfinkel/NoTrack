//This is my attempt at making a script to do Kevin-like smearing on Alex's ntuples.
//Version 5 includes pt and eta hists, makes better plots, and should eventually have higher subdivision in phi.

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
#include "TLatex.h"
#include "TColor.h"

#include "../src/ZEffTree.h"

const int nSegments = 19; //WARNING: GLOBAL VARIABLE! Must be 2^n + 3, n<=4 for symmetry reasons

std::vector<double> doSmearing(std::vector<double>, std::vector<double>, TFile*, TFile*, std::vector<TH1*>&, bool plot=false);//returns chi-squared vector of smeared mass hists

int Smearing05(double mean_low=0, double mean_high=0.1, double sigma_low=0, double sigma_high=0.1, int Nstep=10)
{
	if(mean_low==0 && mean_low==0 && sigma_low==0 && sigma_high==0.1 && Nstep==10)
	{
		std::cout<<"WARNING! Smearing parameters not set. Using default!"<<std::endl;
	}
	
	//define some scanning variables; "Nstep" points in range for mean and sigma
	double meanStep = (mean_high-mean_low)/Nstep;
	double sigmaStep = (sigma_high-sigma_low)/Nstep;	
	std::vector<double> currentMean(nSegments,mean_low), theMean(nSegments,mean_low), theMeanError(nSegments,0);
	std::vector<double> theSigma (nSegments,(sigma_high+sigma_low)/2 ), currentSigma(nSegments,(sigma_high+sigma_low)/2), theSigmaError(nSegments,0); //starting with middle-of-range sigma	
	std::vector<double> chi2_new(nSegments,0), theChi2(nSegments,0);
	
	gROOT->SetStyle("Plain");	
	//setTDRStyle();
	gStyle->SetErrorX(0);
	gStyle->SetOptStat(0);
	
	TCanvas* c1 = new TCanvas("C1", "c1", 1200, 900);
	c1->cd();
	
	TLegend* l1 = new TLegend(0.7,0.75,0.9,0.9);
	l1->SetFillColor(10);

	//get an ntuple to smear
	TFile* f1 = new TFile("/local/cms/user/finkel/NoTrack/Ntuple/DE_MC2012_AlexTrgr3/DE_MC2012_AlexTrgr_newerCuts.root");	
	if(f1 == NULL)
	{
		std::cout<<"Failed to open MC file. Exiting."<<std::endl;
		exit(1);
	}
	//and a data ntuple for comparison
	TFile* f2 = new TFile("/local/cms/user/finkel/NoTrack/Ntuple/DE_ReReco_NewerCuts/DE_ReReco2012FULL_Jan13_AlexTrgr_newerCuts_001.root");	
	if(f1 == NULL)
	{
		std::cout<<"Failed to open Data file. Exiting."<<std::endl;
		exit(1);
	}
	//and make a file to store the results
	TFile* file = new TFile("SmearedZ_01.root","recreate");
	
	//Make vectors of needed hists
	std::vector<TH1*> datM, datPt, datEta;
	
	std::cout<<"All Files Found! Proceeding to make data hists!"<<std::endl;
	for(int i=0; i<nSegments; i++)
	{
		char titleM[128],titlePt[128],titleEta[128];
		char nameM[128],namePt[128],nameEta[128];
		sprintf(titleM, "M_{ee}, Segment-%d", i);
		sprintf(titlePt, "NoTrack p_{T}, Segment-%d", i);
		sprintf(titleEta, "NoTrack |#eta|, Segment-%d", i);
		sprintf(nameM, "datM%d", i);
		sprintf(namePt, "datPt%d", i);
		sprintf(nameEta, "datEta%d", i);
		
		datM.push_back( new TH1D(nameM,titleM,75,0,150) );
		datPt.push_back( new TH1D(namePt,titlePt,50,0,100) );
		datEta.push_back( new TH1D(nameEta,titleEta,25,2.5,3) );
	}
	
	//open file for the Data ntuple
	ZEffTree *ze2;
	ze2 = new ZEffTree(*f2,false);
	const std::string NT = "NTLooseElectronId-EtaDet";
	
	//Fill the data histograms:
	int run = true;
	int seg;
	while(run) //fill the data hists
	{
			ze2->Entries();

		if( (fabs(ze2->reco.eta[1])>2.5) && (fabs(ze2->reco.eta[1])<3.0) && (ze2->reco.isSelected(1,NT)) 
			) 
		{			
			datM[nSegments-1]->Fill(ze2->reco.mz);//All NT
			datPt[nSegments-1]->Fill(ze2->reco.pt[1]);
			datEta[nSegments-1]->Fill( fabs(ze2->reco.eta[1]) );
			
			if( ze2->reco.eta[1]>0 )
			{
				datM[nSegments-2]->Fill(ze2->reco.mz);//NT+
				datPt[nSegments-2]->Fill(ze2->reco.pt[1]);
				datEta[nSegments-2]->Fill( fabs(ze2->reco.eta[1]) );
				
				seg = (int)( (ze2->reco.phi[1]+3.14159)*(nSegments-3)/4/3.14159+(nSegments-3)/2 );//note that phi originally in [-Pi,Pi] range!
				datM[seg]->Fill(ze2->reco.mz);
				datPt[seg]->Fill(ze2->reco.pt[1]);
				datEta[seg]->Fill( fabs(ze2->reco.eta[1]) );
			}
			else	//eta1<0
			{
				datM[nSegments-3]->Fill(ze2->reco.mz);//NT-
				datPt[nSegments-3]->Fill(ze2->reco.pt[1]);
				datEta[nSegments-3]->Fill( fabs(ze2->reco.eta[1]) );
				
				seg = (int)( (ze2->reco.phi[1]+3.14159)*(nSegments-3)/4/3.14159 );
				datM[seg]->Fill(ze2->reco.mz);
				datPt[seg]->Fill(ze2->reco.pt[1]);
				datEta[seg]->Fill( fabs(ze2->reco.eta[1]) );
			}
		}
		
		run = ze2->GetNextEvent();
	}//close while loop
	
	std::cout<<"Data Hists Filled!"<<std::endl;
	
	std::vector<TH1*> Chi2Mean, Chi2Sigma;
	for(int i=0; i<nSegments;i++)
	{
		char chi2MeanName[128], chi2MeanTitle[128];
		char chi2SigmaName[128], chi2SigmaTitle[128];
		sprintf(chi2MeanName,"Chi2Mean%d",i);
		sprintf(chi2MeanTitle,"#chi^{2} vs. Mean, Segment %d;#mu;#chi^{2}/ndof",i);
		sprintf(chi2SigmaName,"Chi2Sigma%d",i);
		sprintf(chi2SigmaTitle,"#chi^{2} vs. Sigma, Segment %d;#sigma;#chi^{2}/ndof",i);
		Chi2Mean.push_back( new TH1D(chi2MeanName,chi2MeanTitle,Nstep+2,mean_low,mean_high+meanStep) );
		Chi2Sigma.push_back( new TH1D(chi2SigmaName,chi2SigmaTitle,Nstep+2,sigma_low,sigma_high+sigmaStep) );
	}

	
	//scan for best Mean for each hist
	std::cout<<"Scanning for best Mean."<<std::endl;
	theChi2 = doSmearing(theMean, theSigma, f1, file, datM);
	for(int i=1;i<(Nstep+1);i++)
	{
		for(int j=0;j<nSegments;j++)	{currentMean[j] += meanStep;} //itterate the mean
		
		chi2_new = doSmearing(currentMean, currentSigma, f1, file, datM);	//compute chi-squareds			
		
		for(int j=0;j<nSegments;j++)
		{
			Chi2Mean[j]->SetBinContent(i+1,chi2_new[j]); //enter the chi-squareds into their respective hists
			if(chi2_new[j] < theChi2[j]) 
			{
				theChi2[j] = chi2_new[j];
				theMean[j] = currentMean[j];
			}
		}
		std::cout<<"Yow-"<<i<<"!"<<std::endl;
	}
	
	//now scan for best sigma	
	std::cout<<"Scanning for best Sigma."<<std::endl;
	for(int i=0;i<(Nstep+1);i++)
	{
		for(int j=0;j<nSegments;j++)	{currentSigma[j] = sigma_low+i*sigmaStep;}
		
		chi2_new = doSmearing(theMean, currentSigma, f1, file, datM);
		
		for(int j=0;j<nSegments;j++)
		{
			Chi2Sigma[j]->SetBinContent(i+1,chi2_new[j]);
			if(chi2_new[j] < theChi2[j]) 
			{
				theChi2[j] = chi2_new[j];
				theSigma[j] = currentSigma[j];
			}
		}
		std::cout<<"Ping-"<<i<<"!"<<std::endl;
	}
	
	
	//fitting function:
	TF1* fun = new TF1("fun","[0]+[1]*pow( x-[2] ,2)",0,1);
	
	for(int j=0; j<nSegments; j++)
	{
		//fit vs. mean
		char fileNameMean[124];
		sprintf(fileNameMean, "SmearingByEighth/Chi2VsMean/Chi2vsMean_Q%d.png",j);
		fun->SetParameter(0,1);
		fun->SetParLimits(0,0,1000);
		fun->SetParameter(1,1);
		fun->SetParLimits(1,0,1e9);
		fun->SetParameter(2,theMean[j]);
		fun->SetParLimits(2,-0.5,1.5*mean_high);	
		Chi2Mean[j]->Fit(fun,"","",mean_low,mean_high);
		c1->Update();
		c1->Print(fileNameMean);
		
		theMean[j] = fun->GetParameter(2);
		theMeanError[j] = fun->GetParError(2);
		
		//now vs. sigma
		char fileNameSigma[124];
		sprintf(fileNameSigma, "SmearingByEighth/Chi2VsSigma/Chi2vsSigma_Q%d.png",j);
		fun->SetParameter(0,1);
		fun->SetParLimits(0,0,1000);
		fun->SetParameter(1,1);
		fun->SetParLimits(1,1,1e9);
		fun->SetParameter(2,theSigma[j]);
		fun->SetParLimits(2,0,1.5*sigma_high);
		Chi2Sigma[j]->Fit(fun,"","",sigma_low,sigma_high);
		//c1->Update();
		c1->Print(fileNameSigma);
		c1->Clear();
		
		theSigma[j] = fun->GetParameter(2);		
		theSigmaError[j] = fun->GetParError(2);		
		//std::cout<<"Best parabolic chi-squared for Quart-"<<j<<" is "<<theChi2[j]<<" at mean of "<<theMean[j]<<" and sigma of "<<theSigma[j]<<std::endl;
		//make a graph of the best fit:
	}
	theChi2 = doSmearing(theMean, theSigma, f1, file, datM, true);	
	c1->Close();
	
//====================================================================

	//Poster smearing diagrams:
	TCanvas* c3 = new TCanvas("C3", "c3", 1200, 900);
	c3->cd();
	c3->SetMargin(0.12,0.02,0.1,0.02);
	TLegend* l2 = new TLegend(0.75,0.8,0.98,0.98);
	l2->SetFillColor(10);
	
	int randSeed=123456; //using constant seed for reproducibility
	TRandom3* rand = new TRandom3(randSeed);

	TLorentzVector elec1, elec2, theZ;
	float pt, eta, phi, E;

	TH1D* oldM = new TH1D("oldM","Unsmeared M",75,0,150);
	TH1D* newM = new TH1D("newM","Smeared M",75,0,150);
	TH1D* oldPt = new TH1D("oldPt","Unsmeared Pt",50,0,100);
	TH1D* newPt = new TH1D("newPt","Smeared Pt",50,0,100);
	TH1D* Eta = new TH1D("eta","NoTrack #eta",25,2.5,3);
	TH1D* theDatM = (TH1D*)datM[nSegments-1]->Clone();
	TH1D* theDatPt = (TH1D*)datPt[nSegments-1]->Clone();
	TH1D* theDatEta = (TH1D*)datEta[nSegments-1]->Clone();
	
	TLatex title;
	title.SetTextSize(0.04);
	title.SetTextFont(42);
	title.SetNDC(true);
	
	//make a tree that accesses the MC ntuple
	ZEffTree *ze1;
	ze1 = new ZEffTree(*f1,false);
	run = true;
	while (run) //fill MC hist and do smearing on the events
	{    	
		ze1->Entries();
		//NOTE! Selecting anything that's WP80 and is from NT region
		if( (fabs(ze1->reco.eta[1])>2.5) && (fabs(ze1->reco.eta[1])<3.0) && (ze1->reco.isSelected(1,NT))
			)//range not limited for MC, since it does not know about crystal badness
		{
			pt = ze1->reco.pt[0];    	
			eta = ze1->reco.eta[0];
			phi = ze1->reco.phi[0];
			E = pt*cosh(eta);
			elec1.SetPtEtaPhiE(pt,eta,phi,E);
			
			pt = ze1->reco.pt[1]; //believe it or not, HERE is the smearing part:
			pt *= ( 1+ rand->Gaus(theMean[nSegments-1],theSigma[nSegments-1]) );//note: using PROPORTIONAL correction now    	
			eta = ze1->reco.eta[1];
			phi = ze1->reco.phi[1];    
			E = pt*cosh(eta);	
			elec2.SetPtEtaPhiE(pt,eta,phi,E);
			
			theZ = elec1+elec2;
			
			oldM->Fill(ze1->reco.mz);
			newM->Fill(theZ.M());
			oldPt->Fill(ze1->reco.pt[1]);
			newPt->Fill(pt);
			Eta->Fill(eta);
		}
		
		run = ze1->GetNextEvent();
	}
	
	theDatM->SetTitle(";M_{ee} (GeV);Events/2GeV");
	theDatM->GetXaxis()->SetRangeUser(60,130);
	theDatM->GetYaxis()->SetRangeUser(0, 1.2*theDatM->GetBinContent(theDatM->GetMaximumBin()) );
	theDatM->GetYaxis()->SetTitleOffset(1.5);
	theDatM->SetMarkerStyle(20);
	newM->Scale(theDatM->Integral(42,52)/newM->Integral(42,52));
	newM->SetLineWidth(4);
	newM->SetLineColor(4);
	newM->SetTitle("");
	theDatM->Draw("E");
	newM->Draw("SAME HIST");
	
	title.DrawLatex(0.13,0.9,"CMS internal, preliminary");
	l2->AddEntry(newM,"MC");
	l2->AddEntry(theDatM,"Data");
	l2->Draw();
	c3->Print("SmearingByEighth/SmearPosterM.pdf");
	c3->Print("SmearingByEighth/SmearPosterM.png");
	c3->Clear();
	l2->Clear();
	
	theDatPt->SetTitle(";NoTrack p_{T} (GeV);Events/2GeV");
	theDatPt->GetXaxis()->SetRangeUser(15,80);
	theDatPt->GetYaxis()->SetRangeUser(0, 1.2*theDatPt->GetBinContent(theDatPt->GetMaximumBin()) );
	theDatPt->SetMarkerStyle(20);
	theDatPt->GetYaxis()->SetTitleOffset(1.5);
	newPt->Scale(theDatPt->Integral()/newPt->Integral());
	newPt->SetLineWidth(4);
	newPt->SetLineColor(4);
	theDatPt->Draw("E");
	newPt->Draw("SAME HIST");

	title.DrawLatex(0.13,0.9,"CMS internal, preliminary");	
	l2->AddEntry(newPt,"MC");
	l2->AddEntry(theDatPt,"Data");
	l2->Draw();
	c3->Print("SmearingByEighth/SmearPosterPt.pdf");
	c3->Print("SmearingByEighth/SmearPosterPt.png");
	c3->Clear();
	l2->Clear();
	
	theDatEta->SetTitle(";NoTrack #eta;Events");
	theDatEta->GetYaxis()->SetRangeUser(0, 1.3*theDatEta->GetBinContent(theDatEta->GetMaximumBin()) );
	theDatEta->SetMarkerStyle(20);
	theDatEta->GetYaxis()->SetTitleOffset(1.5);
	Eta->Scale(theDatEta->Integral()/Eta->Integral());
	Eta->SetLineWidth(4);
	Eta->SetLineColor(4);
	theDatEta->Draw("E");
	Eta->Draw("SAME HIST");

	title.DrawLatex(0.13,0.9,"CMS internal, preliminary");
	l2->AddEntry(Eta,"MC");	
	l2->AddEntry(theDatEta,"Data");
	l2->Draw();
	c3->Print("SmearingByEighth/SmearPosterEta.pdf");
	c3->Print("SmearingByEighth/SmearPosterEta.png");
	l2->Clear();
	c3->Clear();
	
	//make mean and sigma vs. segment plots
	TH1D* MeanBySegment = new TH1D("MeanBySegment","Smearing #mu vs. Segment",nSegments,1,nSegments+1);
	TH1D* SigmaBySegment = new TH1D("SigmaBySegment","Smearing #sigma vs. Segment",nSegments,1,nSegments+1);
	
	for(int i=1; i<=nSegments; i++)
	{
		MeanBySegment->SetBinContent(i,theMean[i-1]);
		MeanBySegment->SetBinError(i,theMeanError[i-1]);
		SigmaBySegment->SetBinContent(i,theSigma[i-1]);
		SigmaBySegment->SetBinError(i,theSigmaError[i-1]);
	}
	MeanBySegment->GetXaxis()->SetBinLabel(nSegments,"All NT");
	MeanBySegment->GetXaxis()->SetBinLabel(nSegments-1,"NT+");
	MeanBySegment->GetXaxis()->SetBinLabel(nSegments-2,"NT-");
	SigmaBySegment->GetXaxis()->SetBinLabel(nSegments,"All NT");
	SigmaBySegment->GetXaxis()->SetBinLabel(nSegments-1,"NT+");
	SigmaBySegment->GetXaxis()->SetBinLabel(nSegments-2,"NT-");
	MeanBySegment->SetTitle(";Segment;#mu");
	SigmaBySegment->SetTitle(";Segment;#sigma");
	MeanBySegment->SetMarkerStyle(20);
	SigmaBySegment->SetMarkerStyle(20);
	MeanBySegment->SetMarkerColor(kBlue+1);
	SigmaBySegment->SetMarkerColor(kBlue+1);
	MeanBySegment->SetMarkerSize(1);
	SigmaBySegment->SetMarkerSize(1);
	//MeanBySegment->SetLineWidth(2);
	//SigmaBySegment->SetLineWidth(2);
	MeanBySegment->GetYaxis()->SetRangeUser(0, 1.2*MeanBySegment->GetBinContent(MeanBySegment->GetMaximumBin()) );
	SigmaBySegment->GetYaxis()->SetRangeUser(0,1.2*SigmaBySegment->GetBinContent(SigmaBySegment->GetMaximumBin()));
	MeanBySegment->Draw("E1");
	c3->Print("SmearingByEighth/MeanBySegment.png");
	c3->Clear();
	SigmaBySegment->Draw("E1");
	c3->Print("SmearingByEighth/SigmaBySegment.png");
	c3->Close();
	
	
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
	std::vector<double> Chi2(nSegments,0);
	
	for(int j=0; j<nSegments; j++)
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
			if( (fabs(ze1->reco.eta[1])>2.5) && (fabs(ze1->reco.eta[1])<3.0) && (ze1->reco.isSelected(1,NT))
				)//range not limited for MC, since it does not know about crystal badness
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

		nDat = (int)(datM[j]->Integral());
		nMC = (int)(newM->Integral());

		newM->Scale( datM[j]->Integral(42,53)/newM->Integral(42,53) );

		for(int i=42;i<53;i++) //WARNING: hard-coded bins!
		{
			Chi2[j] += pow((newM->GetBinContent(i) - datM[j]->GetBinContent(i)),2)/pow(datM[j]->GetBinError(i),2);
		}    
		Chi2[j] /= 10;
			
		if(plot)
		{		
			TCanvas* c2 = new TCanvas("C2", "c2", 1200, 900);
			c2->cd();	
			TLegend* l1 = new TLegend(0.65,0.7,0.9,0.9);
			l1->SetFillColor(10);			
			
			char dataLeg[128];
			char mcLeg[128];
			char fileName[128];
			char title[128];
			sprintf(dataLeg,"Data, %d events", nDat);
			sprintf(mcLeg,"MC, %d events", nMC);
			sprintf(fileName,"SmearingByEighth/SmearedM/SmearingGaussOptimal_Q%d.png",j);
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

			//c2->Update();
			c2->Print(fileName);	
			c2->Close();			
		}
		oldM->Delete();
		newM->Delete();
		
		file->Write();
		
	}
	
	return Chi2;
}
