//this is a script for optimizing shower-shape cuts by making Z-mass plots in sections of each variable

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
//#include "TGraph.h"
//#include <TLorentzVector.h>
//#include "TRandom3.h"
#include <math.h>
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
//#include <TMath.h>
//#include "TGraphAsymmErrors.h"
//#include "TGraphPainter.h"

#include "../src/ZEffTree.h"


int SSVSections()
{
	gROOT->SetStyle("Plain");	
	gStyle->SetErrorX(0);
	gStyle->SetOptStat(0);
	
	TCanvas* c1 = new TCanvas("C1", "c1", 1200, 900);
	c1->cd();
	
	TLegend* l1 = new TLegend(0.65,0.7,0.9,0.9);
	l1->SetFillColor(10);
	
	//grab a data ntuple
	TFile* f2 = new TFile("/local/cms/user/finkel/NoTrack/Ntuple/DE_ReReco_NewerCuts/DE_ReReco2012FULL_Jan13_AlexTrgr_newerCuts_001.root");	
	if(f2 == NULL)
	{
		std::cout<<"Failed to open Data file. Exiting."<<std::endl;
		exit(1);
	}
	
	//make variable segmentation
	double lowBound[5], highBound[5], step[5]; //order is R9, sIeIe, H/EM, hIso, eIso
	std::vector<int> nSlices(5,10);//10 slices in each var. by default
	std::vector<TH1D*> uncutVars;
	
	lowBound[0] = 0.85;highBound[0]=1.05;
	lowBound[1] = 0.005;highBound[1]=0.035;
	lowBound[2] = 0;highBound[2]=0.05;
	lowBound[3] = 0;highBound[3]=0.15;
	lowBound[4] = 0;highBound[4]=0.05;
	for(int i=0; i<5;i++)
	{
		step[i]=(highBound[i]-lowBound[i])/nSlices[i];
	}	
	
	//make hist collections:
	std::vector<TH1D*> MvsVar[5];
	for(int var =0; var<5; var++)
	{
		char name[100], title[200];
		switch (var)
		{
			case 0:
				sprintf(name,"rNine");
				sprintf(title,"M vs. R9");
				break;
			case 1:
				sprintf(name,"sIeIe");
				sprintf(title,"M vs. #sigma_{i#eta,i#eta}");
				break;
			case 2:
				sprintf(name,"HoEM");
				sprintf(title,"M vs. H/EM");
				break;
			case 3:
				sprintf(name,"hIso");
				sprintf(title,"M vs. HcalIso");
				break;
			case 4:
				sprintf(name,"eIso");
				sprintf(title,"M vs. EcalIso");
				break;
		}
		uncutVars.push_back( new TH1D(name,name,50,lowBound[var],highBound[var]) );
		
		for(int i=0; i<nSlices[var]; i++)
		{			
			char varName[100], varTitle[200]; 
			sprintf(varName,"%s%d",name,i+1);
			sprintf(varTitle,"%s, [%f,%f];M_{ee}; Events/2GeV",title,lowBound[var]+i*step[var], lowBound[var]+(i+1)*step[var]);

			MvsVar[var].push_back( new TH1D(varName,varTitle,50,50,150) );
		}
	}
	TH1D* allM = new TH1D("allM","All M_{ee}",50,50,150);
	
	//make the ntuple:
	ZEffTree* ze2 = new ZEffTree(*f2,false);
	
	double vars[5];
	//Fill the data histograms:
	for(int event=0; event<ze2->Entries(); event++ ) //fill the data hists
	{
		if( (fabs(ze2->reco.eta[1])>2.5) && (fabs(ze2->reco.eta[1])<3.0) ) 
		{	//fill uncut variable hists
			vars[0] = ze2->reco.RNine[1];
			vars[1] = ze2->reco.Sieie[1];
			vars[2] = ze2->reco.HoEM[1];
			vars[3] = ze2->reco.Hiso[1];
			vars[4] = ze2->reco.Eiso[1];
			//fill mass slice hists
			allM->Fill(ze2->reco.mz);
			for(int var=0; var<5; var++)
			{
				uncutVars[var]->Fill(vars[var]);
				for(int i=0; i<nSlices[var]; i++)
				{				
					if( ( vars[var] > lowBound[var]+i*step[var] ) && ( vars[var] < lowBound[var]+(i+1)*step[var] ) )
					{
						MvsVar[var][i]->Fill(ze2->reco.mz);
						break;
					}
				}
			}
			//fill n-1 cut hists -- coming soon!
		}
		ze2->GetNextEvent();
	}
	
	for(int var=0; var<5; var++)
	{
		char filename[100];
		if(var==0) sprintf(filename,"Efficiencies/ShowerShapeVars/RnineAll.png");
		if(var==1) sprintf(filename,"Efficiencies/ShowerShapeVars/SieieAll.png");
		if(var==2) sprintf(filename,"Efficiencies/ShowerShapeVars/HoEMAll.png");
		if(var==3) sprintf(filename,"Efficiencies/ShowerShapeVars/HisoAll.png");
		if(var==4) sprintf(filename,"Efficiencies/ShowerShapeVars/EisoAll.png");
		uncutVars[var]->SetLineColor(4);
		uncutVars[var]->SetLineWidth(4);
		uncutVars[var]->Draw("hist");
		c1->Print(filename);
		c1->Clear();
		for(int i=0; i<nSlices[var]; i++)
		{
			char filename2[200];
			if(var==0) sprintf(filename2,"Efficiencies/ShowerShapeVars/MvsRnine%d.png",i);
			if(var==1) sprintf(filename2,"Efficiencies/ShowerShapeVars/MvsSieie%d.png",i);
			if(var==2) sprintf(filename2,"Efficiencies/ShowerShapeVars/MvsHoEM%d.png",i);
			if(var==3) sprintf(filename2,"Efficiencies/ShowerShapeVars/MvsHiso%d.png",i);
			if(var==4) sprintf(filename2,"Efficiencies/ShowerShapeVars/MvsEiso%d.png",i);
			MvsVar[var][i]->SetLineColor(4);
			MvsVar[var][i]->SetLineWidth(4);
			MvsVar[var][i]->Draw("hist");
			c1->Print(filename2);
			c1->Clear();
		}
	}
	allM->SetLineColor(4);
	allM->SetLineWidth(4);
	allM->Draw("hist");
	c1->Print("Efficiencies/ShowerShapeVars/AllM.png");
	c1->Close();
	
	return 0;
}
