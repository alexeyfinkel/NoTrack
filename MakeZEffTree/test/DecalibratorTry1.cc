//this is a copy of the MC projected cal consts script, only ran on data.

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iomanip>
#include <utility> //pairs live here
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <math.h>
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TColor.h"
#include "TProfile.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TGraph.h"

#include "../src/ZEffTree.h"



Double_t voigt(Double_t *x, Double_t *par)
{
    return (par[0]*(TMath::Voigt((x[0]-par[1]),par[2],par[3])) + par[4] + par[5]*(x[0]-80));
}

int decalibrateData()
{
    gROOT->SetStyle("Plain");	
	gStyle->SetErrorX(0);
	gStyle->SetOptStat(0);
	
	TCanvas* c1 = new TCanvas("C1", "c1", 1200, 800);
	c1->cd();
	TLegend* l1 = new TLegend(0.65,0.7,0.9,0.9);
	l1->SetFillColor(10);
	
	//prepare hists:
	TH2D *oldCalConstsN = new TH2D("oldCalConstsN","2012-06-ConstsN",40,30,70,40,30,70);
	TH2D *oldCalConstsP = new TH2D("oldCalConstsP","2012-06-ConstsP",40,30,70,40,30,70);
	TH1D *originalMZ = new TH1D("originalMZ","Original Z Mass;M_{ee}, GeV",30,60,120);
	TH1D *decalibratedMZ = new TH1D("decalibratedMZ","Decalibrated Z Mass;M_{ee}, GeV",30,60,120);

    char line[1000];
    int cix, ciy, iz;
    float calConst, sigma;
    ifstream oldConstsFile;
    //read in cumulative consts file
    oldConstsFile.open("/afs/cern.ch/work/a/afinkel/public/NoTrack/IC_FILES/EcalIntercalibConstants_V20120620_piZPhiSEtaScale2012_IOV2_AlphaStudies.txt");
	if(!oldConstsFile.is_open())
	{
		std::cout<<"Failed to open cumul. cal. constants file. Existing"<<std::endl;
		return 1;
	}
    do
    {
        oldConstsFile.getline(line, 1000)>>cix>>ciy>>iz>>calConst>>sigma;
        if(iz==0)continue;
        if(  (((cix-50.5)*(cix-50.5)+(ciy-50.5)*(ciy-50.5))  < 132.25) || (((cix-50.5)*(cix-50.5)+(ciy-50.5)*(ciy-50.5)) > 342.25)  ) continue;
        if(iz==-1) 
        {
            oldCalConstsN->SetBinContent(cix-30,ciy-30,calConst);
        }
        else if(iz==1)
        {
            oldCalConstsP->SetBinContent(cix-30,ciy-30,calConst);
        }
        else 
        {
            std::cout<<"Error reading in a cal. const. Exiting"<<std::endl;
            return 1;
        }
        
    } while(!oldConstsFile.eof());    
			
	//grab a data ntuple
	TFile* f2 = new TFile("/afs/cern.ch/work/a/afinkel/public/NoTrack/Ntuples/Data2012/DE_ReReco_2012Full_WithClusters_001.root");	
	if(f2 == NULL)
	{
		std::cout<<"Failed to open Data file. Exiting."<<std::endl;
		return 1;
	}
    
	//make data ntuple:
	ZEffTree* ze2 = new ZEffTree(*f2,false);
	//Fill the data histograms:
	int ix,iy;
	int hix, hiy;
	double correction, observedPt;
	
	TLorentzVector elec1, elec2, theZ;
	float pt, eta, phi, E;
    
	for(int event=0; event<ze2->Entries(); event++ ) //fill the data hists
	{
		ix=ze2->reco.ix[1];//seed ix
		iy=ze2->reco.iy[1];//seed iy
        
		if( (fabs(ze2->reco.eta[1])>2.5) && (fabs(ze2->reco.eta[1])<3.0) 
			 &&(ix>29) &&(ix<71) &&(iy>29) &&(iy<71)
			 && ze2->reco.isSelected(1,"NTLooseElectronId-EtaDet") //using only events that pass selection now!
		  ) 
		{	
			if( (((ix-50.5)*(ix-50.5)+(iy-50.5)*(iy-50.5))  < 132.25) || (((ix-50.5)*(ix-50.5)+(iy-50.5)*(iy-50.5)) > 342.25) )
			{
				ze2->GetNextEvent();
				continue;
			}

			//compute the correction factor for observed energy
			correction=0;
            for(unsigned int k=0; k<(ze2->ixs->size()); k++ )
            {	
                hix = ze2->ixs->at(k);
                hiy = ze2->iys->at(k);
                if( (((hix-50.5)*(hix-50.5)+(hiy-50.5)*(hiy-50.5))  < 132.25) || (((hix-50.5)*(hix-50.5)+(hiy-50.5)*(hiy-50.5)) > 342.25) ) 
                //outside NT; assume adequatelycalibrated already
                {
                    correction += ze2->hitEnergyFractions->at(k);
                    continue;
                }
                if(ze2->reco.eta[1]>0) correction += ze2->hitEnergyFractions->at(k)/oldCalConstsP->GetBinContent(hix-30,hiy-30);
                else correction += ze2->hitEnergyFractions->at(k)/oldCalConstsN->GetBinContent(hix-30,hiy-30);
            }
			observedPt = ze2->reco.pt[1]*correction;
						
			//recalculate Z mass using found correction
			pt = ze2->reco.pt[0];    	//tracked electron
			eta = ze2->reco.eta[0];
			phi = ze2->reco.phi[0];
			E = pt*cosh(eta);
			elec1.SetPtEtaPhiE(pt,eta,phi,E);				
			eta = ze2->reco.eta[1];
			phi = ze2->reco.phi[1];    
			E = observedPt*cosh(eta);	
			elec2.SetPtEtaPhiE(observedPt,eta,phi,E);
			theZ = elec1+elec2;		//new Z vector             
			
			//store new and old Z mass:
			originalMZ->Fill(ze2->reco.mz);  
			decalibratedMZ->Fill(theZ.M());
			//std::cout<<theZ.M()<<std::endl;
	    }
		ze2->GetNextEvent();
	}//end loop over events in ntuple

    //Make the plots!!!
    originalMZ->SetMarkerSize(2);
    originalMZ->SetMarkerStyle(20);
    originalMZ->SetMarkerColor(kBlue+1);
    originalMZ->Draw("P");
    decalibratedMZ->SetMarkerSize(2);
    decalibratedMZ->SetMarkerStyle(20);
    decalibratedMZ->SetMarkerColor(kRed+2);
    decalibratedMZ->Draw("same P");
    l1->AddEntry(originalMZ,"Original","p");
    l1->AddEntry(decalibratedMZ,"Decalibrated","p");
    l1->Draw();
    c1->Print("DecalibrationTest_01.png");
    
    c1->Close();
    return 0;
}












