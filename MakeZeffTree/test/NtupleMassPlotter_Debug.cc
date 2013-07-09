//This will be a test of whether I can use the Ntuple that Alexe's code made
//to make efficiency plots

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <TFile.h>
#include <TH1.h>
#include "TROOT.h"
#include "TCanvas.h"
//#include "TStyle.h"
#include "TLegend.h"

#include "../src/ZEffTree.h"
//#include "tdrstyle.C"


int NtupleMassPlotterD()
{
	gROOT->SetStyle("Plain");	
	//setTDRStyle();
	
	TCanvas* c1 = new TCanvas("C1", "c1", 800, 600);
	
	//TLegend* l1 = new TLegend(0.75,0.8,1,1);
	//l1->SetFillColor(10);
	
	TFile* f1 = new TFile("/local/cms/user/finkel/NoTrack/Ntuple/DE_MC2012B/makezefftree/DE_MC2012B_001.root");	
	if(f1 == NULL)
	{
		std::cout<<"Failed to open MC file. Exiting."<<std::endl;
		exit(1);
	}	
	//TFile* f2 = new TFile("temp.root","new");
	
	TH1* h1 = new TH1F("M_allNT", "EBEB-debug", 50, 0, 5 );
	
	//ZEffTree *ze1;
    //ze1 = new ZEffTree(*f1,false);
    
    h1->FillRandom("gaus", 10000);
    
    //run over first ntuple
    /*bool run = true;
    while (run)
    {
        ze1->Entries();
        
        h1->Fill( ze1->reco.mz );
        
    	run = ze1->GetNextEvent();
    }*/    
    
    
    
    //h1->ResetAttMarker();
    h1->SetMarkerStyle(21);
    h1->SetFillColor(kBlue);
    h1->SetFillStyle(3001);
    h1->SetMarkerColor(4);
    //h1->SetMarkerSize(1);
    h1->Draw("HIST");
    c1->Update();
	
	return 0;
}























