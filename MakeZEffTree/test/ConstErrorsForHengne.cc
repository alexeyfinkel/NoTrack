//this will make the fit uncertainty / value of my cal. consts for Hengne

#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <TH1.h>
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"

int plot()
{
    gROOT->SetStyle("Plain");	
	gStyle->SetErrorX(0);
	gStyle->SetOptStat(1111);
	
	TCanvas* c1 = new TCanvas("C1", "c1", 1200, 900);
	c1->cd();
    
    TH1D *relUncertsAll = new TH1D("relUncertsAll","Relative Fit Uncertainty;#deltaC/C, %;crystals",100,0,5);
    TH1D *relUncertsP = new TH1D("relUncertsP","Relative Fit Uncertainty, NT+;#deltaC/C, %;crystals",100,0,5);
    TH1D *relUncertsN = new TH1D("relUncertsN","Relative Fit Uncertainty, NT-;#deltaC/C, %;crystals",100,0,5);
    
    char line[1000];
    int cix, ciy, iz;
    float calConst, sigma;
    ifstream constsFile;
    //read in ODDS cal. consts file (read into "EVENS" hists!):
	constsFile.open("Calibration/CalConsts_Data_Full.txt");
	if(!constsFile.is_open())
	{
		std::cout<<"Failed to open cal. constants file. Existing"<<std::endl;
		return 1;
	}
    
    do
    {
        constsFile.getline(line, 1000)>>cix>>ciy>>iz>>calConst>>sigma;
        relUncertsAll->Fill(sigma/calConst*100);
        if(iz==1) relUncertsP->Fill(sigma/calConst*100);
        if(iz==-1) relUncertsN->Fill(sigma/calConst*100);
        
    } while(!constsFile.eof());
    
    relUncertsAll->SetMarkerStyle(20);
    relUncertsAll->Draw("P");
    //c1->SetLogy(true);
    c1->Print("Calibration/Hengne/FitRelErrors.png");
    c1->Clear();
    relUncertsP->SetMarkerStyle(20);
    relUncertsP->Draw("P");
    //c1->SetLogy(true);
    c1->Print("Calibration/Hengne/FitRelErrorsP.png");
    c1->Clear();
    relUncertsN->SetMarkerStyle(20);
    relUncertsN->Draw("P");
    //c1->SetLogy(true);
    c1->Print("Calibration/Hengne/FitRelErrorsN.png");
    c1->Close();
    
    return 0;
}