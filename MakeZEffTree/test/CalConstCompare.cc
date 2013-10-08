//this should compare the cal. consts from the EVEN and ODD datasets.

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



int calConstCompare()
{
    gROOT->SetStyle("Plain");	
	gStyle->SetErrorX(0);
	gStyle->SetOptStat(0);
	
	TCanvas* c1 = new TCanvas("C1", "c1", 1200, 900);
	c1->cd();
    
    TH2D *constsDiffN = new TH2D("ODDSmeansBeforeN","#deltaC, NT-;i_{x};i_{y}",40,30,70,40,30,70);
	TH2D *constsDiffP = new TH2D("ODDSmeansBeforeP","#deltaC, NT+;i_{x};i_{y}",40,30,70,40,30,70);  
    
    TH2D* ODDScalConstsN = new TH2D("ODDScalConstsN","Calib. Consts for Odd Evts. (found from Even), NT-;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D* ODDScalConstsP = new TH2D("ODDScalConstsP","Calib. Consts for Odd Evts. (found from Even), NT+;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D* EVENScalConstsN = new TH2D("EVENScalConstsN","Calib. Consts for Even Evts. (found from Odd), NT-;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D* EVENScalConstsP = new TH2D("EVENScalConstsP","Calib. Consts for Even Evts. (found from Odd), NT+;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D* FullCalConstsN = new TH2D("EVENScalConstsN","Calib. Consts for All Data, NT-;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D* FullCalConstsP = new TH2D("EVENScalConstsP","Calib. Consts for All Data, NT+;i_{x};i_{y}",40,30,70,40,30,70);
    
    char line[1000];
    int cix, ciy, iz;
    float calConst;
    ifstream oddsFile;
    //read in ODDS cal. consts file (read into "EVENS" hists!):
	oddsFile.open("Calibration/ODDS/CalConsts_Data_ODDS.txt");
	if(!oddsFile.is_open())
	{
		std::cout<<"Failed to open evens cal. constants file. Existing"<<std::endl;
		return 1;
	}
    do
    {
        oddsFile.getline(line, 1000)>>cix>>ciy>>iz>>calConst;
        if(iz==-1) EVENScalConstsN->SetBinContent(cix-30,ciy-30,calConst);
        else if(iz==1)EVENScalConstsP->SetBinContent(cix-30,ciy-30,calConst);
        else 
        {
            std::cout<<"Error readin in a cal. const. Exiting"<<std::endl;
            return 1;
        }
        
    } while(!oddsFile.eof());
    
    //now read in EVENS (use for ODDS!):
    ifstream evensFile;
	evensFile.open("Calibration/EVENS/CalConsts_Data_EVENS.txt");
	if(!evensFile.is_open())
	{
		std::cout<<"Failed to open odds cal. constants file. Existing"<<std::endl;
		return 1;
	}
    do
    {
        evensFile.getline(line, 1000)>>cix>>ciy>>iz>>calConst;
        if(iz==-1) ODDScalConstsN->SetBinContent(cix-30,ciy-30,calConst);
        else if(iz==1)ODDScalConstsP->SetBinContent(cix-30,ciy-30,calConst);
        else 
        {
            std::cout<<"Error reading in a cal. const. Exiting"<<std::endl;
            return 1;
        }
    } while(!evensFile.eof());
    //original cal. consts file:
    ifstream constsFile;
	constsFile.open("Calibration/CalConsts_Data.txt");
	if(!evensFile.is_open())
	{
		std::cout<<"Failed to open original cal. constants file. Existing"<<std::endl;
		return 1;
	}
    do
    {
        constsFile.getline(line, 1000)>>cix>>ciy>>iz>>calConst;
        if(iz==-1) FullCalConstsN->SetBinContent(cix-30,ciy-30,calConst);
        else if(iz==1)FullCalConstsP->SetBinContent(cix-30,ciy-30,calConst);
        else 
        {
            std::cout<<"Error reading in a cal. const. Exiting"<<std::endl;
            return 1;
        }
    } while(!constsFile.eof());
    //now make the finalized cal. consts output file:
    ofstream combinedFile;
	combinedFile.open("Calibration/CombinedCalConsts_Data.txt",std::ios::trunc);
	if(!combinedFile.is_open())
	{
		std::cout<<"Failed to create cal. constants file. Existing"<<std::endl;
		return 1;
	}
	combinedFile<<std::setprecision(7)<<std::fixed;
    
    for(int i=30;i<71;i++)
    {
        for(int j=30;j<71;j++)
        {
            if( (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5))  < 132.25) || (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5)) > 342.25) ) continue;
            constsDiffN->SetBinContent(i-30,j-30, fabs(ODDScalConstsN->GetBinContent(i-30,j-30)-EVENScalConstsN->GetBinContent(i-30,j-30))  );
            constsDiffP->SetBinContent(i-30,j-30, fabs(ODDScalConstsP->GetBinContent(i-30,j-30)-EVENScalConstsP->GetBinContent(i-30,j-30))  );  
            
            combinedFile<<i<<"\t"<<j<<"\t-1\t"<<(float)(FullCalConstsN->GetBinContent(i-30,j-30))<<"\t"<<(float)(constsDiffN->GetBinContent(i-30,j-30))/sqrt(2)<<"\n";
            combinedFile<<i<<"\t"<<j<<"\t1\t"<<(float)(FullCalConstsP->GetBinContent(i-30,j-30))<<"\t"<<(float)(constsDiffP->GetBinContent(i-30,j-30))/sqrt(2)<<"\n";
        }
    }
    
    constsDiffN->GetZaxis()->SetRangeUser(0.00001,0.1);
	constsDiffN->GetZaxis()->SetLabelSize(0.03);
	constsDiffN->Draw("colz");
	c1->Print("Calibration/OddsMinusEvensN.png");
	c1->Clear();
    
    constsDiffP->GetZaxis()->SetRangeUser(0.00001,0.1);
	constsDiffP->GetZaxis()->SetLabelSize(0.03);
	constsDiffP->Draw("colz");
	c1->Print("Calibration/OddsMinusEvensP.png");
	c1->Clear();
    
    c1->Close();
    return 0;
}