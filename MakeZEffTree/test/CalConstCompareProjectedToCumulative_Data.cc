//This is meant to plot the projected cal consts in MC to the "truth" from the miscalibration file.

#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGraph.h"

int compareConstsData()
{
    gROOT->SetStyle("Plain");	
	gStyle->SetErrorX(0);
	gStyle->SetOptStat(0);
	
	TCanvas* c1 = new TCanvas("C1", "c1", 1200, 900);
	c1->cd();
    
    TH2D *projConstsN = new TH2D("projConstsN","Projected Consts, Data, NT+;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *projConstsP = new TH2D("projConstsP","Projected Consts, Data, NT-;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *cumulConstsN = new TH2D("cumulConstsN","Cumulative Consts, MC, NT+;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *cumulConstsP = new TH2D("cumulConstsP","Cumulative Consts, MC, NT-;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *projSlopesN = new TH2D("projSlopesN","Projection Slopes, Data, NT+;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *projSlopesP = new TH2D("projSlopesP","Projection Slopes, Data, NT-;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *ratiosP = new TH2D("ratiosP","C_{p}/C_{tot}, Data, NT+;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *ratiosN = new TH2D("ratiosN","C_{p}/C_{tot}, Data, NT-;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *wtdSlopesN = new TH2D("wtdSlopesN","Weighted Slopes, Data, NT+;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *wtdSlopesP = new TH2D("wtdSlopesP","Weighted Slopes, Data, NT-;i_{x};i_{y}",40,30,70,40,30,70);
    //TH2D *projUncN = new TH2D("projUncN","Projected Consts. Uncert, Data, NT-;i_{x};i_{y}",40,30,70,40,30,70);
    //TH2D *projUncP = new TH2D("projUncN","Projected Consts. Uncert, Data, NT+;i_{x};i_{y}",40,30,70,40,30,70);
    //TH2D *cumulUncN = new TH2D("cumulUncM","Cumulative Consts. Uncert, Data, NT-;i_{x};i_{y}",40,30,70,40,30,70);
    //TH2D *cumulUncP = new TH2D("cumulUncP","Cumulative Consts. Uncert, Data, NT+;i_{x};i_{y}",40,30,70,40,30,70);
    
    //add a correlation plot for MC projected vs "truth"
    TH2D *projCumulCorrelation = new TH2D("projCumulCorrelation","Projected vs. Cumulative;Cumulative Const;Projected Const",80,0.80,1.2,80,0.80,1.25);
    
    TH1D *cumulDist = new TH1D("cumulDist","Cumulative Cal. Consts, Data",40,0.7,1.3);
    TH1D *projDist = new TH1D("projDist","Projected Cal. Consts, Data",40,0.7,1.3);
    TH1D *slopeDist = new TH1D("slopeDist","Projection Slopes, Data",40,-0.01,0.005);
    
    TH1D *cumulPull = new TH1D("cumulPull","Cumuative Consts Pull",40,-40,40);
    //TH1D *projPull = new TH1D("projPull","Projected Consts Pull",40,-1,1);
    
    //arrays to make phi graphs
    std::vector<double> phi[8], wtdSlopeInRingN[8], wtdSlopeInRingP[8], cumulConstInRingN[8], cumulConstInRingP[8], projConstInRingN[8], projConstInRingP[8];
    
    //note: proj and MC are the projected consts
    ifstream projConstsFile, cumulConstsFile, slopesFile;
       
    char line[1000];
    int cix, ciy, iz;
    float calConst, sigma;
    
    projConstsFile.open("Calibration/PU/CalConsts_PUregress_Data.txt");
	if(!projConstsFile.is_open())
	{
		std::cout<<"Failed to open Data constants file. Existing"<<std::endl;
		return 1;
	}
    do
    {
        projConstsFile.getline(line, 1000)>>cix>>ciy>>iz>>calConst>>sigma;
        if( (((cix-50.5)*(cix-50.5)+(ciy-50.5)*(ciy-50.5))  < 132.25) || (((cix-50.5)*(cix-50.5)+(ciy-50.5)*(ciy-50.5)) > 342.25) ) continue;
        projDist->Fill(calConst);
        if(iz==1) 
        {
            projConstsN->SetBinContent(cix-30,ciy-30,calConst);
            projConstsN->SetBinError(cix-30,ciy-30,sigma);
        }
        else if(iz==-1) 
        {
            projConstsP->SetBinContent(cix-30,ciy-30,calConst);
            projConstsP->SetBinError(cix-30,ciy-30,sigma);
        }
        else {std::cout<<"Error reading in a cal const. Exiting."<<std::endl;return 1;}
        
    } while(!projConstsFile.eof());
    
    cumulConstsFile.open("Calibration/CalConsts_Data.txt");
	if(!cumulConstsFile.is_open())
	{
		std::cout<<"Failed to open MC constants file. Existing"<<std::endl;
		return 1;
	}
    do
    {
        cumulConstsFile.getline(line, 1000)>>cix>>ciy>>iz>>calConst>>sigma;
        if( (((cix-50.5)*(cix-50.5)+(ciy-50.5)*(ciy-50.5))  < 132.25) || (((cix-50.5)*(cix-50.5)+(ciy-50.5)*(ciy-50.5)) > 342.25) ) continue;
        cumulDist->Fill(calConst);
        if(iz==1) 
        {
            cumulConstsN->SetBinContent(cix-30,ciy-30,calConst);
            cumulConstsN->SetBinError(cix-30,ciy-30,sigma);
        }
        else if(iz==-1) 
        {
            cumulConstsP->SetBinContent(cix-30,ciy-30,calConst);
            cumulConstsP->SetBinError(cix-30,ciy-30,sigma);
        }
        else {std::cout<<"Error reading in a cal const. Exiting."<<std::endl;return 1;}
        
    } while(!cumulConstsFile.eof());
    
    slopesFile.open("Calibration/PU/CalConstSlopes_PUregress_Data.txt");
	if(!slopesFile.is_open())
	{
		std::cout<<"Failed to open Data constants file. Existing"<<std::endl;
		return 1;
	}
    do
    {
        slopesFile.getline(line, 1000)>>cix>>ciy>>iz>>calConst>>sigma;
        if( (((cix-50.5)*(cix-50.5)+(ciy-50.5)*(ciy-50.5))  < 132.25) || (((cix-50.5)*(cix-50.5)+(ciy-50.5)*(ciy-50.5)) > 342.25) ) continue;
        slopeDist->Fill(calConst);
        if(iz==1) 
        {
            projSlopesN->SetBinContent(cix-30,ciy-30,calConst);
            projSlopesN->SetBinError(cix-30,ciy-30,sigma);
        }
        else if(iz==-1) 
        {
            projSlopesP->SetBinContent(cix-30,ciy-30,calConst);
            projSlopesP->SetBinError(cix-30,ciy-30,sigma);
        }
        else {std::cout<<"Error reading in a slope. Exiting."<<std::endl;return 1;}
        
    } while(!slopesFile.eof());
        
    double constsMeanCumul = projDist->GetMean();
    //double constsMeanProj = cumulDist->GetMean();
    
    int ring;
    for(int i=30;i<71;i++)
	{
		for(int j=30;j<71;j++)
		{
            if( (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5))  < 132.25) || (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5)) > 342.25) ) 
            {
                wtdSlopesN->SetBinContent(i-30,j-30, -999999 );
                wtdSlopesP->SetBinContent(i-30,j-30, -999999 );
                projSlopesN->SetBinContent(i-30,j-30, -999999 );
                projSlopesP->SetBinContent(i-30,j-30, -999999 );
                continue;
            }
            ring = (int)(sqrt( (i-50.5)*(i-50.5)+(j-50.5)*(j-50.5) )-11.5);
            phi[ring].push_back(atan2(j-50.5,i-50.5));
            cumulConstInRingN[ring].push_back( cumulConstsN->GetBinContent(i-30,j-30) );
            cumulConstInRingP[ring].push_back( cumulConstsP->GetBinContent(i-30,j-30) );
            projConstInRingN[ring].push_back( projConstsN->GetBinContent(i-30,j-30) );
            projConstInRingP[ring].push_back( projConstsP->GetBinContent(i-30,j-30) );
            //projPull->Fill( (projConstsN->GetBinContent(i-30,j-30)-constsMeanProj)/projConstsN->GetBinError(i-30,j-30) );
            //projPull->Fill( (projConstsP->GetBinContent(i-30,j-30)-constsMeanProj)/projConstsP->GetBinError(i-30,j-30) );
            projCumulCorrelation->Fill( cumulConstsN->GetBinContent(i-30,j-30), projConstsN->GetBinContent(i-30,j-30) );
            projCumulCorrelation->Fill( cumulConstsP->GetBinContent(i-30,j-30), projConstsP->GetBinContent(i-30,j-30) );
            cumulPull->Fill( (cumulConstsN->GetBinContent(i-30,j-30)-constsMeanCumul)/cumulConstsN->GetBinError(i-30,j-30) );
            cumulPull->Fill( (cumulConstsP->GetBinContent(i-30,j-30)-constsMeanCumul)/cumulConstsP->GetBinError(i-30,j-30) );
            ratiosP->SetBinContent(i-30,j-30, projConstsP->GetBinContent(i-30,j-30)/cumulConstsP->GetBinContent(i-30,j-30) );
            ratiosN->SetBinContent(i-30,j-30, projConstsN->GetBinContent(i-30,j-30)/cumulConstsN->GetBinContent(i-30,j-30) );
            wtdSlopesN->SetBinContent(i-30,j-30, projSlopesN->GetBinContent(i-30,j-30)/projConstsN->GetBinContent(i-30,j-30) );
            wtdSlopesP->SetBinContent(i-30,j-30, projSlopesP->GetBinContent(i-30,j-30)/projConstsP->GetBinContent(i-30,j-30) );
            wtdSlopeInRingN[ring].push_back( wtdSlopesN->GetBinContent(i-30,j-30) );
            wtdSlopeInRingP[ring].push_back( wtdSlopesP->GetBinContent(i-30,j-30) );
        }
    }
    
    TF1 *constsFit = new TF1("fit","gaus",-20,20);
    constsFit->SetLineWidth(4);
    constsFit->SetLineColor(kBlue+2);
    
    char title[1000], filename[1000];
    TH1D* dummy = new TH1D("dummy","",1,-3.5,3.5);
    
    //make TGraphs for eta rings
    //std::vector<TGraph*> cumulByPhiN, cumulByPhiP, projByPhiN, projByPhiP, wtdSlopeByPhiN, wtdSlopeByPhiP;
    for(int i=0;i<7;i++)
    {
        TGraph* cumulByPhiN = new TGraph( phi[i].size(),&phi[i][0],&cumulConstInRingN[i][0]);
        TGraph* cumulByPhiP = new TGraph( phi[i].size(),&phi[i][0],&cumulConstInRingP[i][0]);
        TGraph* projByPhiN = new TGraph( phi[i].size(),&phi[i][0],&projConstInRingN[i][0]);
        TGraph* projByPhiP = new TGraph( phi[i].size(),&phi[i][0],&projConstInRingP[i][0]);
        TGraph* wtdSlopeByPhiN = new TGraph( phi[i].size(),&phi[i][0],&wtdSlopeInRingN[i][0]);
        TGraph* wtdSlopeByPhiP = new TGraph( phi[i].size(),&phi[i][0],&wtdSlopeInRingP[i][0]);
        
        dummy->GetYaxis()->SetRangeUser(0.6,1.4);
        sprintf(title,"Cumulative Constants in Ring %d NT-;#phi;C",i+1);
        sprintf(filename,"Calibration/ProjCumulCompare_Data/Rings/CumulConstsInRing%dN.png",i+1);
        cumulByPhiN->SetMarkerStyle(20);
        cumulByPhiN->SetMarkerColor(kGreen+3);
        dummy->SetTitle(title);
        dummy->Draw("");
        cumulByPhiN->Draw("SAME P");
        c1->Update();
        c1->Print(filename);
        c1->Clear();
        sprintf(title,"Cumulative Constants in Ring %d NT+;#phi;C",i+1);
        sprintf(filename,"Calibration/ProjCumulCompare_Data/Rings/CumulConstsInRing%dP.png",i+1);
        cumulByPhiP->SetMarkerStyle(20);
        cumulByPhiP->SetMarkerColor(kGreen+3);
        dummy->SetTitle(title);
        dummy->Draw("");
        cumulByPhiP->Draw("SAME P");
        c1->Update();
        c1->Print(filename);
        c1->Clear();
        sprintf(title,"Projected Constants in Ring %d NT-;#phi;C",i+1);
        sprintf(filename,"Calibration/ProjCumulCompare_Data/Rings/ProjConstsInRing%dN.png",i+1);
        projByPhiN->SetMarkerStyle(20);
        projByPhiN->SetMarkerColor(kBlue+2);
        dummy->SetTitle(title);
        dummy->Draw("");
        projByPhiN->Draw("SAME P");
        c1->Update();
        c1->Print(filename);
        c1->Clear();
        sprintf(title,"Projected Constants in Ring %d NT+;#phi;C",i+1);
        sprintf(filename,"Calibration/ProjCumulCompare_Data/Rings/ProjConstsInRing%dP.png",i+1);
        projByPhiP->SetMarkerStyle(20);
        projByPhiP->SetMarkerColor(kBlue+2);
        dummy->SetTitle(title);
        dummy->Draw("");
        projByPhiP->Draw("SAME P");
        c1->Update();
        c1->Print(filename);
        c1->Clear();
        dummy->GetYaxis()->SetRangeUser(-0.02,0.02);
        sprintf(title,"Weighted Slopes in Ring %d NT-;#phi;C",i+1);
        sprintf(filename,"Calibration/ProjCumulCompare_Data/Rings/WeightedSlopesInRing%dN.png",i+1);
        wtdSlopeByPhiN->SetMarkerStyle(20);
        wtdSlopeByPhiN->SetMarkerColor(kRed+3);
        dummy->SetTitle(title);
        dummy->Draw("");
        wtdSlopeByPhiN->Draw("SAME P");
        c1->Update();
        c1->Print(filename);
        c1->Clear();
        sprintf(title,"Weighted Slopes in Ring %d NT+;#phi;C",i+1);
        sprintf(filename,"Calibration/ProjCumulCompare_Data/Rings/WeightedSlopesInRing%dP.png",i+1);
        wtdSlopeByPhiP->SetMarkerStyle(20);
        wtdSlopeByPhiP->SetMarkerColor(kRed+3);
        dummy->SetTitle(title);
        dummy->Draw("");
        wtdSlopeByPhiP->Draw("SAME P");
        c1->Update();
        c1->Print(filename);
        c1->Clear();
        
        delete cumulByPhiN;
        delete cumulByPhiP;
        delete projByPhiN;
        delete projByPhiP;
        delete wtdSlopeByPhiN;
        delete wtdSlopeByPhiP;
    }
    
    //make all the plots here! Put into Calibration/ProjCumulCompare_Data/ folder
    c1->SetMargin(0.08,0.02,0.08,0.08);
    TLegend* l1 = new TLegend(0.75,0.7,0.98,0.92);
	l1->SetFillColor(10);
    
    //proj-cumul comparison plot
    projDist->GetYaxis()->SetRangeUser(0,120);
    projDist->SetMarkerStyle(20);
    projDist->SetMarkerSize(1.3);
    projDist->SetMarkerColor(kBlue+1);
    projDist->SetLineColor(kBlue+1);
    projDist->SetLineWidth(4);
    projDist->Draw("P");
    projDist->Fit(constsFit,"QLMN","",0.8,1.1);
    constsFit->SetLineWidth(4);
    constsFit->SetLineColor(kBlue+1);
    constsFit->DrawCopy("SAMEC");
    l1->AddEntry(projDist,"Projected","l");
    cumulDist->SetMarkerColor(kRed);
    cumulDist->SetMarkerSize(1.3);
    cumulDist->SetMarkerStyle(21);
    cumulDist->SetLineColor(kRed);
    cumulDist->SetLineWidth(4);
    cumulDist->Draw("SAME P");
    cumulDist->Fit(constsFit,"QLMN","",0.8,1.1);
    constsFit->SetLineColor(kRed+1);
    constsFit->DrawCopy("SAMEC");
    l1->AddEntry(cumulDist,"Cumulative","l");
    l1->Draw();
    c1->Print("Calibration/ProjCumulCompare_Data/CumulativeVsProjectedDist.png");
    l1->Clear();
    c1->Clear();
    
    //projPull->SetMarkerStyle(20);
    //projPull->Fit(constsFit,"QLM","",-1,1);
    //projPull->Draw("P");
    //sprintf(title,"#sigma = %f",(float)constsFit->GetParameter(2) );
    //l1->AddEntry(projPull,title,"");
    //l1->Draw();
    //c1->Print("Calibration/ProjCumulCompare_Data/ProjPulls.png");
    //c1->Clear();
    //l1->Clear();
    cumulPull->SetMarkerStyle(20);
    cumulPull->Fit(constsFit,"QLM","",-20,20);
    cumulPull->Draw("P");
    sprintf(title,"#sigma = %f",(float)constsFit->GetParameter(2) );
    l1->AddEntry(cumulPull,title,"");
    l1->Draw();
    c1->Print("Calibration/ProjCumulCompare_Data/CumulPulls.png");
    c1->Clear();
    l1->Clear();
    projDist->SetMarkerStyle(20);
    projDist->Fit(constsFit,"QLM","",0.8,1.15);
    projDist->Draw("P");
    sprintf(title,"#mu = %f",(float)constsFit->GetParameter(1) );
    l1->AddEntry(projDist,title,"");
    sprintf(title,"#sigma = %f",(float)constsFit->GetParameter(2) );
    l1->AddEntry(projDist,title,"");
    l1->Draw();
    c1->Print("Calibration/ProjCumulCompare_Data/ProjConstsDist.png");
    c1->Clear();
    l1->Clear();
    cumulDist->SetMarkerStyle(20);
    cumulDist->Fit(constsFit,"QLM","",0.8,1.15);
    cumulDist->Draw("P");
    sprintf(title,"#mu = %f",(float)constsFit->GetParameter(1) );
    l1->AddEntry(projDist,title,"");
    sprintf(title,"#sigma = %f",(float)constsFit->GetParameter(2) );
    l1->AddEntry(projDist,title,"");
    l1->Draw();
    c1->Print("Calibration/ProjCumulCompare_Data/CumulConstsDist.png");
    c1->Clear();
    l1->Clear();
    
    projConstsN->GetZaxis()->SetRangeUser(0.75,1.25);
	projConstsN->GetZaxis()->SetLabelSize(0.03);
	projConstsN->Draw("colz");
	c1->Print("Calibration/ProjCumulCompare_Data/ProjConstsN.png");
	c1->Clear();
    projConstsP->GetZaxis()->SetRangeUser(0.75,1.25);
	projConstsP->GetZaxis()->SetLabelSize(0.03);
	projConstsP->Draw("colz");
	c1->Print("Calibration/ProjCumulCompare_Data/ProjConstsP.png");
	c1->Clear();
    cumulConstsN->GetZaxis()->SetRangeUser(0.75,1.25);
	cumulConstsN->GetZaxis()->SetLabelSize(0.03);
	cumulConstsN->Draw("colz");
	c1->Print("Calibration/ProjCumulCompare_Data/CumulConstsN.png");
	c1->Clear();
    cumulConstsP->GetZaxis()->SetRangeUser(0.75,1.25);
	cumulConstsP->GetZaxis()->SetLabelSize(0.03);
	cumulConstsP->Draw("colz");
	c1->Print("Calibration/ProjCumulCompare_Data/CumulConstsP.png");
	c1->Clear();
    projSlopesN->GetZaxis()->SetRangeUser(-0.01,0.005);
	projSlopesN->GetZaxis()->SetLabelSize(0.03);
	projSlopesN->Draw("colz");
	c1->Print("Calibration/ProjCumulCompare_Data/SlopesN.png");
	c1->Clear();
    projSlopesP->GetZaxis()->SetRangeUser(-0.01,0.005);
	projSlopesP->GetZaxis()->SetLabelSize(0.03);
	projSlopesP->Draw("colz");
	c1->Print("Calibration/ProjCumulCompare_Data/SlopesP.png");
	c1->Clear();
    wtdSlopesN->GetZaxis()->SetRangeUser(-0.01,0.005);
	wtdSlopesN->GetZaxis()->SetLabelSize(0.03);
	wtdSlopesN->Draw("colz");
	c1->Print("Calibration/ProjCumulCompare_Data/WeightedSlopesN.png");
	c1->Clear();
    wtdSlopesP->GetZaxis()->SetRangeUser(-0.01,0.005);
	wtdSlopesP->GetZaxis()->SetLabelSize(0.03);
	wtdSlopesP->Draw("colz");
	c1->Print("Calibration/ProjCumulCompare_Data/WeightedSlopesP.png");
	c1->Clear();
    ratiosN->GetZaxis()->SetRangeUser(0.95,1.1);
	ratiosN->GetZaxis()->SetLabelSize(0.03);
	ratiosN->Draw("colz");
	c1->Print("Calibration/ProjCumulCompare_Data/RatiosN.png");
	c1->Clear();
    ratiosP->GetZaxis()->SetRangeUser(0.95,1.1);
	ratiosP->GetZaxis()->SetLabelSize(0.03);
	ratiosP->Draw("colz");
	c1->Print("Calibration/ProjCumulCompare_Data/RatiosP.png");
	c1->Clear();
     
	projCumulCorrelation->GetZaxis()->SetLabelSize(0.03);
	projCumulCorrelation->Draw("colz");
    sprintf(title,"C = %f",(float)projCumulCorrelation->GetCorrelationFactor() );
    l1->AddEntry(projCumulCorrelation,title,"");
    l1->Draw();
	c1->Print("Calibration/ProjCumulCompare_Data/ProjectedToTCumulativeCorrelation.png");
	c1->Clear();
        
    c1->Close();
    return 0;
}