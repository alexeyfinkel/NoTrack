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

int compareConsts()
{
    gROOT->SetStyle("Plain");	
	gStyle->SetErrorX(0);
	gStyle->SetOptStat(0);
	
	TCanvas* c1 = new TCanvas("C1", "c1", 1200, 900);
	c1->cd();
    
    TH2D *dataConstsN = new TH2D("dataConstsN","Projected Consts, Data, NT+;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *dataConstsP = new TH2D("dataConstsP","Projected Consts, Data, NT-;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *dataUncN = new TH2D("dataUncN","Projected Const. Uncerts, Data, NT+;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *dataUncP = new TH2D("dataUncP","Projected Const. Uncerts, Data, NT-;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *mcConstsN = new TH2D("mcConstsN","Projected Consts, MC, NT+;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *mcConstsP = new TH2D("mcConstsP","Projected Consts, MC, NT-;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *cumulConstsN = new TH2D("cumulConstsN","Cumulative Consts, MC, NT+;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *cumulConstsP = new TH2D("cumulConstsP","Cumulative Consts, MC, NT-;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *mcUncN = new TH2D("mcUncN","Projected Const. Uncerts, MC, NT+;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *mcUncP = new TH2D("mcUncP","Projected Const. Uncerts, MC, NT-;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *mcTruthN = new TH2D("mcTruthN","MC Smearing, NT+;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *mcTruthP = new TH2D("mcTruthP","MC Smearing, NT-;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *corrFactorN = new TH2D("corrFactorN","Correction Factors, NT+;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *corrFactorP = new TH2D("corrFactorP","Correction Factors, NT-;i_{x};i_{y}",40,30,70,40,30,70);
    
    //add a correlation plot for MC projected vs "truth"
    TH2D *projTruthCorrelation = new TH2D("projTruthCorrelation","MC Smearing vs. Projected;Projected Cnst;MC Smearing",60,0.90,1.06,60,0.9,1.1);
    TH2D *projCumulCorrelation = new TH2D("projCumulCorrelation","Projected vs. Cumulative Consts;Projected Cnst;Cumulative Consts",60,0.9,1.2,60,0.9,1.2);
    TH2D *cumulTruthCorrelation = new TH2D("cumulTruthCorrelation","MC Smearing vs. Cumulative Consts;Cumul Cnst;MC Truth",60,0.9,1.04,60,0.9,1.1);
    TH2D *MCtoDataCorrelation = new TH2D("MCtoDataCorrelation","Data vs. MC;MC Const;Data Const",80,0.75,1.25,60,0.7,1.3);
    
    TH1D *dataDist = new TH1D("dataDist","Data Projected Consts",40,0.7,1.3);
    TH1D *mcDist = new TH1D("mcDist","MC Projected Consts",40,0.7,1.3);
    TH1D *truthDist = new TH1D("truthDist","MC Smearing Consts.",40,0.7,1.3);
    TH1D *correctionDist = new TH1D("correctionDist","1D Correction Factor Distribution",40,0.6,1.4);
    TH1D *cumulDist = new TH1D("truthDist","Cumulative Consts.",40,0.7,1.3);
    
    TH1D *dataPull = new TH1D("dataPull","Data Consts Pull",20,-10,10);
    TH1D *mcPull = new TH1D("mcPull","MC Consts Pull",20,-10,10);
    
    //note: data and MC are the projected consts
    ifstream dataConstsFile, mcConstsFile, mcTruthFile, cumulConstsFile;
    ofstream correctionsFile;
    correctionsFile<<std::setprecision(7)<<std::fixed;
    
    char line[1000];
    int cix, ciy, iz;
    float calConst, sigma;
    
    correctionsFile.open("Calibration/DataToMCcorrectionFactors.txt");
	if(!correctionsFile.is_open())
	{
		std::cout<<"Failed to create Corr. Factors file. Existing"<<std::endl;
		return 1;
	}
    
    dataConstsFile.open("Calibration/PU/CalConsts_PUregress_Data.txt");
	if(!dataConstsFile.is_open())
	{
		std::cout<<"Failed to open Data constants file. Existing"<<std::endl;
		return 1;
	}
    do
    {
        dataConstsFile.getline(line, 1000)>>cix>>ciy>>iz>>calConst>>sigma;
        if( (((cix-50.5)*(cix-50.5)+(ciy-50.5)*(ciy-50.5))  < 132.25) || (((cix-50.5)*(cix-50.5)+(ciy-50.5)*(ciy-50.5)) > 342.25) ) continue;
        dataDist->Fill(calConst);
        if(iz==1) 
        {
            dataConstsN->SetBinContent(cix-30,ciy-30,calConst);
            dataUncN->SetBinContent(cix-30,ciy-30,sigma);
        }
        else if(iz==-1) 
        {
            dataConstsP->SetBinContent(cix-30,ciy-30,calConst);
            dataUncP->SetBinContent(cix-30,ciy-30,sigma);
        }
        else {std::cout<<"Error reading in a cal const. Exiting."<<std::endl;return 1;}
        
    } while(!dataConstsFile.eof());
    
    mcConstsFile.open("Calibration/MC/newMC/PU/CalConsts_PUregress_MC.txt");
	if(!mcConstsFile.is_open())
	{
		std::cout<<"Failed to open MC constants file. Existing"<<std::endl;
		return 1;
	}
    do
    {
        mcConstsFile.getline(line, 1000)>>cix>>ciy>>iz>>calConst>>sigma;
        if( (((cix-50.5)*(cix-50.5)+(ciy-50.5)*(ciy-50.5))  < 132.25) || (((cix-50.5)*(cix-50.5)+(ciy-50.5)*(ciy-50.5)) > 342.25) ) continue;
        mcDist->Fill(calConst);
        if(iz==1) 
        {
            mcConstsN->SetBinContent(cix-30,ciy-30,calConst);
            mcUncN->SetBinContent(cix-30,ciy-30,sigma);
        }
        else if(iz==-1) 
        {
            mcConstsP->SetBinContent(cix-30,ciy-30,calConst);
            mcUncP->SetBinContent(cix-30,ciy-30,sigma);
        }
        else {std::cout<<"Error reading in a cal const. Exiting."<<std::endl;return 1;}
        
    } while(!mcConstsFile.eof());
    
    mcTruthFile.open("Calibration/MC/MC_SmearingConsts.txt");
	if(!mcTruthFile.is_open())
	{
		std::cout<<"Failed to open MC constants file. Existing"<<std::endl;
		return 1;
	}
    do
    {
        mcTruthFile.getline(line, 1000)>>cix>>ciy>>iz>>calConst;
        if( (((cix-50.5)*(cix-50.5)+(ciy-50.5)*(ciy-50.5))  < 132.25) || (((cix-50.5)*(cix-50.5)+(ciy-50.5)*(ciy-50.5)) > 342.25) ) continue;
        truthDist->Fill(calConst);
        if(iz==1) mcTruthN->SetBinContent(cix-30,ciy-30,calConst);
        else if(iz==-1) mcTruthP->SetBinContent(cix-30,ciy-30,calConst);
        else {std::cout<<"Error reading truth file. Exiting."<<std::endl;return 1;}
        
    } while(!mcTruthFile.eof());
    cumulConstsFile.open("Calibration/MC/newMC/PU/CumulativeCalConsts_MC.txt");
	if(!cumulConstsFile.is_open())
	{
		std::cout<<"Failed to open MC constants file. Existing"<<std::endl;
		return 1;
	}
    do
    {
        cumulConstsFile.getline(line, 1000)>>cix>>ciy>>iz>>calConst;
        if( (((cix-50.5)*(cix-50.5)+(ciy-50.5)*(ciy-50.5))  < 132.25) || (((cix-50.5)*(cix-50.5)+(ciy-50.5)*(ciy-50.5)) > 342.25) ) continue;
        cumulDist->Fill(calConst);
        if(iz==1) cumulConstsN->SetBinContent(cix-30,ciy-30,calConst);
        else if(iz==-1) cumulConstsN->SetBinContent(cix-30,ciy-30,calConst);
        else {std::cout<<"Error reading cumuls file. Exiting."<<std::endl;return 1;}
        
    } while(!cumulConstsFile.eof());
    
    double constsMeanData = dataDist->GetMean();
    double constsMeanMC = mcDist->GetMean();
    
    for(int i=30;i<71;i++)
	{
		for(int j=30;j<71;j++)
		{
            if( (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5))  < 132.25) || (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5)) > 342.25) ) continue;
            corrFactorN->SetBinContent(i-30,j-30, mcConstsN->GetBinContent(i-30,j-30)/dataConstsN->GetBinContent(i-30,j-30) );
            corrFactorP->SetBinContent(i-30,j-30, mcConstsP->GetBinContent(i-30,j-30)/dataConstsP->GetBinContent(i-30,j-30) );
            correctionDist->Fill(mcConstsN->GetBinContent(i-30,j-30)/dataConstsN->GetBinContent(i-30,j-30));
            correctionDist->Fill(mcConstsP->GetBinContent(i-30,j-30)/dataConstsP->GetBinContent(i-30,j-30));
            correctionsFile<<i<<"\t"<<j<<"\t-1\t"<<(float)(corrFactorN->GetBinContent(i-30,j-30) )<<"\n";
            correctionsFile<<i<<"\t"<<j<<"\t1\t"<<(float)(corrFactorP->GetBinContent(i-30,j-30) )<<"\n";
            dataPull->Fill( (dataConstsN->GetBinContent(i-30,j-30)-constsMeanData)/dataUncN->GetBinContent(i-30,j-30) );
            dataPull->Fill( (dataConstsP->GetBinContent(i-30,j-30)-constsMeanData)/dataUncP->GetBinContent(i-30,j-30) );
            mcPull->Fill( (mcConstsN->GetBinContent(i-30,j-30)-constsMeanMC)/mcUncN->GetBinContent(i-30,j-30) );
            mcPull->Fill( (mcConstsP->GetBinContent(i-30,j-30)-constsMeanMC)/mcUncP->GetBinContent(i-30,j-30) );
            projTruthCorrelation->Fill( mcConstsN->GetBinContent(i-30,j-30), mcTruthN->GetBinContent(i-30,j-30) );
            projTruthCorrelation->Fill( mcConstsP->GetBinContent(i-30,j-30), mcTruthP->GetBinContent(i-30,j-30) );
            projCumulCorrelation->Fill( mcConstsN->GetBinContent(i-30,j-30), cumulConstsN->GetBinContent(i-30,j-30) );
            projCumulCorrelation->Fill( mcConstsP->GetBinContent(i-30,j-30), cumulConstsP->GetBinContent(i-30,j-30) );
            cumulTruthCorrelation->Fill( cumulConstsN->GetBinContent(i-30,j-30), mcTruthN->GetBinContent(i-30,j-30) );
            cumulTruthCorrelation->Fill( cumulConstsP->GetBinContent(i-30,j-30), mcTruthP->GetBinContent(i-30,j-30) );
            MCtoDataCorrelation->Fill( mcConstsN->GetBinContent(i-30,j-30), dataConstsN->GetBinContent(i-30,j-30) );
            MCtoDataCorrelation->Fill( mcConstsP->GetBinContent(i-30,j-30), dataConstsN->GetBinContent(i-30,j-30) );
        }
    }
    
    //make a pull distribution for Data and for MC
    TF1 *constsFit = new TF1("fit","gaus",-10,10);
    constsFit->SetLineWidth(4);
    constsFit->SetLineColor(kBlue+2);
  
    //make all the plots here! Put into Calibration/DataMCcompare/ folder
    TLegend* l1 = new TLegend(0.7,0.7,0.90,0.90);
	l1->SetFillColor(10);
    l1->SetTextSize(0.03);
    char title[1000];
    
    dataPull->SetMarkerStyle(20);
    dataPull->Fit(constsFit,"QLM","",-8,6);
    dataPull->Draw("P");
    sprintf(title,"#sigma = %f",(float)constsFit->GetParameter(2) );
    l1->AddEntry(dataPull,title,"");
    l1->Draw();
    c1->Print("Calibration/DataMCcompare/DataPulls.png");
    c1->Clear();
    l1->Clear();
    mcPull->SetMarkerStyle(20);
    mcPull->Fit(constsFit,"QLM","",-6,4);
    mcPull->Draw("P");
    sprintf(title,"#sigma = %f",(float)constsFit->GetParameter(2) );
    l1->AddEntry(mcPull,title,"");
    l1->Draw();
    c1->Print("Calibration/DataMCcompare/mcPulls.png");
    c1->Clear();
    l1->Clear();
    dataDist->SetMarkerStyle(20);
    dataDist->Fit(constsFit,"QLM","",0.8,1.2);
    dataDist->Draw("P");
    sprintf(title,"#sigma = %f",(float)constsFit->GetParameter(2) );
    l1->AddEntry(dataDist,title,"");
    l1->Draw();
    c1->Print("Calibration/DataMCcompare/DataProjDist.png");
    c1->Clear();
    l1->Clear();
    mcDist->SetMarkerStyle(20);
    mcDist->Fit(constsFit,"QLMN","",0.9,1.05);
    mcDist->Draw("P");
    constsFit->DrawCopy("CSAME");
    sprintf(title,"#sigma = %f",(float)constsFit->GetParameter(2) );
    l1->AddEntry(mcDist,title,"");
    l1->Draw();
    c1->Print("Calibration/DataMCcompare/MCprojDist.png");
    c1->Clear();
    l1->Clear();
    correctionDist->SetMarkerStyle(20);
    correctionDist->Fit(constsFit,"QLM","",0.7,1.3);
    correctionDist->Draw("P");
    sprintf(title,"#sigma = %f",(float)constsFit->GetParameter(2) );
    l1->AddEntry(correctionDist,title,"");
    l1->Draw();
    c1->Print("Calibration/DataMCcompare/CorrectionsDist.png");
    c1->Clear();
    l1->Clear();
    truthDist->SetMarkerStyle(20);
    truthDist->Fit(constsFit,"QLM","",0.9,1.1);
    truthDist->Draw("P");
    sprintf(title,"#sigma = %.3f",(float)constsFit->GetParameter(2) );
    l1->AddEntry(truthDist,title,"");
    l1->Draw();
    c1->Print("Calibration/DataMCcompare/mcTruthDist.png");
    c1->Clear();
    l1->Clear();
    cumulDist->SetMarkerStyle(20);
    cumulDist->Fit(constsFit,"QLMN","",0.91,1.01);
    cumulDist->Draw("P");
    constsFit->DrawCopy("CSAME");
    sprintf(title,"#sigma = %.3f",(float)constsFit->GetParameter(2) );
    l1->AddEntry(cumulDist,title,"");
    l1->Draw();
    c1->Print("Calibration/DataMCcompare/MCcumulDist.png");
    c1->Clear();
    l1->Clear();
    
    truthDist->SetMarkerStyle(20);
    truthDist->SetMarkerColor(kBlue+2);
    truthDist->SetLineColor(kBlue+2);
    truthDist->SetLineWidth(4);
    truthDist->SetTitle("MC Truth vs. Cuml. vs. Proj. Distributions");
    truthDist->Fit(constsFit,"QLMN","",0.9,1.1);
    truthDist->Draw("P");
    constsFit->DrawCopy("SAMEC");
    sprintf(title,"Truth #sigma = %.3f",(float)constsFit->GetParameter(2) );
    l1->AddEntry(truthDist,title,"l");
    cumulDist->SetMarkerStyle(20);
    cumulDist->SetMarkerColor(kRed+2);
    cumulDist->SetLineColor(kRed+2);
    constsFit->SetLineColor(kRed+2);
    cumulDist->SetLineWidth(4);
    cumulDist->Fit(constsFit,"QLMN","",0.91,1.01);
    cumulDist->Draw("Same P");
    constsFit->DrawCopy("SAMEC");
    sprintf(title,"Cuml. #sigma = %.3f",(float)constsFit->GetParameter(2) );
    l1->AddEntry(cumulDist,title,"l");
    mcDist->SetMarkerStyle(20);
    mcDist->SetMarkerColor(kBlack);
    mcDist->SetLineColor(kBlack);
    constsFit->SetLineColor(kBlack);
    mcDist->SetLineWidth(4);
    mcDist->Fit(constsFit,"QLMN","",0.91,1.01);
    mcDist->Draw("Same P");
    constsFit->DrawCopy("SAMEC");
    sprintf(title,"Proj. #sigma = %.3f",(float)constsFit->GetParameter(2) );
    l1->AddEntry(mcDist,title,"l");
    l1->Draw();
    c1->Print("Calibration/DataMCcompare/TruthCumulProjectedDist.png");
    c1->Clear();
    l1->Clear();
    
    dataConstsN->GetZaxis()->SetRangeUser(0.75,1.25);
	dataConstsN->GetZaxis()->SetLabelSize(0.03);
	dataConstsN->Draw("colz");
	c1->Print("Calibration/DataMCcompare/DataConstsN.png");
	c1->Clear();
    dataConstsP->GetZaxis()->SetRangeUser(0.75,1.25);
	dataConstsP->GetZaxis()->SetLabelSize(0.03);
	dataConstsP->Draw("colz");
	c1->Print("Calibration/DataMCcompare/DataConstsP.png");
	c1->Clear();
    mcConstsN->GetZaxis()->SetRangeUser(0.75,1.25);
	mcConstsN->GetZaxis()->SetLabelSize(0.03);
	mcConstsN->Draw("colz");
	c1->Print("Calibration/DataMCcompare/mcConstsN.png");
	c1->Clear();
    mcConstsP->GetZaxis()->SetRangeUser(0.75,1.25);
	mcConstsP->GetZaxis()->SetLabelSize(0.03);
	mcConstsP->Draw("colz");
	c1->Print("Calibration/DataMCcompare/mcConstsP.png");
	c1->Clear();
    mcTruthN->GetZaxis()->SetRangeUser(0.75,1.25);
	mcTruthN->GetZaxis()->SetLabelSize(0.03);
	mcTruthN->Draw("colz");
	c1->Print("Calibration/DataMCcompare/mcTruthN.png");
	c1->Clear();
    mcTruthP->GetZaxis()->SetRangeUser(0.75,1.25);
	mcTruthP->GetZaxis()->SetLabelSize(0.03);
	mcTruthP->Draw("colz");
	c1->Print("Calibration/DataMCcompare/mcTruthP.png");
	c1->Clear();
    corrFactorN->GetZaxis()->SetRangeUser(0.75,1.25);
	corrFactorN->GetZaxis()->SetLabelSize(0.03);
	corrFactorN->Draw("colz");
	c1->Print("Calibration/DataMCcompare/CorrectionFactorsN.png");
	c1->Clear();
    corrFactorP->GetZaxis()->SetRangeUser(0.75,1.25);
	corrFactorP->GetZaxis()->SetLabelSize(0.03);
	corrFactorP->Draw("colz");
	c1->Print("Calibration/DataMCcompare/CorrectionFactorsP.png");
	c1->Clear();
    //projTruthCorrelation->GetZaxis()->SetRangeUser(0.9,1.1);
	projTruthCorrelation->GetZaxis()->SetLabelSize(0.03);
	projTruthCorrelation->Draw("colz");
    sprintf(title,"C = %.3f",(float)projTruthCorrelation->GetCorrelationFactor() );
    l1->AddEntry(projTruthCorrelation,title,"");
    l1->Draw();
	c1->Print("Calibration/DataMCcompare/ProjectedToTruthCorrelation.png");
	c1->Clear();
    l1->Clear();
    projCumulCorrelation->GetZaxis()->SetLabelSize(0.03);
	projCumulCorrelation->Draw("colz");
    sprintf(title,"C = %.3f",(float)projCumulCorrelation->GetCorrelationFactor() );
    l1->AddEntry(projCumulCorrelation,title,"");
    l1->Draw();
	c1->Print("Calibration/DataMCcompare/ProjectedToCumulativeCorrelation.png");
	c1->Clear();
    l1->Clear();
    cumulTruthCorrelation->GetZaxis()->SetLabelSize(0.03);
	cumulTruthCorrelation->Draw("colz");
    sprintf(title,"C = %.3f",(float)cumulTruthCorrelation->GetCorrelationFactor() );
    l1->AddEntry(cumulTruthCorrelation,title,"");
    l1->Draw();
	c1->Print("Calibration/DataMCcompare/CumulativeToTruthCorrelation.png");
	c1->Clear();
    l1->Clear();
    //corrFactorP->GetZaxis()->SetRangeUser(0.9,1.1);
	MCtoDataCorrelation->GetZaxis()->SetLabelSize(0.03);
	MCtoDataCorrelation->Draw("colz");
    sprintf(title,"C = %.3f",(float)MCtoDataCorrelation->GetCorrelationFactor() );
    l1->AddEntry(MCtoDataCorrelation,title,"");
    l1->Draw();
	c1->Print("Calibration/DataMCcompare/MCtoDataCorrelation.png");
	c1->Clear();
    l1->Clear();
    
    c1->Close();
    return 0;
}