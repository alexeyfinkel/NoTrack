
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

int testConsts()
{
    gROOT->SetStyle("Plain");	
	gStyle->SetErrorX(0);
	gStyle->SetOptStat(0);
	
	TCanvas* c1 = new TCanvas("C1", "c1", 1200, 900);
	c1->cd();
    
    TH2D *cumulConstsN = new TH2D("cumulConstsN","Cumulative Consts, Data, NT+;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *cumulConstsP = new TH2D("cumulConstsP","Cumulative Consts, Data, NT-;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *projConstsN = new TH2D("projConstsN","Projected Consts, Data, NT+;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *projConstsP = new TH2D("projConstsP","Projected Consts, Data, NT-;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *corrConstsN = new TH2D("corrConstsN","Corrected Consts, Data, NT+;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *corrConstsP = new TH2D("corrConstsP","Corrected Consts, Data, NT-;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *corrFactorN = new TH2D("corrFactorN","Correction Factors, NT+;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *corrFactorP = new TH2D("corrFactorP","Correction Factors, NT-;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *zBeforeN = new TH2D("zBeforeN","Z-Peaks Before Calibration, NT+;M_{Z}",40,30,70,40,30,70);
    TH2D *zBeforeP = new TH2D("zBeforeP","Z-Peaks Before Calibration, NT-;M_{Z}",40,30,70,40,30,70);
    TH2D *zAfterN = new TH2D("zAfterN","Z-Peaks After Calibration, NT+;M_{Z}",40,30,70,40,30,70);
    TH2D *zAfterP = new TH2D("zAfterP","Z-Peaks After Calibration, NT-;M_{Z}",40,30,70,40,30,70);
    
    TH1D *mZbefore = new TH1D("mZbefore",";M_{Z} (GeV);Evts/Gev",15,80,110);
    TH1D *mZafterProj = new TH1D("mZafterProj",";M_{Z} (GeV);Evts/Gev",15,80,110);
    TH1D *mZafterCumul = new TH1D("mZafterCumul",";M_{Z} (GeV);Evts/Gev",15,80,110);
    
    char title[256],name[128], filename[128];
    std::vector<TH1D*> rawZbyPU, cumulCorrZByPU, projCorrZByPU, corrCorrZByPU;
    for(int bin=0;bin<4;bin++)
    {
        sprintf(name,"rawZ_pu%d",bin);
		sprintf(title,"Uncalibrated M_{Z}, PU %d-%d;M_{ee};Events",5*(bin+1),5*(bin+2));
        rawZbyPU.push_back(new TH1D(name,title,30,60,120));
        sprintf(name,"cumulCorrZ_pu%d",bin);
		sprintf(title,"M_{Z} After Cumulative Correction, PU %d-%d;M_{ee};Events",5*(bin+1),5*(bin+2));
        cumulCorrZByPU.push_back(new TH1D(name,title,30,60,120));
        sprintf(name,"projCorrZ_pu%d",bin);
		sprintf(title,"M_{Z} After Projected Correction, PU %d-%d;M_{ee};Events",5*(bin+1),5*(bin+2));
        projCorrZByPU.push_back(new TH1D(name,title,30,60,120));
        sprintf(name,"corrCorrZ_pu%d",bin);
		sprintf(title,"M_{Z} After MC-Corrected Correction, PU %d-%d;M_{ee};Events",5*(bin+1),5*(bin+2));
        corrCorrZByPU.push_back(new TH1D(name,title,30,60,120));
    }
    
    //Z-hist map
    std::map< std::pair<int,int>, TH1D* > zPeaksBeforeN, zPeaksBeforeP, zPeaksAfterN, zPeaksAfterP;
    for(int i=30;i<71;i++)
	{
		for(int j=30;j<71;j++)
		{
            if( (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5))  < 132.25) || (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5)) > 342.25) ) continue;
            sprintf(name,"zPeaksBefore_ix%diy%dzP",i,j);
            sprintf(title,"M_{Z} Before, NT+ i_{x} = %d, i_{y} = %d;M_{Z}",i,j);
            zPeaksBeforeP[std::make_pair(i,j)]=new TH1D(name,title,40,70,110);
            sprintf(name,"zPeaksBefore_ix%diy%dzN",i,j);
            sprintf(title,"M_{Z} Before, NT- i_{x} = %d, i_{y} = %d;M_{Z}",i,j);
            zPeaksBeforeN[std::make_pair(i,j)]=new TH1D(name,title,40,70,110);
            sprintf(name,"zPeaksAfter_ix%diy%dzP",i,j);
            sprintf(title,"M_{Z} After, NT+ i_{x} = %d, i_{y} = %d;M_{Z}",i,j);
            zPeaksAfterP[std::make_pair(i,j)]=new TH1D(name,title,40,70,110);
            sprintf(name,"zPeaksAfter_ix%diy%dzN",i,j);
            sprintf(title,"M_{Z} After, NT- i_{x} = %d, i_{y} = %d;M_{Z}",i,j);
            zPeaksAfterN[std::make_pair(i,j)]=new TH1D(name,title,40,70,110);
        }
    }
    
    ifstream cumulConstsFile, projConstsFile, corrFactorFile;
    char line[1000];
    int cix, ciy, iz;
    float calConst, sigma;
    
    projConstsFile.open("Calibration/PU/CalConsts_PUregress_Data.txt");
	if(!projConstsFile.is_open())
	{
		std::cout<<"Failed to open Projected Consts. file. Existing"<<std::endl;
		return 1;
	}
    do
    {
        projConstsFile.getline(line, 1000)>>cix>>ciy>>iz>>calConst>>sigma;
        if( (((cix-50.5)*(cix-50.5)+(ciy-50.5)*(ciy-50.5))  < 132.25) || (((cix-50.5)*(cix-50.5)+(ciy-50.5)*(ciy-50.5)) > 342.25) ) continue;
        if(iz==1) projConstsP->SetBinContent(cix-30,ciy-30,calConst);
        else if(iz==-1) projConstsN->SetBinContent(cix-30,ciy-30,calConst);
        else {std::cout<<"Error reading in a cal const. Exiting."<<std::endl;return 1;}
        
    } while(!projConstsFile.eof());
    
    cumulConstsFile.open("Calibration/CalConsts_Data.txt");
	if(!cumulConstsFile.is_open())
	{
		std::cout<<"Failed to open Cumulative Consts. file. Existing"<<std::endl;
		return 1;
	}
    do
    {
        cumulConstsFile.getline(line, 1000)>>cix>>ciy>>iz>>calConst>>sigma;
        if( (((cix-50.5)*(cix-50.5)+(ciy-50.5)*(ciy-50.5))  < 132.25) || (((cix-50.5)*(cix-50.5)+(ciy-50.5)*(ciy-50.5)) > 342.25) ) continue;
        if(iz==1) cumulConstsP->SetBinContent(cix-30,ciy-30,calConst);
        else if(iz==-1) cumulConstsN->SetBinContent(cix-30,ciy-30,calConst);
        else {std::cout<<"Error reading in a cal const. Exiting."<<std::endl;return 1;}
        
    } while(!cumulConstsFile.eof());
    
    corrFactorFile.open("Calibration/DataToMCcorrectionFactors.txt");
	if(!corrFactorFile.is_open())
	{
		std::cout<<"Failed to open Corr. Factors file. Existing"<<std::endl;
		return 1;
	}
    do
    {
        corrFactorFile.getline(line, 1000)>>cix>>ciy>>iz>>calConst;
        if( (((cix-50.5)*(cix-50.5)+(ciy-50.5)*(ciy-50.5))  < 132.25) || (((cix-50.5)*(cix-50.5)+(ciy-50.5)*(ciy-50.5)) > 342.25) ) continue;
        if(iz==1) corrFactorP->SetBinContent(cix-30,ciy-30,calConst);
        else if(iz==-1) corrFactorN->SetBinContent(cix-30,ciy-30,calConst);
        else {std::cout<<"Error creating corrections file. Exiting."<<std::endl;return 1;}
        
    } while(!corrFactorFile.eof());
    
    //make corrected consts. Note: That's Projected Data*Corr. factor to match the MC consts
    for(int i=30;i<71;i++)
	{
		for(int j=30;j<71;j++)
		{
            if( (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5))  < 132.25) || (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5)) > 342.25) ) continue;
            corrConstsN->SetBinContent(i-30,j-30, projConstsN->GetBinContent(i-30,j-30)*corrFactorN->GetBinContent(i-30,j-30) );
            corrConstsP->SetBinContent(i-30,j-30, projConstsP->GetBinContent(i-30,j-30)*corrFactorP->GetBinContent(i-30,j-30) );
        }
    }
    
    //grab a data ntuple
	TFile* f2 = new TFile("/local/cms/user/finkel/NoTrack/Ntuple/DE_ReReco_2012Full_WithClusters_001.root");	
	if(f2 == NULL)
	{
		std::cout<<"Failed to open Data file. Exiting."<<std::endl;
		return 1;
	}
    
    //tracking averages:
    std::vector<int> evtsInBin(4,0), nVertInBin(4,0);
    //Loop over events and apply different calibrations: cumulative, Projected, and Corrected Projected
    //make data ntuple:
    ZEffTree* ze2 = new ZEffTree(*f2,false);
    //Fill the data histograms:
    int ix,iy;
    int nPu, puBin;
    int hix, hiy;
    double correction1, correction2, correction3;//cumulative, projected, and corrected-projected

    TLorentzVector elec1, elec2, theZ;
    float pt, eta, phi, E;

    for(int event=0; event<ze2->Entries(); event++ ) //fill the data hists
    {
        ix=ze2->reco.ix[1];//seed ix
        iy=ze2->reco.iy[1];//seed iy
        nPu=ze2->reco.nverts-1;
        if(nPu<0 || nPu>100) std::cout<<"Warning in event # "<<event<<"!"<<std::endl;
        if(nPu>=25) puBin=3;
        else if(nPu<=5 ) puBin=0;
        else puBin=(int)((nPu-5)/5);

        if(event>0)
        {
            nVertInBin[puBin]+= nPu;
            evtsInBin[puBin]++;
        }

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
            else {correction1=0;correction2=0;correction3=0;}
            for(unsigned int k=0; k<(ze2->ixs->size()); k++ )
            {
                hix = ze2->ixs->at(k);
                hiy = ze2->iys->at(k);
                if( (((hix-50.5)*(hix-50.5)+(hiy-50.5)*(hiy-50.5))  < 132.25) || (((hix-50.5)*(hix-50.5)+(hiy-50.5)*(hiy-50.5)) > 342.25) )
                {
                    correction1+=ze2->hitEnergyFractions->at(k);
                    correction2+=ze2->hitEnergyFractions->at(k);
                    correction3+=ze2->hitEnergyFractions->at(k);
                    continue;
                }
                if(ze2->reco.eta[1]>0)//positive endcap
                {
                    correction1+=(cumulConstsP->GetBinContent(hix-30,hiy-30))*(ze2->hitEnergyFractions->at(k)); 
                    correction2+= (projConstsP->GetBinContent(hix-30,hiy-30))*(ze2->hitEnergyFractions->at(k));
                    correction3+= (corrConstsP->GetBinContent(hix-30,hiy-30))*(ze2->hitEnergyFractions->at(k));
                }

                else //negative endcap
                {
                    correction1+=(cumulConstsN->GetBinContent(hix-30,hiy-30))*(ze2->hitEnergyFractions->at(k)); 
                    correction2+= (projConstsN->GetBinContent(hix-30,hiy-30))*(ze2->hitEnergyFractions->at(k));
                    correction3+= (corrConstsN->GetBinContent(hix-30,hiy-30))*(ze2->hitEnergyFractions->at(k));
                }
            }
            
            rawZbyPU[puBin]->Fill(ze2->reco.mz);
            mZbefore->Fill(ze2->reco.mz);
            //recalculate Z mass using found corrections and fill hists
            pt = ze2->reco.pt[0];    	//tracked electron
            eta = ze2->reco.eta[0];
            phi = ze2->reco.phi[0];
            E = pt*cosh(eta);
            elec1.SetPtEtaPhiE(pt,eta,phi,E);				
            pt = ze2->reco.pt[1]; 		//now untracked
            pt *= correction1;	
            eta = ze2->reco.eta[1];
            phi = ze2->reco.phi[1];    
            E = pt*cosh(eta);	
            elec2.SetPtEtaPhiE(pt,eta,phi,E);
            theZ = elec1+elec2;		//new Z vector
            cumulCorrZByPU[puBin]->Fill(theZ.M());//cumulative correction
            mZafterCumul->Fill(theZ.M());
            //fill Z-before and -after hists
            if( (ze2->reco.mz>80) && (ze2->reco.mz<110) )
            {
                if(ze2->reco.eta[1]>0) 
                {
                    zPeaksBeforeP[std::make_pair(ix,iy)]->Fill(ze2->reco.mz);
                    zPeaksAfterP[std::make_pair(ix,iy)]->Fill(theZ.M());
                }
                else 
                {
                    zPeaksBeforeN[std::make_pair(ix,iy)]->Fill(ze2->reco.mz);
                    zPeaksAfterN[std::make_pair(ix,iy)]->Fill(theZ.M());
                }
            }
            pt = ze2->reco.pt[1]; 
            pt *= correction2;	    
            E = pt*cosh(eta);	
            elec2.SetPtEtaPhiE(pt,eta,phi,E);
            theZ = elec1+elec2;		//new Z vector
            projCorrZByPU[puBin]->Fill(theZ.M());// projected correction
            mZafterProj->Fill(theZ.M());
            
            pt = ze2->reco.pt[1];
            pt *= correction3;	  
            E = pt*cosh(eta);	
            elec2.SetPtEtaPhiE(pt,eta,phi,E);
            theZ = elec1+elec2;		//new Z vector
            corrCorrZByPU[puBin]->Fill(theZ.M());//mc-corrected projected correction

        }
        ze2->GetNextEvent();
    }//end loop over events in ntuple
    
    //pre-reqs for Z-plots:
    std::vector<double> binAveragePU(4,0), zPeaksCumul(4,0), zPeaksProj(4,0), zPeaksCorr(4,0);
    TF1 *voigtFit1 = new TF1("fit1",voigt,75,105,6);
    TF1 *bkgd = new TF1("bkgd","[0]+[1]*(x-80)",80,110 );
    TF1 *puFit = new TF1("fit","pol1",5,25);
    TF1 *gausFit = new TF1("fit","gaus",80,100);
    
    //fit all the Z-hists
    //voigtFit1->SetParameter(1,90);
    //voigtFit1->SetParLimits(1,90,100);
    //voigtFit1->SetParameter(2,5);
    //voigtFit1->SetParLimits(2,1,20);
    //voigtFit1->FixParameter(3,4.9904);
    //voigtFit1->SetParameter(4,1000);
    //voigtFit1->SetParameter(5,-10);
    for(int i=30;i<71;i++)
	{
		for(int j=30;j<71;j++)
		{
            if( (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5))  < 132.25) || (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5)) > 342.25) ) continue;
            zPeaksBeforeP[std::make_pair(i,j)]->Fit(gausFit,"QLMN","",80,100);
            //zPeaksBeforeP[std::make_pair(i,j)]->Draw("P");
            //c1->Update();
            //c1->Clear();
            zBeforeP->SetBinContent( i-30,j-30,gausFit->GetParameter(1) );
            zPeaksBeforeN[std::make_pair(i,j)]->Fit(gausFit,"QLMN","",80,100);
            zBeforeN->SetBinContent( i-30,j-30,gausFit->GetParameter(1) );
            zPeaksAfterP[std::make_pair(i,j)]->Fit(gausFit,"QLMN","",80,100);
            zAfterP->SetBinContent( i-30,j-30,gausFit->GetParameter(1) );
            zPeaksAfterN[std::make_pair(i,j)]->Fit(gausFit,"QLMN","",80,100);
            zAfterN->SetBinContent( i-30,j-30,gausFit->GetParameter(1) );
        }
    }
    
    
    TLegend* l1 = new TLegend(0.55,0.70,0.98,0.94);
	l1->SetFillColor(10);
    l1->SetTextSize(0.032);
    
    c1->SetRightMargin(0.02); 
    c1->SetLeftMargin(0.05);
    c1->SetTopMargin(0.06);
    //make lots of Z-plots
    //Full-smaple plots:
    double maxY = mZbefore->GetBinContent(mZbefore->GetMaximumBin())*1.4;
    TH1D *dummy = new TH1D("dummy","",1,75,115);
    dummy->SetTitle("Z-Peak Before and After Calibration;M_{Z} (GeV)");
    dummy->GetYaxis()->SetRangeUser(0,maxY);
    dummy->Draw(""); 
    mZbefore->SetMarkerStyle(21);
    mZbefore->SetMarkerColor(kRed+2);
    mZbefore->SetLineColor(kRed+2);
    mZbefore->SetLineWidth(4);
    mZbefore->SetMarkerSize(1.1);
    voigtFit1->SetParameter(1,90);
    voigtFit1->SetParLimits(1,90,100);
    voigtFit1->SetParameter(2,5);
    voigtFit1->SetParLimits(2,1,20);
    voigtFit1->FixParameter(3,4.9904);
    voigtFit1->SetParameter(4,1000);
    voigtFit1->SetParameter(5,-10);
    voigtFit1->SetLineWidth(4);
    voigtFit1->SetLineColor(kRed+2);
    mZbefore->Fit(voigtFit1,"QLMN","",80,110);
    sprintf(title,"Uncalibrated: #mu = %.3g, #sigma = %.3g",(float)voigtFit1->GetParameter(1),(float)voigtFit1->GetParameter(2) );
    l1->AddEntry(mZbefore,title,"l");
    voigtFit1->DrawCopy("CSAME")->Draw("SAME");
    mZbefore->DrawCopy("SAME P");
    /*mZafterCumul->SetMarkerStyle(21); //needed a graph with only original and projected consts.
    mZafterCumul->SetMarkerColor(kBlue+1);
    mZafterCumul->SetLineColor(kBlue+1);
    mZafterCumul->SetLineWidth(4);
    mZafterCumul->SetMarkerSize(1.1);
    voigtFit1->SetParameter(1,90);
    voigtFit1->SetParLimits(1,90,100);
    voigtFit1->SetParameter(2,5);
    voigtFit1->SetParLimits(2,1,20);
    voigtFit1->FixParameter(3,4.9904);
    voigtFit1->SetParameter(4,1000);
    voigtFit1->SetParameter(5,-10);
    voigtFit1->SetLineWidth(4);
    voigtFit1->SetLineColor(kBlue+1);
    mZafterCumul->Fit(voigtFit1,"QLMN","",80,110);
    sprintf(title,"Cumulative: #mu = %.3g, #sigma = %.3g",(float)voigtFit1->GetParameter(1),(float)voigtFit1->GetParameter(2) );
    l1->AddEntry(mZafterCumul,title,"l");
    voigtFit1->DrawCopy("CSAME")->Draw("SAME");
    mZafterCumul->DrawCopy("SAME P");*/
    mZafterProj->SetMarkerStyle(21);
    mZafterProj->SetMarkerColor(kBlue+2);
    mZafterProj->SetLineColor(kBlue+2);
    mZafterProj->SetLineWidth(4);
    mZafterProj->SetMarkerSize(1.1);
    voigtFit1->SetParameter(1,90);
    voigtFit1->SetParLimits(1,90,100);
    voigtFit1->SetParameter(2,5);
    voigtFit1->SetParLimits(2,1,20);
    voigtFit1->FixParameter(3,4.9904);
    voigtFit1->SetParameter(4,1000);
    voigtFit1->SetParameter(5,-10);
    voigtFit1->SetLineWidth(4);
    voigtFit1->SetLineColor(kBlue+2);
    mZafterProj->Fit(voigtFit1,"QLMN","",80,110);
    sprintf(title,"Calibrated: #mu = %.3g, #sigma = %.3g",(float)voigtFit1->GetParameter(1),(float)voigtFit1->GetParameter(2) );
    l1->AddEntry(mZafterProj,title,"l");
    voigtFit1->DrawCopy("CSAME")->Draw("SAME");
    mZafterProj->DrawCopy("SAME P");    
    l1->Draw();
    c1->Print("Calibration/ProjCumulCompare_Data/UncalAndProjZPeak.png");
    l1->Clear();
    c1->Clear(); 
    
    c1->SetRightMargin(0.1); 
    //plot Zs
    zBeforeP->GetZaxis()->SetRangeUser(80,100);
	zBeforeP->GetZaxis()->SetLabelSize(0.03);
	zBeforeP->Draw("colz");
	c1->Print("Calibration/ProjCumulCompare_Data/ZpeaksBeforeCumulP.png");
	c1->Clear();
    zBeforeN->GetZaxis()->SetRangeUser(80,100);
	zBeforeN->GetZaxis()->SetLabelSize(0.03);
	zBeforeN->Draw("colz");
	c1->Print("Calibration/ProjCumulCompare_Data/ZpeaksBeforeCumulN.png");
	c1->Clear();
    zAfterP->GetZaxis()->SetRangeUser(80,100);
	zAfterP->GetZaxis()->SetLabelSize(0.03);
	zAfterP->Draw("colz");
	c1->Print("Calibration/ProjCumulCompare_Data/ZpeaksAfterCumulP.png");
	c1->Clear();
    zAfterN->GetZaxis()->SetRangeUser(80,100);
	zAfterN->GetZaxis()->SetLabelSize(0.03);
	zAfterN->Draw("colz");
	c1->Print("Calibration/ProjCumulCompare_Data/ZpeaksAfterCumulN.png");
	c1->Clear();
    
    c1->SetTopMargin(0.02);
    //print cal consts for sanity check
    cumulConstsN->SetTitle(";ix;iy");
    cumulConstsN->GetZaxis()->SetRangeUser(0.75,1.25);
	cumulConstsN->GetZaxis()->SetLabelSize(0.035);
	cumulConstsN->Draw("colz");
	c1->Print("Calibration/ProjCumulCompare_Data/CumulativsConstsN.png");
	c1->Clear();
    cumulConstsP->SetTitle(";ix;iy");
    cumulConstsP->GetZaxis()->SetRangeUser(0.75,1.25);
	cumulConstsP->GetZaxis()->SetLabelSize(0.035);
	cumulConstsP->Draw("colz");
	c1->Print("Calibration/ProjCumulCompare_Data/CumulativsConstsP.png");
	c1->Clear();
        
    c1->SetRightMargin(0.05); 
    c1->SetTopMargin(0.07);
    //PU-binned plots:
    for(int k=0;k<4;k++)
    {
        binAveragePU[k] = (double)nVertInBin[k]/(double)evtsInBin[k];
        
        sprintf(title,"M_{Z} with Cumulative Corrections, PU %d-%d;M_{ee} (GeV);Events/GeV",5*(k+1),5*(k+2) );
        cumulCorrZByPU[k]->SetTitle(title);
        cumulCorrZByPU[k]->GetXaxis()->SetRangeUser(75,110);
        cumulCorrZByPU[k]->SetMarkerStyle(21);
        cumulCorrZByPU[k]->SetLineColor(kBlue);
        cumulCorrZByPU[k]->SetMarkerColor(kBlue+3);
        cumulCorrZByPU[k]->SetMarkerSize(1.3);
        voigtFit1->SetParameter(1,90);
        voigtFit1->SetParLimits(1,90,100);
        voigtFit1->SetParameter(2,5);
        voigtFit1->SetParLimits(2,1,20);
        voigtFit1->FixParameter(3,4.9904);
        voigtFit1->SetParameter(4,1000);
        voigtFit1->SetParameter(5,-10);
        voigtFit1->SetLineWidth(4);
        voigtFit1->SetLineColor(kBlue);
        cumulCorrZByPU[k]->Fit(voigtFit1,"QLM","",80,105);
        sprintf(title,"#mu = %.3g, #sigma = %.3g",(float)voigtFit1->GetParameter(1),(float)voigtFit1->GetParameter(2) );
        l1->AddEntry(voigtFit1,title,"");
        sprintf(filename,"Calibration/PUtest/ZmassCumulative_pu%d.png",k);
        cumulCorrZByPU[k]->Draw("P");
        l1->Draw();
        bkgd->SetParameter(0,voigtFit1->GetParameter(4));
        bkgd->SetParameter(1,voigtFit1->GetParameter(5));
        bkgd->SetLineStyle(2);
        bkgd->SetLineWidth(3);
        bkgd->SetLineColor(kRed+3);
        bkgd->Draw("same");
        c1->Print(filename);
        l1->Clear();
        c1->Clear();        
        zPeaksCumul[k]=voigtFit1->GetParameter(1);
        std::cout<<"bin "<<k<<" avPU "<<binAveragePU[k]<<", total nvert in bin "<<nVertInBin[k]<<", events in bin "<<evtsInBin[k]<<std::endl;
        
        sprintf(title,"M_{Z} with Projected Corrections, PU %d-%d;M_{ee} (GeV);Events/GeV",5*(k+1),5*(k+2) );
        projCorrZByPU[k]->SetTitle(title);
        projCorrZByPU[k]->GetXaxis()->SetRangeUser(75,110);
        projCorrZByPU[k]->SetMarkerStyle(21);
        projCorrZByPU[k]->SetLineColor(kBlue);
        projCorrZByPU[k]->SetMarkerColor(kBlue+3);
        projCorrZByPU[k]->SetMarkerSize(1.3);
        voigtFit1->SetParameter(1,90);
        voigtFit1->SetParLimits(1,90,100);
        voigtFit1->SetParameter(2,5);
        voigtFit1->SetParLimits(2,1,20);
        voigtFit1->FixParameter(3,4.9904);
        voigtFit1->SetParameter(4,1000);
        voigtFit1->SetParameter(5,-10);
        voigtFit1->SetLineWidth(4);
        voigtFit1->SetLineColor(kBlue);
        projCorrZByPU[k]->Fit(voigtFit1,"QLM","",80,105);
        sprintf(title,"#mu = %.3g, #sigma = %.3g",(float)voigtFit1->GetParameter(1),(float)voigtFit1->GetParameter(2) );
        l1->AddEntry(voigtFit1,title,"");
        sprintf(filename,"Calibration/PUtest/ZmassProjected_pu%d.png",k);
        projCorrZByPU[k]->Draw("P");
        l1->Draw();
        bkgd->SetParameter(0,voigtFit1->GetParameter(4));
        bkgd->SetParameter(1,voigtFit1->GetParameter(5));
        bkgd->SetLineStyle(2);
        bkgd->SetLineWidth(3);
        bkgd->SetLineColor(kRed+3);
        bkgd->Draw("same");
        c1->Print(filename);
        l1->Clear();
        c1->Clear();        
        zPeaksProj[k]=voigtFit1->GetParameter(1);
        
        sprintf(title,"M_{Z} with MC-Corrected Projected Corrections, PU %d-%d;M_{ee} (GeV);Events/GeV",5*(k+1),5*(k+2) );
        corrCorrZByPU[k]->SetTitle(title);
        corrCorrZByPU[k]->GetXaxis()->SetRangeUser(75,110);
        corrCorrZByPU[k]->SetMarkerStyle(21);
        corrCorrZByPU[k]->SetLineColor(kBlue);
        corrCorrZByPU[k]->SetMarkerColor(kBlue+3);
        corrCorrZByPU[k]->SetMarkerSize(1.3);
        voigtFit1->SetParameter(1,90);
        voigtFit1->SetParLimits(1,90,100);
        voigtFit1->SetParameter(2,5);
        voigtFit1->SetParLimits(2,1,20);
        voigtFit1->FixParameter(3,4.9904);
        voigtFit1->SetParameter(4,1000);
        voigtFit1->SetParameter(5,-10);
        voigtFit1->SetLineWidth(4);
        voigtFit1->SetLineColor(kBlue);
        corrCorrZByPU[k]->Fit(voigtFit1,"QLM","",80,105);
        sprintf(title,"#mu = %.3g, #sigma = %.3g",(float)voigtFit1->GetParameter(1),(float)voigtFit1->GetParameter(2) );
        l1->AddEntry(voigtFit1,title,"");
        sprintf(filename,"Calibration/PUtest/ZmassCorrected_pu%d.png",k);
        corrCorrZByPU[k]->Draw("P");
        l1->Draw();
        bkgd->SetParameter(0,voigtFit1->GetParameter(4));
        bkgd->SetParameter(1,voigtFit1->GetParameter(5));
        bkgd->SetLineStyle(2);
        bkgd->SetLineWidth(3);
        bkgd->SetLineColor(kRed+3);
        bkgd->Draw("same");
        c1->Print(filename);
        l1->Clear();
        c1->Clear();
        zPeaksCorr[k]=voigtFit1->GetParameter(1);
    }
    //make the vs. PU plots:
    TLegend* l2 = new TLegend(0.65,0.65,0.95,0.93);
	l2->SetFillColor(10);
    
    //for make stupid TGraph work
    TH2D* dummy2 = new TH2D("dummy","Z-peak location, Cumulative Cal. Consts, Vs. Average PU;nVert-1;M_{Z}",1,0,25,1,90,95);    
    dummy2->Draw();
    puFit->SetParameter(0,90);
    puFit->SetParameter(1,0);
    puFit->SetParLimits(0,80,100);
    puFit->SetParLimits(1,0,10);
    puFit->SetLineColor(kBlue);
    TGraph* zCumulByPu = new TGraph(4,&binAveragePU[0],&zPeaksCumul[0]);
    zCumulByPu->SetMarkerStyle(20);
    zCumulByPu->SetMarkerSize(1.3);
    zCumulByPu->SetMarkerColor(kBlue+2);
    zCumulByPu->Fit(puFit,"QM","",0,25);
    zCumulByPu->GetXaxis()->SetRangeUser(0,25);
    zCumulByPu->GetYaxis()->SetRangeUser(88,94);
    zCumulByPu->GetYaxis()->SetTitleOffset(1.1);
    zCumulByPu->Draw("SAME P");
    sprintf(title,"Intercept = %.3f",(float)puFit->GetParameter(0) );
    l2->AddEntry(zCumulByPu,title,"");
    sprintf(title,"Slope = %.3f",(float)puFit->GetParameter(1) );
    l2->AddEntry(zCumulByPu,title,"");
    l2->Draw();
    c1->Print("Calibration/PUtest/ZmassCumulativeVsAveragePU.png");
    c1->Clear(); 
    l2->Clear();
    
    dummy2->SetTitle("Z-peak location, Projected Cal. Consts, Vs. Average PU;nVert-1;M_{Z}");    
    dummy2->Draw();
    puFit->SetParameter(0,90);
    puFit->SetParameter(1,0);
    puFit->SetParLimits(0,80,100);
    puFit->SetParLimits(1,0,10);
    puFit->SetLineColor(kBlue);
    TGraph* zProjByPu = new TGraph(4,&binAveragePU[0],&zPeaksProj[0]);
    zProjByPu->SetMarkerStyle(20);
    zProjByPu->SetMarkerSize(1.3);
    zProjByPu->SetMarkerColor(kBlue+2);
    zProjByPu->Fit(puFit,"QM","",0,25);
    zProjByPu->GetXaxis()->SetRangeUser(0,25);
    zProjByPu->GetYaxis()->SetRangeUser(88,94);
    zProjByPu->GetYaxis()->SetTitleOffset(1.1);
    zProjByPu->Draw("SAME P");
    sprintf(title,"Intercept = %.3f",(float)puFit->GetParameter(0) );
    l2->AddEntry(zProjByPu,title,"");
    sprintf(title,"Slope = %.3f",(float)puFit->GetParameter(1) );
    l2->AddEntry(zProjByPu,title,"");
    l2->Draw();
    c1->Print("Calibration/PUtest/ZmassProjectedVsAveragePU.png");
    c1->Clear();
    l2->Clear();
    
    dummy2->SetTitle("Z-peak location, MC-Corrected Cal. Consts, Vs. Average PU;nVert-1;M_{Z}");    
    dummy2->Draw();
    puFit->SetParameter(0,90);
    puFit->SetParameter(1,0);
    puFit->SetParLimits(0,80,100);
    puFit->SetParLimits(1,0,10);
    puFit->SetLineColor(kBlue);
    TGraph* zCorrByPu = new TGraph(4,&binAveragePU[0],&zPeaksCorr[0]);
    zCorrByPu->SetMarkerStyle(20);
    zCorrByPu->SetMarkerSize(1.3);
    zCorrByPu->SetMarkerColor(kBlue+2);
    zCorrByPu->Fit(puFit,"QM","",0,25);
    zCorrByPu->GetXaxis()->SetRangeUser(0,25);
    zCorrByPu->GetYaxis()->SetRangeUser(88,94);
    zCorrByPu->GetYaxis()->SetTitleOffset(1.1);
    zCorrByPu->Draw("SAME P");
    sprintf(title,"Intercept = %.3f",(float)puFit->GetParameter(0) );
    l2->AddEntry(zCorrByPu,title,"");
    sprintf(title,"Slope = %.3f",(float)puFit->GetParameter(1) );
    l2->AddEntry(zCorrByPu,title,"");
    l2->Draw();
    c1->Print("Calibration/PUtest/ZmassCorrectedVsAveragePU.png");
    c1->Clear();
    l2->Clear();
    
    c1->Close();
    return 0;
}
