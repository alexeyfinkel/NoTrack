//this script is for comparing MC and Data Calibration plots
//
//NOTE: To save run time and space, it DOES NOT SAVE any of the individual crystal-level hists!
//Also makes no hitmap-type hists, since can'tcompare on those!
//


#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
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

#include "../src/ZEffTree.h"

int crystalCalibCombined()
{
	gROOT->SetStyle("Plain");	
	gStyle->SetErrorX(0);
	gStyle->SetOptStat(0);
	
	TCanvas* c1 = new TCanvas("C1", "c1", 1200, 900);
	c1->cd();

	TLegend* l1 = new TLegend(0.65,0.7,0.9,0.9);
	l1->SetFillColor(10);
	
	//let's be fancy and use a map!
	std::map< pair<int,int>, TH1D* > histMapDataP, histMapDataN, ratioMapDataP, ratioMapDataN, histMapMCP, histMapMCN, ratioMapMCP, ratioMapMCN;
	char nameMDataP[128], nameMDataN[128], nameRDataN[128], nameRDataP[128], nameMCP[128], nameMCN[128], nameRmcN[128], nameRmcP[128];	
	
	for(int i=30;i<71;i++)
	{
		for(int j=30;j<71;j++)
		{
			//data
			if( (((i-50)*(i-50)+(j-50)*(j-50))  < 100) || (((i-50)*(i-50)+(j-50)*(j-50)) > 324) ) continue;
			sprintf(nameMDataP,"mZDatix%diy%dzP",i,j);
			sprintf(nameMDataN,"mZDatix%diy%dzN",i,j);
			sprintf(nameRDataP,"RatioDatIx%diy%dzP",i,j);
			sprintf(nameRDataN,"RatioDatIx%diy%dzN",i,j);
			histMapDataP[std::make_pair(i,j)]= new TH1D(nameMDataP,"",30,60,120);
			histMapDataN[std::make_pair(i,j)]= new TH1D(nameMDataN,"",30,60,120);
			ratioMapDataP[std::make_pair(i,j)]= new TH1D(nameRDataP,"",20,0.4,1.6);
			ratioMapDataN[std::make_pair(i,j)]= new TH1D(nameRDataN,"",20,0.4,1.6);
			//MC
			sprintf(nameMCP,"mZix%diy%dzP",i,j);
			sprintf(nameMCN,"mZix%diy%dzN",i,j);
			sprintf(nameRmcP,"RatioIx%diy%dzP",i,j);
			sprintf(nameRmcN,"RatioIx%diy%dzN",i,j);
			histMapMCP[std::make_pair(i,j)]= new TH1D(nameMCP,"",30,60,120);
			histMapMCN[std::make_pair(i,j)]= new TH1D(nameMCN,"",30,60,120);
			ratioMapMCP[std::make_pair(i,j)]= new TH1D(nameRmcP,"",20,0.4,1.6);
			ratioMapMCN[std::make_pair(i,j)]= new TH1D(nameRmcN,"",20,0.4,1.6);
			
		}
	}
		
	//getting real ambitious now: gonna do eta rings!
	std::vector<TH1D*> DataMCpeakByEtaP, DataMCpeakByEtaN, DataMCmeanByEtaP, DataMCmeanByEtaN, MCpeakByEtaP, MCpeakByEtaN, MCmeanByEtaP, MCmeanByEtaN;
	char name[128], title[256];// nameMC[128], titleMC[256];
	for(int i=0;i<7;i++)
	{
		//data
		sprintf(name,"peaksDatPEtaRing%d",i+1);
		sprintf(title,"Gauss-fit Peaks, Data, NT+, Eta Ring #%d;GF Peak;Events",i+1);
		DataMCpeakByEtaP.push_back( new TH1D(name,title,30,60,120) );
		sprintf(name,"peaksDatNEtaRing%d",i+1);
		sprintf(title,"Gauss-fit Peaks, Data NT-, Eta Ring #%d;GF Peak;Events",i+1);
		DataMCpeakByEtaN.push_back( new TH1D(name,title,30,60,120) );
		sprintf(name,"meansDatPEtaRing%d",i+1);
		sprintf(title,"Mean Ratios, Data, NT+, Eta Ring #%d;Mean Ratio;Events",i+1);
		DataMCmeanByEtaP.push_back( new TH1D(name,title,40,0,2) );
		sprintf(name,"meansDatNEtaRing%d",i+1);
		sprintf(title,"Mean Ratios, Data, NT-, Eta Ring #%d;Mean Ratio;Events",i+1);
		DataMCmeanByEtaN.push_back( new TH1D(name,title,40,0,2) );
		//MC	
		sprintf(name,"peaksEtaPRing%d",i+1);
		sprintf(title,"Gauss-fit Peaks, NT+, Eta Ring #%d;GF Peak;Events",i+1);
		MCpeakByEtaP.push_back( new TH1D(name,title,30,60,120) );
		sprintf(name,"peaksEtaNRing%d",i+1);
		sprintf(title,"Gauss-fit Peaks, NT-, Eta Ring #%d;GF Peak;Events",i+1);
		MCpeakByEtaN.push_back( new TH1D(name,title,30,60,120) );
		sprintf(name,"meansEtaPRing%d",i+1);
		sprintf(title,"Mean Ratios, NT+, Eta Ring #%d;Mean Ratio;Events",i+1);
		MCmeanByEtaP.push_back( new TH1D(name,title,40,0,2) );
		sprintf(name,"meansEtaNRing%d",i+1);
		sprintf(title,"Mean Ratios, NT-, Eta Ring #%d;Mean Ratio;Events",i+1);
		MCmeanByEtaN.push_back( new TH1D(name,title,40,0,2) );	
	}
		
	//grab an MC ntuple
	TFile* f1 = new TFile("/local/cms/user/finkel/NoTrack/Ntuple/DE_MC2012_AlexTrgr3/DE_MC2012_AlexTrgr_newerCuts.root");	
	if(f1 == NULL)
	{
		std::cout<<"Failed to open MC file. Exiting."<<std::endl;
		return 1;
	}
	
	//make an ntuple:
	ZEffTree* zMC = new ZEffTree(*f1,false);
	//Fill the data histograms:
	double expectedPtMC;
	int ix,iy;
	const double Mz = 91.1876; //nominal Z mass
	for(int event=0; event<zMC->Entries(); event++ ) //fill the data hists
	{
		ix=zMC->reco.ix[1];
		iy=zMC->reco.iy[1];
		if( (fabs(zMC->reco.eta[1])>2.5) && (fabs(zMC->reco.eta[1])<3.0) 
			 &&(ix>29) &&(ix<71) &&(iy>29) &&(iy<71)
			 && zMC->reco.isSelected(1,"NTLooseElectronId-EtaDet") //using only events that pass selection now!
		  ) 
		{	
			if( (((ix-50)*(ix-50)+(iy-50)*(iy-50))  < 100) || (((ix-50)*(ix-50)+(iy-50)*(iy-50)) > 324) )
			{
				zMC->GetNextEvent();
				continue;
			}
			
			//"expected energy
			expectedPtMC = Mz*Mz / ( 2*zMC->reco.pt[1]*( cosh(zMC->reco.eta[1]-zMC->reco.eta[0]) - cos(zMC->reco.phi[1]-zMC->reco.phi[0]) ) );
			//NOTE: this is the energy of the entire CLUSTER, not just SEED!
			//Though most of it still comes from the seed crystal... so good enough for now.
			
			if(zMC->reco.eta[1]>0)//positive endcap
			{
				histMapMCP[std::make_pair(ix,iy)]->Fill(zMC->reco.mz);
				if( (zMC->reco.mz>75) && (zMC->reco.mz<105) )	ratioMapMCP[std::make_pair(ix,iy)]->Fill(expectedPtMC/zMC->reco.pt[1]);  //NOTE: Using Expected/Observed!
			}
			else//negative endcap
			{
				histMapMCN[std::make_pair(ix,iy)]->Fill(zMC->reco.mz);
				if( (zMC->reco.mz>75) && (zMC->reco.mz<105) )	ratioMapMCN[std::make_pair(ix,iy)]->Fill(expectedPtMC/zMC->reco.pt[1]);
			}
		}
		zMC->GetNextEvent();
	}
		
	//grab a data ntuple
	TFile* f2 = new TFile("/local/cms/user/finkel/NoTrack/Ntuple/DE_ReReco_NewerCuts/DE_ReReco2012FULL_Jan13_AlexTrgr_newerCuts_001.root");	
	if(f2 == NULL)
	{
		std::cout<<"Failed to open Data file. Exiting."<<std::endl;
		return 1;
	}
	
	//make data ntuple:
	ZEffTree* zDat = new ZEffTree(*f2,false);
	//Fill the data histograms:
	double expectedPtDat;
	for(int event=0; event<zDat->Entries(); event++ ) //fill the data hists
	{
		ix=zDat->reco.ix[1];
		iy=zDat->reco.iy[1];
		if( (fabs(zDat->reco.eta[1])>2.5) && (fabs(zDat->reco.eta[1])<3.0) 
			 &&(ix>29) &&(ix<71) &&(iy>29) &&(iy<71)
			 && zDat->reco.isSelected(1,"NTLooseElectronId-EtaDet") //using only events that pass selection now!
		  ) 
		{	
			if( (((ix-50)*(ix-50)+(iy-50)*(iy-50))  < 100) || (((ix-50)*(ix-50)+(iy-50)*(iy-50)) > 324) )
			{
				zDat->GetNextEvent();
				continue;
			}
			
			//"expected energy
			expectedPtDat = Mz*Mz / ( 2*zDat->reco.pt[1]*( cosh(zDat->reco.eta[1]-zDat->reco.eta[0]) - cos(zDat->reco.phi[1]-zDat->reco.phi[0]) ) );
			//NOTE: this is the energy of the entire CLUSTER, not just SEED!
			//Though most of it still comes from the seed crystal... so good enough for now.
			
			if(zDat->reco.eta[1]>0)//positive endcap
			{
				histMapDataP[std::make_pair(ix,iy)]->Fill(zDat->reco.mz);
				if( (zDat->reco.mz>75) && (zDat->reco.mz<105) )	ratioMapDataP[std::make_pair(ix,iy)]->Fill(expectedPtDat/zDat->reco.pt[1]);  //NOTE: Using Expected/Observed!
			} 
			else//negative endcap
			{
				histMapDataN[std::make_pair(ix,iy)]->Fill(zDat->reco.mz);
				if( (zDat->reco.mz>75) && (zDat->reco.mz<105) )	ratioMapDataN[std::make_pair(ix,iy)]->Fill(expectedPtDat/zDat->reco.pt[1]);
			}
		}
		zDat->GetNextEvent();
	}

	std::cout<<"Hists filled. Proceeding to makeing plots."<<std::endl;
	//gonna try fitting now!
	TF1 *ratioFit = new TF1("fit","gaus",0.4,1.6);
	TF1 *massFit = new TF1("fit","gaus",60,120);
	
	//hists to contain the means and the peaks:
	TH1D *allMeansDat = new TH1D("allMeansDat","Distribution of Means, all NT;mean",40,0.4,1.6);
	TH1D *meansDatN = new TH1D("meansDatN","Distribution of Means, NT-;mean",40,0.4,1.6);
	TH1D *meansDatP = new TH1D("meansDatP","Distribution of Means, NT+;mean",40,0.4,1.6);
	TH1D *allPeaksDat = new TH1D("allPeaksDat","Distribution of GF Peaks, all NT;mean",30,80,110);
	TH1D *peaksDatN = new TH1D("peaksDatN","Distribution of GF Peaks, NT-;mean",30,80,110);
	TH1D *peaksDatP = new TH1D("peaksDatP","Distribution of GF Peaks, NT+;mean",30,80,110);
	TH1D *ringAveragePeakDatN = new TH1D("ringAveragePeakDatN","NT- Ring-averaged GF Peaks;Eta Ring;GF Peak",8,0,8);
	TH1D *ringAveragePeakDatP = new TH1D("ringAveragePeakDatP","NT+ Ring-averaged GF Peaks;Eta Ring;GF Peak",8,0,8);
	TH1D *ringAverageMeanDatN = new TH1D("ringAverageMeanDatN","NT- Ring-averaged Mean Ratios",8,0,8);
	TH1D *ringAverageMeanDatP = new TH1D("ringAverageMeanDatP","NT+ Ring-averaged Mean Ratios",8,0,8);
	
	TH1D *allMeansMC = new TH1D("allMeansMC","",40,0.4,1.6);
	TH1D *meansMCN = new TH1D("meansMCN","",40,0.4,1.6);
	TH1D *meansMCP = new TH1D("meansMCP","",40,0.4,1.6);
	TH1D *allPeaksMC = new TH1D("allPeaksMC","",30,80,110);
	TH1D *peaksMCN = new TH1D("peaksMCN","",30,80,110);
	TH1D *peaksMCP = new TH1D("peaksMCP","",30,80,110);
	TH1D *ringAveragePeakMCN = new TH1D("ringAveragePeakMCN","",8,0,8);
	TH1D *ringAveragePeakMCP = new TH1D("ringAveragePeakMCP","",8,0,8);
	TH1D *ringAverageMeanMCN = new TH1D("ringAverageMeanMCN","",8,0,8);
	TH1D *ringAverageMeanMCP = new TH1D("ringAverageMeanMCP","",8,0,8);
	
	
	//MC-Data discrepancy maps:
	TH2D* meanDiscrepMapN = new TH2D("meanDiscrepMapN","#mu_{Data} - #mu_{MC} NT-;i_{x};i_{y}",40,30,70,40,30,70);
	TH2D* meanDiscrepMapP = new TH2D("meanDiscrepMapP","#mu_{Data} - #mu_{MC} NT+;i_{x};i_{y}",40,30,70,40,30,70);
	TH2D* peakDiscrepMapN = new TH2D("peakDiscrepMapN","GFPeak_{Data} - GFPeak_{MC} NT-;i_{x};i_{y}",40,30,70,40,30,70);
	TH2D* peakDiscrepMapP = new TH2D("peakDiscrepMapP","GFPeak_{Data} - GFPeak_{MC} NT+;i_{x};i_{y}",40,30,70,40,30,70);
	
	char filenameData[128];
	int ring;
	double datMean, datPeak;
	for(int i=30;i<71;i++)
	{
		for(int j=30;j<71;j++)
		{
			//if (massCount >= 10) break;
			if( (((i-50)*(i-50)+(j-50)*(j-50))  < 121) || (((i-50)*(i-50)+(j-50)*(j-50)) >= 324) )
			{	
				peakDiscrepMapP->SetBinContent(i-30,j-30,-1e9);
				peakDiscrepMapN->SetBinContent(i-30,j-30,-1e9);
				meanDiscrepMapP->SetBinContent(i-30,j-30,-1e9);
				meanDiscrepMapN->SetBinContent(i-30,j-30,-1e9);
				continue;			
			}
			ring = (int)sqrt( (i-50)*(i-50)+(j-50)*(j-50) )-11;
			
			if( histMapDataP[std::make_pair(i,j)]->Integral() > 50 )
			{				
				histMapDataP[std::make_pair(i,j)]->Fit(massFit,"QLM","",80,110);
				datPeak = massFit->GetParameter(1);
				DataMCpeakByEtaP[ring]->Fill(massFit->GetParameter(1));
				allPeaksDat->Fill(massFit->GetParameter(1));
				peaksDatP->Fill(massFit->GetParameter(1));
				histMapMCP[std::make_pair(i,j)]->Fit(massFit,"QML","",80,110);
				MCpeakByEtaP[ring]->Fill(massFit->GetParameter(1));
				allPeaksMC->Fill(massFit->GetParameter(1));
				peaksMCP->Fill(massFit->GetParameter(1));
				peakDiscrepMapP->SetBinContent(i-30,j-30,datPeak - (double)massFit->GetParameter(1));
			}
			if( histMapDataN[std::make_pair(i,j)]->Integral() > 50 )
			{
				histMapDataN[std::make_pair(i,j)]->Fit(massFit,"QLM","",80,110);
				datPeak = massFit->GetParameter(1);
				DataMCpeakByEtaN[ring]->Fill(massFit->GetParameter(1));
				allPeaksDat->Fill(massFit->GetParameter(1));
				peaksDatN->Fill(massFit->GetParameter(1));
				histMapMCN[std::make_pair(i,j)]->Fit(massFit,"QML","",80,110);
				MCpeakByEtaN[ring]->Fill(massFit->GetParameter(1));
				allPeaksMC->Fill(massFit->GetParameter(1));
				peaksMCN->Fill(massFit->GetParameter(1));
				peakDiscrepMapN->SetBinContent(i-30,j-30,datPeak - (double)massFit->GetParameter(1));
			} 
			if( ratioMapDataP[std::make_pair(i,j)]->Integral() > 50 )
			{
				ratioMapDataP[std::make_pair(i,j)]->Fit(ratioFit,"QLM","",0.4,1.6);
				datMean = ratioFit->GetParameter(1);
				allMeansDat->Fill(ratioFit->GetParameter(1));
				meansDatP->Fill(ratioFit->GetParameter(1));
				DataMCmeanByEtaP[ring]->Fill(ratioFit->GetParameter(1));
				ratioMapMCP[std::make_pair(i,j)]->Fit(ratioFit,"QML","",0.4,1.6);
				allMeansMC->Fill(ratioFit->GetParameter(1));
				meansMCP->Fill(ratioFit->GetParameter(1));
				MCmeanByEtaP[ring]->Fill(ratioFit->GetParameter(1));
				meanDiscrepMapP->SetBinContent(i-30,j-30,datMean - (double)ratioFit->GetParameter(1));
			}
			if( ratioMapDataN[std::make_pair(i,j)]->Integral() > 50 )
			{
				ratioMapDataN[std::make_pair(i,j)]->Fit(ratioFit,"QLM","",0.4,1.6);
				datMean = ratioFit->GetParameter(1);
				allMeansDat->Fill(ratioFit->GetParameter(1));
				meansDatN->Fill(ratioFit->GetParameter(1));
				DataMCmeanByEtaN[ring]->Fill(ratioFit->GetParameter(1));
				ratioMapMCN[std::make_pair(i,j)]->Fit(ratioFit,"QML","",0.4,1.6);
				allMeansMC->Fill(ratioFit->GetParameter(1));
				meansMCN->Fill(ratioFit->GetParameter(1));
				MCmeanByEtaN[ring]->Fill(ratioFit->GetParameter(1));
				meanDiscrepMapN->SetBinContent(i-30,j-30,datMean - (double)ratioFit->GetParameter(1));
			}
		}
	}
	std::cout<<"peakDiscrepMapP zero is "<<peakDiscrepMapP->GetBinContent(20,20)<<std::endl;
	std::cout<<"peakDiscrepMapN zero is "<<peakDiscrepMapN->GetBinContent(20,20)<<std::endl;
	std::cout<<"meanDiscrepMapP zero is "<<meanDiscrepMapP->GetBinContent(20,20)<<std::endl;
	std::cout<<"meanDiscrepMapN zero is "<<meanDiscrepMapN->GetBinContent(20,20)<<std::endl;
	
	//draw the comparison plots:
	allMeansMC->SetMarkerStyle(21);
	allMeansMC->SetMarkerColor(kRed);
	allMeansMC->Draw("P");
	allMeansDat->SetMarkerStyle(20);
	allMeansDat->SetMarkerColor(kBlue+1);
	allMeansDat->Draw("SAME P");
	l1->AddEntry(allMeansDat,"Data");
	l1->AddEntry(allMeansMC,"MC");
	l1->Draw();
	c1->Print("Efficiencies/CrystalCalibCombined/allMeans.png");
	c1->Clear();
	l1->Clear();

	allPeaksMC->SetMarkerStyle(21);
	allPeaksMC->SetMarkerColor(kRed);
	allPeaksMC->Draw("P");
	allPeaksDat->SetMarkerStyle(20);
	allPeaksDat->SetMarkerColor(kBlue+1);
	allPeaksDat->Draw("SAME P");
	l1->AddEntry(allPeaksDat,"Data");
	l1->AddEntry(allPeaksMC,"MC");
	l1->Draw();
	c1->Print("Efficiencies/CrystalCalibCombined/allPeaks.png");
	c1->Clear();
	l1->Clear();
	
	meansMCP->SetMarkerStyle(21);
	meansMCP->SetMarkerColor(kRed);
	meansMCP->Draw("P");
	meansDatP->SetMarkerStyle(20);
	meansDatP->SetMarkerColor(kBlue+1);
	meansDatP->Draw("SAME P");
	l1->AddEntry(meansDatP,"Data");
	l1->AddEntry(meansMCP,"MC");
	l1->Draw();
	c1->Print("Efficiencies/CrystalCalibCombined/meansP.png");
	c1->Clear();
	l1->Clear();
	
	peaksMCP->SetMarkerStyle(21);
	peaksMCP->SetMarkerColor(kRed);
	peaksMCP->Draw("P");
	peaksDatP->SetMarkerStyle(20);
	peaksDatP->SetMarkerColor(kBlue+1);
	peaksDatP->Draw("SAME P");
	l1->AddEntry(peaksDatP,"Data");
	l1->AddEntry(peaksMCP,"MC");
	l1->Draw();
	c1->Print("Efficiencies/CrystalCalibCombined/peaksP.png");
	c1->Clear();
	l1->Clear();
	
	meansMCN->SetMarkerStyle(21);
	meansMCN->SetMarkerColor(kRed);
	meansMCN->Draw("P");
	meansDatN->SetMarkerStyle(20);
	meansDatN->SetMarkerColor(kBlue+1);
	meansDatN->Draw("SAME P");
	l1->AddEntry(meansDatN,"Data");
	l1->AddEntry(meansMCN,"MC");
	l1->Draw();
	c1->Print("Efficiencies/CrystalCalibCombined/meansN.png");
	c1->Clear();
	l1->Clear();
	
	peaksMCN->SetMarkerStyle(21);
	peaksMCN->SetMarkerColor(kRed);
	peaksMCN->Draw("P");
	peaksDatN->SetMarkerStyle(20);
	peaksDatN->SetMarkerColor(kBlue+1);
	peaksDatN->Draw("SAME P");
	l1->AddEntry(peaksDatN,"Data");
	l1->AddEntry(peaksMCN,"MC");
	l1->Draw();
	c1->Print("Efficiencies/CrystalCalibCombined/peaksN.png");
	c1->Clear();
	l1->Clear();
	
	//discrepancy maps:
	peakDiscrepMapP->GetZaxis()->SetRangeUser(-10,10);
	peakDiscrepMapP->Draw("COLZ");
	c1->Print("Efficiencies/CrystalCalibCombined/DataMCPeakDiscrepMapP.png");
	c1->Clear();
	peakDiscrepMapN->GetZaxis()->SetRangeUser(-10,10);
	peakDiscrepMapN->Draw("COLZ");
	c1->Print("Efficiencies/CrystalCalibCombined/DataMCPeakDiscrepMapN.png");
	c1->Clear();
	meanDiscrepMapP->GetZaxis()->SetRangeUser(-0.4,0.4);
	meanDiscrepMapP->Draw("COLZ");
	c1->Print("Efficiencies/CrystalCalibCombined/DataMCMeanDiscrepMapP.png");
	c1->Clear();
	meanDiscrepMapN->GetZaxis()->SetRangeUser(-0.4,0.4);
	meanDiscrepMapN->Draw("COLZ");
	c1->Print("Efficiencies/CrystalCalibCombined/DataMCMeanDiscrepMapN.png");
	c1->Clear();
	

	
	for(int i=0;i<7;i++)
	{
		sprintf(filenameData,"Efficiencies/CrystalCalibCombined/MassPeakPEtaRing%d.png",i);
		MCpeakByEtaP[i]->SetMarkerStyle(21);
		MCpeakByEtaP[i]->SetMarkerColor(kRed);
		MCpeakByEtaP[i]->Draw("P");
		MCpeakByEtaP[i]->Fit(massFit,"QMLN","",80,100);
		massFit->SetLineWidth(0.5);
		massFit->SetLineColor(kRed);
		massFit->DrawCopy("same");
		ringAveragePeakMCP->SetBinContent(i+1,massFit->GetParameter(1));
		ringAveragePeakMCP->SetBinError(i+1,massFit->GetParError(1));
		DataMCpeakByEtaP[i]->SetMarkerStyle(20);
		DataMCpeakByEtaP[i]->SetMarkerColor(kBlue+1);
		DataMCpeakByEtaP[i]->Draw("SAME P");
		DataMCpeakByEtaP[i]->Fit(massFit,"QLMN","",80,100);
		massFit->SetLineColor(kBlue+1);
		massFit->DrawCopy("same");
		ringAveragePeakDatP->SetBinContent(i+1,massFit->GetParameter(1));
		ringAveragePeakDatP->SetBinError(i+1,massFit->GetParError(1));
		l1->AddEntry(DataMCpeakByEtaP[i],"Data");
		l1->AddEntry(MCpeakByEtaP[i],"MC");
		l1->Draw();
		c1->Print(filenameData);
		c1->Clear();
		l1->Clear();
		
		sprintf(filenameData,"Efficiencies/CrystalCalibCombined/MassPeakNEtaRing%d.png",i);
		MCpeakByEtaN[i]->SetMarkerStyle(21);
		MCpeakByEtaN[i]->SetMarkerColor(kRed);
		MCpeakByEtaN[i]->Draw("P");
		MCpeakByEtaN[i]->Fit(massFit,"QMLN","",80,100);
		massFit->SetLineWidth(0.5);
		massFit->SetLineColor(kRed);
		massFit->DrawCopy("same");
		ringAveragePeakMCN->SetBinContent(i+1,massFit->GetParameter(1));
		ringAveragePeakMCN->SetBinError(i+1,massFit->GetParError(1));
		DataMCpeakByEtaN[i]->SetMarkerStyle(20);
		DataMCpeakByEtaN[i]->SetMarkerColor(kBlue+1);
		DataMCpeakByEtaN[i]->Draw("SAME P");
		DataMCpeakByEtaN[i]->Fit(massFit,"QLMN","",80,100);
		massFit->SetLineColor(kBlue+1);
		massFit->DrawCopy("same");
		ringAveragePeakDatN->SetBinContent(i+1,massFit->GetParameter(1));
		ringAveragePeakDatN->SetBinError(i+1,massFit->GetParError(1));
		l1->AddEntry(DataMCpeakByEtaN[i],"Data");
		l1->AddEntry(MCpeakByEtaN[i],"MC");
		l1->Draw();
		c1->Print(filenameData);
		c1->Clear();
		l1->Clear();
		
		sprintf(filenameData,"Efficiencies/CrystalCalibCombined/MeanPEtaRing%d.png",i);
		MCmeanByEtaP[i]->SetMarkerStyle(21);
		MCmeanByEtaP[i]->SetMarkerColor(kRed);
		MCmeanByEtaP[i]->Draw("P");
		MCmeanByEtaP[i]->Fit(ratioFit,"QMLN","",0.4,1.6);
		ratioFit->SetLineWidth(0.5);
		ratioFit->SetLineColor(kRed);
		ratioFit->DrawCopy("same");
		ringAverageMeanMCP->SetBinContent(i+1,ratioFit->GetParameter(1));
		ringAverageMeanMCP->SetBinError(i+1,ratioFit->GetParError(1));
		DataMCmeanByEtaP[i]->SetMarkerStyle(20);
		DataMCmeanByEtaP[i]->SetMarkerColor(kBlue+1);
		DataMCmeanByEtaP[i]->Draw("SAME P");
		DataMCmeanByEtaP[i]->Fit(ratioFit,"QLMN","",0.4,1.6);
		ratioFit->SetLineColor(kBlue+1);
		ratioFit->DrawCopy("same");
		ringAverageMeanDatP->SetBinContent(i+1,ratioFit->GetParameter(1));
		ringAverageMeanDatP->SetBinError(i+1,ratioFit->GetParError(1));
		l1->AddEntry(DataMCmeanByEtaP[i],"Data");
		l1->AddEntry(MCmeanByEtaP[i],"MC");
		l1->Draw();
		c1->Print(filenameData);
		c1->Clear();
		l1->Clear();
		
		sprintf(filenameData,"Efficiencies/CrystalCalibCombined/MeanNEtaRing%d.png",i);
		MCmeanByEtaN[i]->SetMarkerStyle(21);
		MCmeanByEtaN[i]->SetMarkerColor(kRed);
		MCmeanByEtaN[i]->Draw("P");
		MCmeanByEtaN[i]->Fit(ratioFit,"QMLN","",0.4,1.6);
		ratioFit->SetLineWidth(0.5);
		ratioFit->SetLineColor(kRed);
		ratioFit->DrawCopy("same");
		ringAverageMeanMCN->SetBinContent(i+1,ratioFit->GetParameter(1));
		ringAverageMeanMCN->SetBinError(i+1,ratioFit->GetParError(1));
		DataMCmeanByEtaN[i]->SetMarkerStyle(20);
		DataMCmeanByEtaN[i]->SetMarkerColor(kBlue+1);
		DataMCmeanByEtaN[i]->Draw("SAME P");
		DataMCmeanByEtaN[i]->Fit(ratioFit,"QLMN","",0.4,1.6);
		ratioFit->SetLineColor(kBlue+1);
		ratioFit->DrawCopy("same");
		ringAverageMeanDatN->SetBinContent(i+1,ratioFit->GetParameter(1));
		ringAverageMeanDatN->SetBinError(i+1,ratioFit->GetParError(1));
		l1->AddEntry(DataMCmeanByEtaN[i],"Data");
		l1->AddEntry(MCmeanByEtaN[i],"MC");
		l1->Draw();
		c1->Print(filenameData);
		c1->Clear();
		l1->Clear();
	}
	
	ringAveragePeakDatP->SetMarkerStyle(20);
	ringAveragePeakDatP->SetMarkerColor(kBlue+1);
	ringAveragePeakDatP->GetYaxis()->SetRangeUser(85,105);
	ringAveragePeakDatP->Draw("E");
	ringAveragePeakMCP->SetMarkerStyle(21);
	ringAveragePeakMCP->SetMarkerColor(kRed);
	ringAveragePeakMCP->GetYaxis()->SetRangeUser(85,100);
	ringAveragePeakMCP->Draw("SAME E");
	l1->AddEntry(ringAveragePeakDatP,"Data");
	l1->AddEntry(ringAveragePeakMCP,"MC");
	l1->Draw();
	c1->Print("Efficiencies/CrystalCalibCombined/AveragePeakVsEtaRingP.png");
	c1->Clear();
	l1->Clear();
	
	ringAveragePeakDatN->SetMarkerStyle(20);
	ringAveragePeakDatN->SetMarkerColor(kBlue+1);
	ringAveragePeakDatN->GetYaxis()->SetRangeUser(85,105);
	ringAveragePeakDatN->Draw("E");
	ringAveragePeakMCN->SetMarkerStyle(21);
	ringAveragePeakMCN->SetMarkerColor(kRed);
	ringAveragePeakMCN->GetYaxis()->SetRangeUser(85,100);
	ringAveragePeakMCN->Draw("SAME E");
	l1->AddEntry(ringAveragePeakDatN,"Data");
	l1->AddEntry(ringAveragePeakMCN,"MC");
	l1->Draw();
	c1->Print("Efficiencies/CrystalCalibCombined/AveragePeakVsEtaRingN.png");
	c1->Clear();
	l1->Clear();
	
	ringAverageMeanDatP->SetMarkerStyle(20);
	ringAverageMeanDatP->SetMarkerColor(kBlue+1);
	ringAverageMeanDatP->GetYaxis()->SetRangeUser(0.8,1.2);
	ringAverageMeanDatP->Draw("E");
	ringAverageMeanMCP->SetMarkerStyle(21);
	ringAverageMeanMCP->SetMarkerColor(kRed);
	ringAverageMeanMCP->GetYaxis()->SetRangeUser(0.8,1.2);
	ringAverageMeanMCP->Draw("SAME E");
	l1->AddEntry(ringAverageMeanDatP,"Data");
	l1->AddEntry(ringAverageMeanMCP,"MC");
	l1->Draw();
	c1->Print("Efficiencies/CrystalCalibCombined/AverageMeanVsEtaRingP.png");
	c1->Clear();
	l1->Clear();	
	
	ringAverageMeanDatN->SetMarkerStyle(20);
	ringAverageMeanDatN->SetMarkerColor(kBlue+1);
	ringAverageMeanDatN->GetYaxis()->SetRangeUser(0.8,1.2);
	ringAverageMeanDatN->Draw("E");
	ringAverageMeanMCN->SetMarkerStyle(21);
	ringAverageMeanMCN->SetMarkerColor(kRed);
	ringAverageMeanMCN->GetYaxis()->SetRangeUser(0.8,1.2);
	ringAverageMeanMCN->Draw("SAME E");
	l1->AddEntry(ringAverageMeanDatN,"Data");
	l1->AddEntry(ringAverageMeanMCN,"MC");
	l1->Draw();
	c1->Print("Efficiencies/CrystalCalibCombined/AverageMeanVsEtaRingN.png");
	c1->Clear();
	l1->Clear();
		
	c1->Close();
	return 0;
}























