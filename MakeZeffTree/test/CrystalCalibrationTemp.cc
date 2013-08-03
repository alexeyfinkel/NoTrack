//this script is going to attempt making Z mass plots for each NT xtal.


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


#include "../src/ZEffTree.h"

int crystalCalibTemp()
{
	gROOT->SetStyle("Plain");	
	gStyle->SetErrorX(0);
	gStyle->SetOptStat(0);
	
	TCanvas* c1 = new TCanvas("C1", "c1", 1200, 900);
	c1->cd();
	//c1->SetLeftMargin(0.08);
	//c1->SetRightMargin(0.12);
	
	TLegend* l1 = new TLegend(0.65,0.7,0.9,0.9);
	l1->SetFillColor(10);
	
	//let's be fancy and use a map!
	std::map< pair<int,int>, TH1D* > histMapP, histMapN, ratioMapP, ratioMapN;
	char nameP[128], titleP[256],nameN[128], titleN[256],nameRN[128], titleRN[256], nameRP[128], titleRP[256];
	for(int i=30;i<71;i++)
	{
		for(int j=30;j<71;j++)
		{
			if( (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5))  < 132.25) || (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5)) > 342.25) ) continue;
			sprintf(nameP,"mZix%diy%dzP",i,j);
			sprintf(titleP,"M_{ee}, NT+ i_{x} = %d, i_{y} = %d;M_{ee};Events/2GeV",i,j);
			sprintf(nameN,"mZix%diy%dzN",i,j);
			sprintf(titleN,"M_{ee}, NT- i_{x} = %d, i_{y} = %d;M_{ee};Events/2GeV",i,j);
			sprintf(nameRP,"RatioIx%diy%dzP",i,j);
			sprintf(titleRP,"E_{exp}/E_{obs}, NT+ i_{x} = %d, i_{y} = %d;E_{exp}/E_{obs};Events",i,j);
			sprintf(nameRN,"RatioIx%diy%dzN",i,j);
			sprintf(titleRN,"E_{exp}/E_{obs}, NT- i_{x} = %d, i_{y} = %d;E_{exp}/E_{obs};Events",i,j);
			histMapP[std::make_pair(i,j)]= new TH1D(nameP,titleP,30,60,120);
			histMapN[std::make_pair(i,j)]= new TH1D(nameN,titleN,30,60,120);
			ratioMapP[std::make_pair(i,j)]= new TH1D(nameRP,titleRP,20,0.4,1.6);
			ratioMapN[std::make_pair(i,j)]= new TH1D(nameRN,titleRN,20,0.4,1.6);
		}
	}
	std::cout<<"Map made; proceeding to make ntuple."<<std::endl;
		
	//getting real ambitious now: gonna do eta rings!
	std::vector<TH1D*> peakByEtaP, peakByEtaN, meanByEtaP, meanByEtaN;
	char ringName[128], ringTitle[256];
	for(int i=0;i<7;i++)
	{
		sprintf(ringName,"peaksEtaPRing%d",i+1);
		sprintf(ringTitle,"Gauss-fit Peaks, NT+, Eta Ring #%d;GF Peak;Events",i+1);
		peakByEtaP.push_back( new TH1D(ringName,ringTitle,30,60,120) );
		sprintf(ringName,"peaksEtaNRing%d",i+1);
		sprintf(ringTitle,"Gauss-fit Peaks, NT-, Eta Ring #%d;GF Peak;Events",i+1);
		peakByEtaN.push_back( new TH1D(ringName,ringTitle,30,60,120) );
		sprintf(ringName,"meansEtaPRing%d",i+1);
		sprintf(ringTitle,"Mean Ratios, NT+, Eta Ring #%d;Mean Ratio;Events",i+1);
		meanByEtaP.push_back( new TH1D(ringName,ringTitle,40,0,2) );
		sprintf(ringName,"meansEtaNRing%d",i+1);
		sprintf(ringTitle,"Mean Ratios, NT-, Eta Ring #%d;Mean Ratio;Events",i+1);
		meanByEtaN.push_back( new TH1D(ringName,ringTitle,40,0,2) );		
	}	
	
	//file to store calibration constants
	ofstream calFile;
	calFile.open("CalConstsMeanRatioData.txt",ios::trunc);
	if(!calFile.is_open())
	{
		std::cout<<"Failed to create cal. constants file. Existing"<<std::endl;
		return 1;
	}
	//st precision
	calFile<<std::setprecision(7)<<std::fixed;
			
	//grab a data ntuple
	TFile* f2 = new TFile("/afs/cern.ch/user/a/afinkel/public/NoTrack/CMSSW_5_3_8_patch3/src/NoTrack/MakeZeffTree/test_newCuts.root");	
	if(f2 == NULL)
	{
		std::cout<<"Failed to open Data file. Exiting."<<std::endl;
		return 1;
	}
	
	//make data ntuple:
	ZEffTree* ze2 = new ZEffTree(*f2,false);
	//Fill the data histograms:
	int ix,iy;
	const double Mz = 91.1876; //nominal Z mass
	double expectedPt;
	for(int event=0; event<ze2->Entries(); event++ ) //fill the data hists
	{
		ix=ze2->reco.ix[1];
		iy=ze2->reco.iy[1];
		if( (fabs(ze2->reco.eta[1])>2.5) && (fabs(ze2->reco.eta[1])<3.0) 
			 &&(ix>29) &&(ix<71) &&(iy>29) &&(iy<71)
			 && ze2->reco.isSelected(1,"NTLooseElectronId-EtaDet") //using only events that pass selection now!
		  ) 
		{
				std::cout<<"Ping!"<<std::endl;	
			if( (((ix-50.5)*(ix-50.5)+(iy-50.5)*(iy-50.5))  < 132.25) || (((ix-50.5)*(ix-50.5)+(iy-50.5)*(iy-50.5)) > 342.25) )
			{
				ze2->GetNextEvent();
				continue;
			}
			
			//"expected energy
			expectedPt = Mz*Mz / ( 2*ze2->reco.pt[1]*( cosh(ze2->reco.eta[1]-ze2->reco.eta[0]) - cos(ze2->reco.phi[1]-ze2->reco.phi[0]) ) );
			//NOTE: this is the energy of the entire CLUSTER, not just SEED!
			//Though most of it still comes from the seed crystal... so good enough for now.
			
			if(ze2->reco.eta[1]>0)//positive endcap
			{
				histMapP[std::make_pair(ix,iy)]->Fill(ze2->reco.mz);
				if( (ze2->reco.mz>75) && (ze2->reco.mz<105) )
				{
					std::cout<<"\n\nNEXT HIT ("<<ze2->ixs->size()<<" xtl hits):"<<std::endl;
					for(unsigned int k=0; k<(ze2->ixs->size()); k++ )
					{
						
						std::cout<<"Xtal ("<<(ze2->ixs->at(k))<<","<<(ze2->iys->at(k))<<"), fraction = "<<(ze2->hitEnergyFractions->at(k))<<std::endl;
						//ratioMapP[std::make_pair(hitIx,hitIy)]->Fill(expectedPt/ze2->reco.pt[1],ze2->hitEnergyFractions);  //NOTE: Using Expected/Observed!
					}
				}
			}
			else//negative endcap
			{
				histMapN[std::make_pair(ix,iy)]->Fill(ze2->reco.mz);
				if( (ze2->reco.mz>75) && (ze2->reco.mz<105) )
				{
					std::cout<<"\n\nNEXT HIT ("<<ze2->ixs->size()<<" xtl hits):"<<std::endl;
					for(unsigned int k=0; k<(ze2->ixs->size()); k++ )
					{
						
						std::cout<<"Xtal ("<<(ze2->ixs->at(k))<<","<<(ze2->iys->at(k))<<"), fraction = "<<(ze2->hitEnergyFractions->at(k))<<std::endl;
						//ratioMapN[std::make_pair(hitIx,hitIy)]->Fill(expectedPt/ze2->reco.pt[1],ze2->hitEnergyFractions);  //NOTE: Using Expected/Observed!
					}
				}
			}
		}
		ze2->GetNextEvent();
	}

	std::cout<<"Hists filled. Proceeding to makeing plots."<<std::endl;
	//gonna try fitting now!
	/*TF1 *ratioFit = new TF1("fit","gaus",0.4,1.6);
	TF1 *massFit = new TF1("fit","gaus",60,120);
	
	//hists to contain the means and the peaks:
	TH1D *allMeans = new TH1D("allMeans","Distribution of Means, all NT;mean",40,0,2);
	TH1D *meansN = new TH1D("meansN","Distribution of Means, NT-;mean",40,0,2);
	TH1D *meansP = new TH1D("meansP","Distribution of Means, NT+;mean",40,0,2);
	TH2D *meansMapN = new TH2D("meansMapN","Mean, NT-;i_{x};i_{y}",40,30,70,40,30,70);
	TH2D *meansMapP = new TH2D("meansMapP","Mean, NT+;i_{x};i_{y}",40,30,70,40,30,70);
	TH2D *peaksMapN = new TH2D("peaksMapN","Gauss-fit Z peak, NT-;i_{x};i_{y}",40,30,70,40,30,70);
	TH2D *peaksMapP = new TH2D("peaksMapP","Gauss-fit Z peak, NT+;i_{x};i_{y}",40,30,70,40,30,70);
	TH1D *ringAveragePeakN = new TH1D("ringAveragePeaksN","Ring-averaged GF Peaks, NT-;#eta-ring;GF Peak",8,0,8);
	TH1D *ringAveragePeakP = new TH1D("ringAveragePeaksP","Ring-averaged GF Peaks, NT+;#eta-ring;GF Peak",8,0,8);
	TH1D *ringAverageMeanN = new TH1D("ringAverageMeansN","Ring-averaged Mean Ratios, NT-;#eta-ring;<E_{exp}/E_{obs}>",8,0,8);
	TH1D *ringAverageMeanP = new TH1D("ringAverageMeansP","Ring-averaged Mean Ratios, NT+;#eta-ring;<E_{exp}/E_{obs}>",8,0,8);
	
	//std::cout<<"Ping!"<<std::endl;
	
	//char filename[128];//,legend1[256],legend2[256];
	int massCount=0;
	int ratioCount=0;
	int ring;
	for(int i=30;i<71;i++)
	{
		for(int j=30;j<71;j++)
		{
			//if (massCount >= 10) break;
			if( (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5))  < 132.25) || (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5)) > 342.25) ) continue;
			
			ring = (int)(sqrt( (i-50.5)*(i-50.5)+(j-50.5)*(j-50.5) )-11.5);
			
			if( ring<0 || ring>7)
			{
				std::cout<<"Invalid Ring number: ring "<<ring<<std::endl;
			}
			
			if( histMapP[std::make_pair(i,j)]->Integral() > 50 )
			{				
				histMapP[std::make_pair(i,j)]->Fit(massFit,"QLMN","",80,110);
				
				peaksMapP->SetBinContent(i-30,j-30,massFit->GetParameter(1));
				peakByEtaP[ring]->Fill(massFit->GetParameter(1));
				massCount++;
			}
			if( histMapN[std::make_pair(i,j)]->Integral() > 50 )
			{
				histMapN[std::make_pair(i,j)]->Fit(massFit,"QLMN","",80,110);
				peaksMapN->SetBinContent(i-30,j-30,massFit->GetParameter(1));
				peakByEtaN[ring]->Fill(massFit->GetParameter(1));
				massCount++;
			} 
			if( ratioMapP[std::make_pair(i,j)]->Integral() > 50 )
			{
				ratioMapP[std::make_pair(i,j)]->Fit(ratioFit,"QLMN","",0.4,1.6);
				allMeans->Fill(ratioFit->GetParameter(1));
				meansP->Fill(ratioFit->GetParameter(1));
				meansMapP->SetBinContent(i-30,j-30,ratioFit->GetParameter(1));
				meanByEtaP[ring]->Fill(ratioFit->GetParameter(1));
				calFile<<i<<"\t"<<j<<"\t1\t"<<(float)ratioFit->GetParameter(1)<<"\t"<<(float)ratioFit->GetParError(1)<<"\n";
				ratioCount++;				
			}
			else
			{
				calFile<<i<<"\t"<<j<<"\t1\t-1\t999\n";
			}
			if( ratioMapN[std::make_pair(i,j)]->Integral() > 50 )
			{
				ratioMapN[std::make_pair(i,j)]->Fit(ratioFit,"QLMN","",0.4,1.6);
				allMeans->Fill(ratioFit->GetParameter(1));
				meansN->Fill(ratioFit->GetParameter(1));
				meansMapN->SetBinContent(i-30,j-30,ratioFit->GetParameter(1));	
				meanByEtaN[ring]->Fill(ratioFit->GetParameter(1));
				calFile<<i<<"\t"<<j<<"\t-1\t"<<(float)ratioFit->GetParameter(1)<<"\t"<<(float)ratioFit->GetParError(1)<<"\n";
				ratioCount++;
			}
			else
			{
				calFile<<i<<"\t"<<j<<"\t-1\t-1\t999\n";
			}
		}
	}
	calFile.close();
	std::cout<<"Looped thought all the crystals."<<std::endl;
	
	//draw histograms:
	allMeans->SetMarkerStyle(20);
	allMeans->Draw("E");
	c1->Print("Efficiencies/CrystalCalibData/allMeans.png");
	c1->Clear();
	meansP->SetMarkerStyle(20);
	meansP->Draw("E");
	c1->Print("Efficiencies/CrystalCalibData/MeansP.png");
	c1->Clear();
	meansN->SetMarkerStyle(20);
	meansN->Draw("E");
	c1->Print("Efficiencies/CrystalCalibData/MeansN.png");
	c1->Clear();
	meansMapP->GetZaxis()->SetRangeUser(0.6,1.4);
	meansMapP->GetZaxis()->SetLabelSize(0.03);
	meansMapP->Draw("colz");
	c1->Print("Efficiencies/CrystalCalibData/MeansMapP.png");
	c1->Clear();
	meansMapN->GetZaxis()->SetRangeUser(0.6,1.4);
	meansMapN->GetZaxis()->SetLabelSize(0.03);
	meansMapN->Draw("colz");
	c1->Print("Efficiencies/CrystalCalibData/MeansMapN.png");
	c1->Clear();
	peaksMapP->GetZaxis()->SetRangeUser(80,105);
	peaksMapP->GetZaxis()->SetLabelSize(0.03);
	peaksMapP->Draw("colz");
	c1->Print("Efficiencies/CrystalCalibData/MassPeaksMapP.png");
	c1->Clear();
	peaksMapN->GetZaxis()->SetRangeUser(80,105);
	peaksMapN->GetZaxis()->SetLabelSize(0.03);
	peaksMapN->Draw("colz");
	c1->Print("Efficiencies/CrystalCalibData/MassPeaksMapN.png");
	c1->Clear();	
	
	for(int i=0;i<7;i++)
	{
		sprintf(filename,"Efficiencies/CrystalCalibData/MassPeakPEtaRing%d.png",i);
		peakByEtaP[i]->SetMarkerStyle(20);
		peakByEtaP[i]->SetMarkerColor(kBlue+2);
		peakByEtaP[i]->Draw("E");
		//peakByEtaP[i]->Fit(massFit,"Q","",80,100);
		//ringAveragePeakP->SetBinContent(i+1,massFit->GetParameter(1));
		//ringAveragePeakP->SetBinError(i+1,massFit->GetParError(1));
		ringAveragePeakP->SetBinContent(i+1,peakByEtaP[i]->GetMean(1));
		ringAveragePeakP->SetBinError(i+1,peakByEtaP[i]->GetMeanError(1));
		c1->Print(filename);
		c1->Clear();
		sprintf(filename,"Efficiencies/CrystalCalibData/MassPeakNEtaRing%dN.png",i);
		peakByEtaN[i]->SetMarkerStyle(20);
		peakByEtaN[i]->SetMarkerColor(kBlue+2);
		peakByEtaN[i]->Draw("E");
		//peakByEtaN[i]->Fit(massFit,"Q","",80,100);
		//ringAveragePeakN->SetBinContent(i+1,massFit->GetParameter(1));
		//ringAveragePeakN->SetBinError(i+1,massFit->GetParError(1));
		ringAveragePeakN->SetBinContent(i+1,peakByEtaN[i]->GetMean(1));
		ringAveragePeakN->SetBinError(i+1,peakByEtaN[i]->GetMeanError(1));
		c1->Print(filename);
		c1->Clear();
		sprintf(filename,"Efficiencies/CrystalCalibData/MeanPEtaRing%d.png",i);
		meanByEtaP[i]->SetMarkerStyle(20);
		meanByEtaP[i]->SetMarkerColor(kBlue+2);
		meanByEtaP[i]->Draw("E");
		//meanByEtaP[i]->Fit(ratioFit,"Q","",0.4,1.6);
		//ringAverageMeanP->SetBinContent(i+1,ratioFit->GetParameter(1));
		//ringAverageMeanP->SetBinError(i+1,ratioFit->GetParError(1));
		ringAverageMeanP->SetBinContent(i+1,meanByEtaP[i]->GetMean(1));
		ringAverageMeanP->SetBinError(i+1,meanByEtaP[i]->GetMeanError(1));
		c1->Print(filename);
		c1->Clear();
		sprintf(filename,"Efficiencies/CrystalCalibData/MeanNEtaRing%dN.png",i);
		meanByEtaN[i]->SetMarkerStyle(20);
		meanByEtaN[i]->SetMarkerColor(kBlue+2);
		meanByEtaN[i]->Draw("E");
		//meanByEtaN[i]->Fit(ratioFit,"Q","",0.4,1.6);
		//ringAverageMeanN->SetBinContent(i+1,ratioFit->GetParameter(1));
		//ringAverageMeanN->SetBinError(i+1,ratioFit->GetParError(1));
		ringAverageMeanN->SetBinContent(i+1,meanByEtaN[i]->GetMean(1));
		ringAverageMeanN->SetBinError(i+1,meanByEtaN[i]->GetMeanError(1));
		c1->Print(filename);
		c1->Clear();
	}
	
	ringAveragePeakP->SetMarkerStyle(20);
	ringAveragePeakP->SetMarkerSize(2);
	ringAveragePeakP->SetMarkerColor(kBlue+1);
	ringAveragePeakP->GetYaxis()->SetRangeUser(85,100);
	ringAveragePeakP->Draw("E");
	c1->Print("Efficiencies/CrystalCalibData/AveragePeakVsEtaRingP.png");
	c1->Clear();
	ringAveragePeakN->SetMarkerStyle(20);
	ringAveragePeakN->SetMarkerSize(2);
	ringAveragePeakN->SetMarkerColor(kBlue+1);
	ringAveragePeakN->GetYaxis()->SetRangeUser(85,100);
	ringAveragePeakN->Draw("E");
	c1->Print("Efficiencies/CrystalCalibData/AveragePeakVsEtaRingN.png");
	c1->Clear();
	ringAverageMeanP->SetMarkerStyle(20);
	ringAverageMeanP->SetMarkerSize(2);
	ringAverageMeanP->SetMarkerColor(kBlue+1);
	ringAverageMeanP->GetYaxis()->SetRangeUser(0.8,1.2);
	ringAverageMeanP->Draw("E");
	c1->Print("Efficiencies/CrystalCalibData/AverageMeanVsEtaRingP.png");
	c1->Clear();
	ringAverageMeanN->SetMarkerStyle(20);
	ringAverageMeanN->SetMarkerSize(2);
	ringAverageMeanN->SetMarkerColor(kBlue+1);
	ringAverageMeanN->GetYaxis()->SetRangeUser(0.8,1.2);
	ringAverageMeanN->Draw("E");
	c1->Print("Efficiencies/CrystalCalibData/AverageMeanVsEtaRingN.png");
	c1->Clear();*/
		
	c1->Close();
	return 0;
}
















