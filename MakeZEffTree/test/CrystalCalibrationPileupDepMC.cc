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

struct triplet
{
	public:
	int first, second, third;
	triplet(int i,int j, int k):first(i),second(j),third(k) {}
	bool operator<(const triplet lhs) const
	{
		if(lhs.first < first) return true;
		else if(lhs.first > first) return false;
		if(lhs.second < second) return true;
		else if(lhs.second > second) return false;
		if(lhs.third < third) return true;
		else if(lhs.third > third) return false;
		return false;
	}
};

int crystalCalibPileupDepMC()
{
	gROOT->SetStyle("Plain");	
	gStyle->SetErrorX(0);
	gStyle->SetOptStat(0);
	
	TCanvas* c1 = new TCanvas("C1", "c1", 1200, 900);
	c1->cd();

	TLegend* l1 = new TLegend(0.65,0.7,0.9,0.9);
	l1->SetFillColor(10);
	
	//let's be fancy and use a map!
	std::map< triplet, TH1D* > histMapP, histMapN, ratioMapP, ratioMapN;
	char nameMeanP[128], nameMeanN[128], nameRatioN[128], nameRatioP[128];	
	
	for(int i=30;i<71;i++)
	{
		for(int j=30;j<71;j++)
		{
			for(int k=0;k<4;k++)
			{
				//data
				if( (((i-50)*(i-50)+(j-50)*(j-50))  < 100) || (((i-50)*(i-50)+(j-50)*(j-50)) > 324) ) continue;
				sprintf(nameMeanP,"mZDatPU%dix%diy%dzP",k,i,j);
				sprintf(nameMeanN,"mZDatPU%dix%diy%dzN",k,i,j);
				sprintf(nameRatioP,"RatioDatPU%dIx%diy%dzP",k,i,j);
				sprintf(nameRatioN,"RatioDatPU%dIx%diy%dzN",k,i,j);
				histMapP[triplet(i,j,k)]= new TH1D(nameMeanP,"",30,60,120);
				histMapN[triplet(i,j,k)]= new TH1D(nameMeanN,"",30,60,120);
				ratioMapP[triplet(i,j,k)]= new TH1D(nameRatioP,"",20,0.4,1.6);
				ratioMapN[triplet(i,j,k)]= new TH1D(nameRatioN,"",20,0.4,1.6);
			}
		}
	}
		
	//getting real ambitious now: gonna do eta rings, with pileup!
	std::vector<TH1D*> peakByEtaP[4], peakByEtaN[4], meanByEtaP[4], meanByEtaN[4];
	char name[128], title[256];
	for(int i=0;i<7;i++)
	{
		for(int j=0;j<4;j++)
		{
			//data
			sprintf(name,"peaksDatPEtaRing%dPU%d",i+1,j);
			sprintf(title,"GF Peaks, Data, NT+, Eta Ring #%d, PU %d-%d;GF Peak;Events",i+1,10*j,10*(j+1));
			peakByEtaP[j].push_back( new TH1D(name,title,30,60,120) );
			sprintf(name,"peaksDatNEtaRing%dPU%d",i+1,j);
			sprintf(title,"GF Peaks, Data NT-, Eta Ring #%d, PU %d-%d;GF Peak;Events",i+1,10*j,10*(j+1));
			peakByEtaN[j].push_back( new TH1D(name,title,30,60,120) );
			sprintf(name,"meansDatPEtaRing%dPU%d",i+1,j);
			sprintf(title,"Mean Ratios, Data, NT+, Eta Ring #%d, PU %d-%d;Mean Ratio;Events",i+1,10*j,10*(j+1));
			meanByEtaP[j].push_back( new TH1D(name,title,40,0,2) );
			sprintf(name,"meansDatNEtaRing%dPU%d",i+1,j);
			sprintf(title,"Mean Ratios, Data, NT-, Eta Ring #%d, PU %d-%d;Mean Ratio;Events",i+1,10*j,10*(j+1));
			meanByEtaN[j].push_back( new TH1D(name,title,40,0,2) );
		}
	}
	std::vector<TH1D*> ringAveragePeaksP, ringAveragePeaksN,ringAverageMeansP, ringAverageMeansN;
	for(int i=0;i<4;i++)
	{
		sprintf(name,"ringAveragePeaksPPU%d",i);
		ringAveragePeaksP.push_back( new TH1D(name,"",8,0,8) );
		sprintf(name,"ringAveragePeaksNPU%d",i);
		ringAveragePeaksN.push_back( new TH1D(name,"",8,0,8) );
		sprintf(name,"ringAverageMeansPPU%d",i);
		ringAverageMeansP.push_back( new TH1D(name,"",8,0,8) );
		sprintf(name,"ringAverageMeansNPU%d",i);
		ringAverageMeansN.push_back( new TH1D(name,"",8,0,8) );
	}
		
	//grab a data ntuple
	TFile* f2 = new TFile("/local/cms/user/finkel/NoTrack/Ntuple/DE_MC2012_AlexTrgr3/DE_MC2012_AlexTrgr_newerCuts.root");	
	if(f2 == NULL)
	{
		std::cout<<"Failed to open Data file. Exiting."<<std::endl;
		return 1;
	}
	
	//make data ntuple:
	ZEffTree* zDat = new ZEffTree(*f2,false);
	//Fill the data histograms:
	double expectedPt;
	int ix,iy;
	int puBin;
	const double Mz = 91.1876; //nominal Z mass
	for(int event=0; event<zDat->Entries(); event++ ) //fill the data hists
	{
		ix=zDat->reco.ix[1];
		iy=zDat->reco.iy[1];
		if(zDat->reco.nverts>30 ) puBin=3;
		else puBin=(int)(zDat->reco.nverts/10);
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
			expectedPt = Mz*Mz / ( 2*zDat->reco.pt[1]*( cosh(zDat->reco.eta[1]-zDat->reco.eta[0]) - cos(zDat->reco.phi[1]-zDat->reco.phi[0]) ) );
			//NOTE: this is the energy of the entire CLUSTER, not just SEED!
			//Though most of it still comes from the seed crystal... so good enough for now.
			
			if(zDat->reco.eta[1]>0)//positive endcap
			{
				histMapP[triplet(ix,iy,puBin)]->Fill(zDat->reco.mz);
				if( (zDat->reco.mz>70) && (zDat->reco.mz<110) )	ratioMapP[triplet(ix,iy,puBin)]->Fill(expectedPt/zDat->reco.pt[1]);  //NOTE: Using Expected/Observed!
			} 	//NOTE: wider range to combat reduced statistics per PU slice!
			else//negative endcap
			{
				histMapN[triplet(ix,iy,puBin)]->Fill(zDat->reco.mz);
				if( (zDat->reco.mz>70) && (zDat->reco.mz<110) )	ratioMapN[triplet(ix,iy,puBin)]->Fill(expectedPt/zDat->reco.pt[1]);
			}
		}
		zDat->GetNextEvent();
	}

	std::cout<<"Hists filled. Proceeding to makeing plots."<<std::endl;
	//gonna try fitting now!
	TF1 *ratioFit = new TF1("ratioFit","gaus",0.4,1.6);
	TF1 *massFit = new TF1("massFit","gaus",60,120);

	int ring;
	for(int i=30;i<71;i++)
	{
		for(int j=30;j<71;j++)
		{
			if( (((i-50)*(i-50)+(j-50)*(j-50))  < 121) || (((i-50)*(i-50)+(j-50)*(j-50)) >= 324) ) continue;
			
			ring = (int)sqrt( (i-50)*(i-50)+(j-50)*(j-50) )-11;

			for(int k=0;k<4;k++)
			{
				std::cout<<histMapP[triplet(i,j,k)]->Integral()<<" events in PeaksP ix "<<i<<" iy "<<j<<" bin "<<k<<std::endl;
				if( histMapP[triplet(i,j,k)]->Integral() > 20 ) //again, lower stats, so lower requirements
				{	
					histMapP[triplet(i,j,k)]->Fit(massFit,"Q","",80,110);					
					peakByEtaP[k][ring]->Fill(massFit->GetParameter(1));
				}
				if( histMapN[triplet(i,j,k)]->Integral() > 20 )
				{
					histMapN[triplet(i,j,k)]->Fit(massFit,"Q","",80,110);
					peakByEtaN[k][ring]->Fill(massFit->GetParameter(1));
				} 
				if( ratioMapP[triplet(i,j,k)]->Integral() > 20 )
				{
					ratioMapP[triplet(i,j,k)]->Fit(ratioFit,"Q","",0.4,1.6);
					meanByEtaP[k][ring]->Fill(ratioFit->GetParameter(1));			
				}
				if( ratioMapN[triplet(i,j,k)]->Integral() > 20 )
				{
					ratioMapN[triplet(i,j,k)]->Fit(ratioFit,"Q","",0.4,1.6);	
					meanByEtaN[k][ring]->Fill(ratioFit->GetParameter(1));
				}
			}
		}
	}

	for(int b=0;b<4;b++)
	{
		for(int r=0;r<7;r++)
		{
			//std::cout<<"Bin "<<b<<", ring "<<r<<" has "<<peakByEtaP[b][r]->Integral()<<" events in NT+ Peaks"<<std::endl;
			if(peakByEtaP[b][r]->Integral() > 10)
			{
				//peakByEtaP[b][r]->Fit(massFit,"Q","",80,100);
				//ringAveragePeaksP[b]->SetBinContent(r+1,massFit->GetParameter(1));
				//ringAveragePeaksP[b]->SetBinError(r+1,massFit->GetParError(1));
				ringAveragePeaksP[b]->SetBinContent(r+1,peakByEtaP[b][r]->GetMean(1));
				ringAveragePeaksP[b]->SetBinError(r+1,peakByEtaP[b][r]->GetMeanError(1));
			}
			//else std::cout<<"Insufficienct Data for NT+ Peaks"
			if(peakByEtaN[b][r]->Integral() > 10)
			{
				//peakByEtaN[b][r]->Fit(massFit,"Q","",80,100);
				//ringAveragePeaksN[b]->SetBinContent(r+1,massFit->GetParameter(1));
				//ringAveragePeaksN[b]->SetBinError(r+1,massFit->GetParError(1));
				ringAveragePeaksN[b]->SetBinContent(r+1,peakByEtaN[b][r]->GetMean(1));
				ringAveragePeaksN[b]->SetBinError(r+1,peakByEtaN[b][r]->GetMeanError(1));
			}
			if(meanByEtaP[b][r]->Integral() > 10)
			{
				//meanByEtaP[b][r]->Fit(ratioFit,"Q","",0.4,1.6);
				//ringAverageMeansP[b]->SetBinContent(r+1,ratioFit->GetParameter(1));
				//ringAverageMeansP[b]->SetBinError(r+1,ratioFit->GetParError(1));
				ringAverageMeansP[b]->SetBinContent(r+1,meanByEtaP[b][r]->GetMean(1));
				ringAverageMeansP[b]->SetBinError(r+1,meanByEtaP[b][r]->GetParError(1));
			}
			if(meanByEtaN[b][r]->Integral() > 10)
			{
				//meanByEtaN[b][r]->Fit(ratioFit,"Q","",0.4,1.6);
				//ringAverageMeansN[b]->SetBinContent(r+1,ratioFit->GetParameter(1));
				//ringAverageMeansN[b]->SetBinError(r+1,ratioFit->GetParError(1));
				ringAverageMeansN[b]->SetBinContent(r+1,meanByEtaN[b][r]->GetMean(1));
				ringAverageMeansN[b]->SetBinError(r+1,meanByEtaN[b][r]->GetParError(1));
			}
		}
	}

	int marker = 20;
	int color = 4;
	char legend[128];
	ringAveragePeaksP[0]->SetTitle("GF Peak vs. #eta-ring, NT+;#eta-ring;<Z-peak>");
	sprintf(legend,"0<nVert<10");
	ringAveragePeaksP[0]->SetMarkerStyle(marker);
	ringAveragePeaksP[0]->SetMarkerColor(color);
	ringAveragePeaksP[0]->GetYaxis()->SetRangeUser(85,100);
	l1->AddEntry(ringAveragePeaksP[0],legend);
	ringAveragePeaksP[0]->Draw("E");
	for(int b=1;b<3;b++)
	{
		marker++;
		color--;
		sprintf(legend,"%d<nVert<%d",10*b,10*(b+1));
		ringAveragePeaksP[b]->SetMarkerStyle(marker);
		ringAveragePeaksP[b]->SetMarkerColor(color);
		l1->AddEntry(ringAveragePeaksP[b],legend);
		ringAveragePeaksP[b]->Draw("same E");
	}
	l1->Draw();
	c1->Print("Efficiencies/CrystalCalibPileupDep_MC/RingAveragePeaksVsEtaP.png");
	c1->Clear();
	l1->Clear();
	
	marker = 20;
	color = 4;
	ringAveragePeaksN[0]->SetTitle("GF Peak vs. #eta-ring, NT-;#eta-ring;<Z-peak>");
	sprintf(legend,"0<nVert<10");
	ringAveragePeaksN[0]->SetMarkerStyle(marker);
	ringAveragePeaksN[0]->SetMarkerColor(color);
	ringAveragePeaksN[0]->GetYaxis()->SetRangeUser(85,100);
	l1->AddEntry(ringAveragePeaksN[0],legend);
	ringAveragePeaksN[0]->Draw("E");
	for(int b=1;b<3;b++)
	{
		marker++;
		color--;
		sprintf(legend,"%d<nVert<%d",10*b,10*(b+1));
		ringAveragePeaksN[b]->SetMarkerStyle(marker);
		ringAveragePeaksN[b]->SetMarkerColor(color);
		l1->AddEntry(ringAveragePeaksN[b],legend);
		ringAveragePeaksN[b]->Draw("same E");
	}
	l1->Draw();
	c1->Print("Efficiencies/CrystalCalibPileupDep_MC/RingAveragePeaksVsEtaN.png");
	c1->Clear();
	l1->Clear();
	
	marker = 20;
	color = 4;
	ringAverageMeansP[0]->SetTitle("Mean Ratio vs. #eta-ring, NT+;#eta-ring;<E_{exp}/E_{obs}>");
	sprintf(legend,"0<nVert<10");
	ringAverageMeansP[0]->SetMarkerStyle(marker);
	ringAverageMeansP[0]->SetMarkerColor(color);
	ringAverageMeansP[0]->GetYaxis()->SetRangeUser(0.8,1.05);
	l1->AddEntry(ringAverageMeansP[0],legend);
	ringAverageMeansP[0]->Draw("E");
	for(int b=1;b<3;b++)
	{
		marker++;
		color--;
		sprintf(legend,"%d<nVert<%d",10*b,10*(b+1));
		ringAverageMeansP[b]->SetMarkerStyle(marker);
		ringAverageMeansP[b]->SetMarkerColor(color);
		l1->AddEntry(ringAverageMeansP[b],legend);
		ringAverageMeansP[b]->Draw("same E");
	}
	l1->Draw();
	c1->Print("Efficiencies/CrystalCalibPileupDep_MC/RingAverageMeansVsEtaP.png");
	c1->Clear();
	l1->Clear();
	
	marker = 20;
	color = 4;
	ringAverageMeansN[0]->SetTitle("Mean Ratio vs. #eta-ring, NT-;#eta-ring;<E_{exp}/E_{obs}>");
	sprintf(legend,"0<nVert<10");
	ringAverageMeansN[0]->SetMarkerStyle(marker);
	ringAverageMeansN[0]->SetMarkerColor(color);
	ringAverageMeansN[0]->GetYaxis()->SetRangeUser(0.8,1.05);
	l1->AddEntry(ringAverageMeansN[0],legend);
	ringAverageMeansN[0]->Draw("E");
	for(int b=1;b<3;b++)
	{
		marker++;
		color--;
		sprintf(legend,"%d<nVert<%d",10*b,10*(b+1));
		ringAverageMeansN[b]->SetMarkerStyle(marker);
		ringAverageMeansN[b]->SetMarkerColor(color);
		l1->AddEntry(ringAverageMeansN[b],legend);
		ringAverageMeansN[b]->Draw("same E");
	}
	l1->Draw();
	c1->Print("Efficiencies/CrystalCalibPileupDep_MC/RingAverageMeansVsEtaN.png");
	c1->Clear();
	l1->Clear();

	c1->Close();

	return 0;
}
