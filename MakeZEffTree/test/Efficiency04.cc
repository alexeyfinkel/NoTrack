// this is an attempt at getting the NoTrack selection efficiency relative to Alex's trigger for trackless electrons
//this version also attempts to do binning by nvertex and probe pt
//this is version 4, which tries to use Root's built-in stat. uncertainty calculation and the ix, iy from new Ntuples
//also, using MC reco, instead of gen, like supposed to!

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include "TGraph.h"
#include <TLorentzVector.h>
#include "TRandom3.h"
#include <math.h>
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include <TMath.h>
#include "TGraphAsymmErrors.h"
#include "TGraphPainter.h"

#include "../src/ZEffTree.h"
#include "EfficiencyStatistics.h"

struct fitStruct
{	
	TH1* datHist;
	fitStruct(TH1*, TH1*);
	TF1 *fitFunct, *bkgdFunct;
	void doFit();
	double getSignalEvt1();
	double getSignalEvt2();
	
	double datEvts, bkgdEvts;
	
	private:
	struct fitFunction
	{
		TH1* sigHist;
		
	    fitFunction(TH1* sh)
	    {
	    	sigHist = (TH1*)sh->Clone();
	        sigHist->Scale(1/sigHist->Integral());
	    }
		Double_t fun(Double_t *x,Double_t *par);
		
		double operator()(double *x, double *par)
		{
			return fun(x, par);
		}
	} fit;
	
	struct bkdgFunction
	{
		Double_t bkgdFun(Double_t *x,Double_t *par);
		
		double operator()(double *x, double *par)
		{
			return bkgdFun(x, par);
		}
	} bkgd;
	
	double integral(TF1*, int ,int);
};

fitStruct::fitStruct(TH1* signal, TH1* data) : fit(signal)
{
	datHist = (TH1*)data->Clone();
	//TF1* fitFunct = new TF1("fitFunct", "fun", 0, 150);
}

Double_t fitStruct::fitFunction::fun(Double_t *x,Double_t *par)
{
	//return exp(par[0]-x[0]*par[1]) + par[2]*sigHist->GetBinContent(sigHist->FindBin(x[0]) + par[3]);
	return ( TMath::Erfc((par[0]-x[0])/par[1])*exp(par[2]-x[0]/par[3]) + par[4]*sigHist->GetBinContent(sigHist->FindBin(x[0])) );// + par[5] );
}

Double_t fitStruct::bkdgFunction::bkgdFun(Double_t *x,Double_t *par)
{
	return ( TMath::Erfc((par[0]-x[0])/par[1])*exp(par[2]-x[0]/par[3]) );// +par[4]);
}

void fitStruct::doFit()
{
	fitFunct = new TF1("fitFunct", fit, 60, 150, 5);
	fitFunct->SetParLimits(0,50,81);
	fitFunct->SetParLimits(1,2,25);
	fitFunct->SetParLimits(2,5,50);
	fitFunct->SetParLimits(3,20,90);
	fitFunct->SetParLimits(4,1,1e9);
	//fitFunct->SetParLimits(5,0,100);
	fitFunct->SetParameter(0,50);
	fitFunct->SetParameter(1,2);
	fitFunct->SetParameter(2,1);
	fitFunct->SetParameter(3,20);
	fitFunct->SetParameter(4,10000);
	//fitFunct->SetParameter(5,10);

	char xTitle[100];

	datHist->GetXaxis()->SetRangeUser(60, 150);
	//datHist->SetTitle(";M_{ee} (GeV);Events");
	datHist->SetMarkerStyle(20);	
	datHist->Fit("fitFunct","RML","E");
	sprintf(xTitle,"M_{ee} (GeV), #chi^{2}=%f",fitFunct->GetChisquare()/datHist->Integral());
	datHist->GetXaxis()->SetTitle(xTitle);
	datHist->GetYaxis()->SetTitle("Events");
	
	bkgdFunct = new TF1("bkgdFunct", bkgd, 60, 150, 4);
	bkgdFunct->SetParameter(0,fitFunct->GetParameter(0));
	bkgdFunct->SetParameter(1,fitFunct->GetParameter(1));
	bkgdFunct->SetParameter(2,fitFunct->GetParameter(2));
	bkgdFunct->SetParameter(3,fitFunct->GetParameter(3));
	//bkgdFunct->SetParameter(4,fitFunct->GetParameter(4));
	
	bkgdFunct->SetLineStyle(2);
	bkgdFunct->Draw("SAME");
	
	datEvts = datHist->Integral(datHist->FindBin(81),datHist->FindBin(103));
	bkgdEvts = integral(bkgdFunct,82,102); 
}

double fitStruct::getSignalEvt1()
{
	return fitFunct->GetParameter(4);
}

double fitStruct::getSignalEvt2()
{
	return datHist->Integral(datHist->FindBin(81),datHist->FindBin(103)) - integral(bkgdFunct,82,102);
}

double fitStruct::integral(TF1* theFunct, int low, int high)
{
	//NOTE: The current bin size for all my mass hists is 2 (GeV).
	int a, b;
	double sum = 0;
	
	a = ( low%2 == 0 )? (low+1) : low; 
	b = ( high%2 == 0 )? (high-1) : high;
	
	for(int i=a; i<=b; i+= 2)
	{
		sum += theFunct->Eval(i);
	}
	
	return sum;
}



int Efficiency04(double mean=0.01696, double sigma=0.0658)
{
	if(mean==0.01696 && sigma==0.0658)
	{
		std::cout<<"NOTE: Using All-NT mean and sigma for signal smearing!"<<std::endl;
	}
	
	int nvLow, nvHigh, nvStep, nvBins, ptBins;
	double ptLow, ptHigh, ptStep;
	
	nvBins = 12; //note: there are actually this +1 nvBins, with one for for overflow.	
	nvLow = 4;
	nvHigh = 28;
	nvStep = (int)((nvHigh-nvLow)/nvBins);
	
	ptBins = 14;
	ptLow = 20.0;
	ptHigh = 90.0;
	ptStep = (ptHigh-ptLow)/ptBins;
	
	gROOT->SetStyle("Plain");	
	gStyle->SetErrorX(0);
	gStyle->SetOptStat(0);
	
	TCanvas* c1 = new TCanvas("C1", "c1", 1200, 900);
	c1->cd();
	
	TLegend* l1 = new TLegend(0.65,0.7,0.9,0.9);
	l1->SetFillColor(10);
	TLegend* l2 = new TLegend(0.65,0.1,0.9,0.3);
	l2->SetFillColor(10);

	//get an ntuple to smear
	TFile* f1 = new TFile("/local/cms/user/finkel/NoTrack/Ntuple/DE_MC2012_AlexTrgr3/DE_MC2012_AlexTrgr_newerCuts.root");	
	if(f1 == NULL)
	{
		std::cout<<"Failed to open MC file. Exiting."<<std::endl;
		exit(1);
	}
	//and a data ntuple to fit
	TFile* f2 = new TFile("/local/cms/user/finkel/NoTrack/Ntuple/DE_ReReco_NewerCuts/DE_ReReco2012FULL_Jan13_AlexTrgr_newerCuts_001.root");	
	if(f2 == NULL)
	{
		std::cout<<"Failed to open Data file. Exiting."<<std::endl;
		exit(1);
	}
	
	//and make a file to store the results
	TFile* file = new TFile("NoTrackEfficiencies_01.root","recreate");
	
	//make relevant hists
	TH1D* unSmearedM = new TH1D("unSmearedM","M_{Z} before smearing, all-NT",75,0,150);
	TH1D* smearedM = new TH1D("smearedM","Signal, all-NT smearing",75,0,150);	
	TH1D* datMbefore = new TH1D("datMbefore","All MZ from Data, no NT cut",75,0,150);
	TH1D* datMafter = new TH1D("datMafter","All MZ from Data, with NT cut",75,0,150);
	
	//diagnostic high-pt hitmaps
	std::vector<TH2D*> highPtHitmap;
	highPtHitmap.push_back( new TH2D("highPtHitmap0", "#eta-#phi for Pt>60 GeV NT-", 20, -3.14, 3.14 , 20, 2.5, 3.0) );
	highPtHitmap.push_back( new TH2D("highPtHitmap1", "#eta-#phi for Pt>60 GeV NT+", 20, -3.14, 3.14 , 20, 2.5, 3.0) );
	highPtHitmap.push_back( new TH2D("highPtHitmap2", "x-y for Pt>60 GeV NT-", 20, -0.6, 0.6 , 20, -0.6, 0.6) );
	highPtHitmap.push_back( new TH2D("highPtHitmap3", "x-y for Pt>60 GeV NT+", 20, -0.6, 0.6 , 20, -0.6, 0.6) );
	highPtHitmap.push_back( new TH2D("highPtHitmap4", "ix-iy for Pt>60 GeV NT-", 40, 30, 70 , 40, 30, 70) );
	highPtHitmap.push_back( new TH2D("highPtHitmap5", "ix-iy for Pt>60 GeV NT+", 40, 30, 70 , 40, 30, 70) );
	
	std::vector<TH1D*> datBeforeByNV, datAfterByNV, datBeforeByPt, datAfterByPt, unSmearedMbyPt, smearedMbyPt;
	
	datBeforeByNV.push_back(datMbefore);
	datAfterByNV.push_back(datMafter);
	datBeforeByPt.push_back(datMbefore);
	datAfterByPt.push_back(datMafter);
	unSmearedMbyPt.push_back(unSmearedM);
	smearedMbyPt.push_back(smearedM);
	
	for(int i=1; i<=nvBins; i++)
	{
		char nameNVbefore[100], titleNVbefore[200], nameNVafter[100], titleNVafter[200];	
		
		sprintf(nameNVbefore,"MbeforeNV%d",i);
		sprintf(nameNVafter,"MafterNV%d",i);
		sprintf(titleNVbefore,"MZ No NT cut, %d<nVert<%d",nvLow+nvStep*(i-1),nvLow+nvStep*i);
		sprintf(titleNVafter,"MZ With NT cut, %d<nVert<%d",nvLow+nvStep*(i-1),nvLow+nvStep*i);

		datBeforeByNV.push_back( new TH1D(nameNVbefore,titleNVbefore,75,0,150) );
		datAfterByNV.push_back( new TH1D(nameNVafter,titleNVafter,75,0,150) );
		//if this works, add also the vs.Pt hists!
	}
	
	for (int j=1; j<=ptBins; j++)
	{
		char namePtBefore[100], titlePtBefore[200], namePtAfter[100], titlePtAfter[200], unSmearedName[200], unSmearedTitle[200], smearedName[200], smearedTitle[200];
		
		sprintf(namePtBefore,"MbeforePt%d",j);
		sprintf(namePtAfter,"MafterPt%d",j);
		sprintf(unSmearedName,"unSmearedMbyPt%d",j);
		sprintf(smearedName,"smearedMbyPt%d",j);
		sprintf(titlePtBefore,"MZ No NT cut, %d<Pt<%d",(int)(ptLow+ptStep*(j-1)),(int)(ptLow+ptStep*j));
		sprintf(titlePtAfter,"MZ With NT cut, %d<Pt<%d",(int)(ptLow+ptStep*(j-1)),(int)(ptLow+ptStep*j));
		sprintf(unSmearedTitle,"Un-smeared MZ, %d<Pt<%d",(int)(ptLow+ptStep*(j-1)),(int)(ptLow+ptStep*j));
		sprintf(smearedTitle,"Smeared MZ, %d<Pt<%d",(int)(ptLow+ptStep*(j-1)),(int)(ptLow+ptStep*j));
		
		datBeforeByPt.push_back( new TH1D(namePtBefore,titlePtBefore,75,0,150) );
		datAfterByPt.push_back( new TH1D(namePtAfter,titlePtAfter,75,0,150) );
		unSmearedMbyPt.push_back( new TH1D(unSmearedName,unSmearedTitle,75,0,150) );
		smearedMbyPt.push_back( new TH1D(smearedName,smearedTitle,75,0,150) );
	}
	//do overflow nvBins separately!
	char nameNVbefore[100], titleNVbefore[200], nameNVafter[100], titleNVafter[200];	
	
	sprintf(nameNVbefore,"MbeforeNV%d",nvBins+1);
	sprintf(nameNVafter,"MafterNV%d",nvBins+1);
	sprintf(titleNVbefore,"MZ No NT cut, %d<nVert",nvHigh);
	sprintf(titleNVafter,"MZ With NT cut, %d<nVert",nvHigh);
	
	datBeforeByNV.push_back( new TH1D(nameNVbefore,titleNVbefore,75,0,150) );
	datAfterByNV.push_back( new TH1D(nameNVafter,titleNVafter,75,0,150) );
	
	//same for Pt
	char namePtBefore[100], titlePtBefore[200], namePtAfter[100], titlePtAfter[200], unSmearedName[200], unSmearedTitle[200], smearedName[200], smearedTitle[200];	
	
	sprintf(namePtBefore,"MbeforePt%d",ptBins+1);
	sprintf(namePtAfter,"MafterPt%d",ptBins+1);
	sprintf(unSmearedName,"UnSmearedMPt%d",ptBins+1);
	sprintf(smearedName,"SmearedMPt%d",ptBins+1);
	sprintf(titlePtBefore,"MZ No NT cut, %d<Pt",(int)ptHigh);
	sprintf(titlePtAfter,"MZ With NT cut, %d<Pt",(int)ptHigh);
	sprintf(unSmearedTitle,"Un-Smeared MZ, %d<Pt",(int)ptHigh);
	sprintf(smearedTitle,"Smeared MZ, %d<Pt",(int)ptHigh);
	
	datBeforeByPt.push_back( new TH1D(namePtBefore,titlePtBefore,75,0,150) );
	datAfterByPt.push_back( new TH1D(namePtAfter,titlePtAfter,75,0,150) );
	unSmearedMbyPt.push_back( new TH1D(unSmearedName,unSmearedTitle,75,0,150) );
	smearedMbyPt.push_back( new TH1D(smearedName,smearedTitle,75,0,150) );
	
	//shower-shape diagrams
	TH1D* eIsoBefore = new TH1D("eIsoBefore","Ecal Iso",50,0,0.1);
	TH1D* hIsoBefore = new TH1D("hIsoBefore","Hcal Iso",50,0,0.25);
	TH1D* HoEMBefore = new TH1D("HoEMBefore","Hardonic/EM",50,0,0.1);
	TH1D* rNineBefore = new TH1D("rNineBefore","R9",50,0.8,1.05);
	TH1D* sIeIeBefore = new TH1D("sIeIeBefore","sIeIe",50,0,0.05);
	
	TH1D* eIsoAfter = new TH1D("eIsoAfter","Ecal Iso",50,0,0.1);
	TH1D* hIsoAfter = new TH1D("hIsoAfter","Hcal Iso",50,0,0.25);
	TH1D* HoEMAfter = new TH1D("HoEMAfter","Hardonic/EM",50,0,0.1);
	TH1D* rNineAfter = new TH1D("rNineAfter","R9",50,0.8,1.05);
	TH1D* sIeIeAfter = new TH1D("sIeIeAfter","sIeIe",50,0,0.05);
	
	std::cout<<"Ready to fill hists!"<<std::endl;
	
	
	//open file for the Data ntuple
	ZEffTree *ze2, *ze1;
	ze1 = new ZEffTree(*f1,false);
	ze2 = new ZEffTree(*f2,false);
	const std::string NT = "NTLooseElectronId-EtaDet";
	
	//for diagnostic hitmaps:
	double x, y, theta, Phi; //not to be confused with "phi" used for TLorentzVector!
	double z = 3.154;
	
	
	TH1D* test = new TH1D("test","test", 10, 50, 150);//diagnostic hist
	
	//Fill the data histograms:
	bool run = true;
	while(run) //fill the data hists
	{
			ze2->Entries();
			test->Fill(ze2->reco.mz);

		if( (fabs(ze2->reco.eta[1])>2.5) && (fabs(ze2->reco.eta[1])<3.0) ) 
		{	
			//fill the "before" hists
			datBeforeByNV[0]->Fill(ze2->reco.mz);
			//if(ze2->reco.nverts < nvLow) continue;
			if(ze2->reco.nverts > nvHigh)
			{
				datBeforeByNV[nvBins+1]->Fill(ze2->reco.mz);
				//continue;
			}
			for(int i=1; i<=nvBins; i++)
			{
				
				if( ( ze2->reco.nverts > nvLow+nvStep*(i-1) ) && ( ze2->reco.nverts < nvLow+nvStep*i ) )
				{
					datBeforeByNV[i]->Fill(ze2->reco.mz);
					break;
				}
			}
			
			datBeforeByPt[0]->Fill(ze2->reco.mz);
			for(int j=1; j<=ptBins; j++)
			{
				if(ze2->reco.pt[1] < ptLow) break;
				if(ze2->reco.pt[1] > ptHigh)
				{
					datBeforeByPt[ptBins+1]->Fill(ze2->reco.mz);
					
					theta = 2*atan( exp( ze2->reco.eta[1]) );
					Phi = ze2->reco.phi[1];
					x = z*tan(theta)*cos(Phi);
					y = z*tan(theta)*sin(Phi);
					if(ze2->reco.eta[1] < 0) 
					{
						highPtHitmap[0]->Fill(ze2->reco.phi[1], -1*ze2->reco.eta[1]);
						highPtHitmap[2]->Fill(x, y);
						highPtHitmap[4]->Fill(ze2->reco.ix[1], ze2->reco.iy[1]);
					}
					else 
					{
						highPtHitmap[1]->Fill(ze2->reco.phi[1], ze2->reco.eta[1]);
						highPtHitmap[3]->Fill(x, y);
						highPtHitmap[5]->Fill(ze2->reco.ix[1], ze2->reco.iy[1]);
					}
					break;
				}
				if( ( ze2->reco.pt[1] > ptLow+ptStep*(j-1) ) && ( ze2->reco.pt[1] < ptLow+ptStep*j ) )
				{
					datBeforeByPt[j]->Fill(ze2->reco.mz);
					break;
				}
			}
			
			//fill the shower-shape hists
			eIsoBefore->Fill(ze2->reco.Eiso[1]);
			hIsoBefore->Fill(ze2->reco.Hiso[1]);
			rNineBefore->Fill(ze2->reco.RNine[1]);
			HoEMBefore->Fill(ze2->reco.HoEM[1]);
			sIeIeBefore->Fill(ze2->reco.Sieie[1]);
				
			//now fill the "after" ones
			if( ze2->reco.isSelected(1,NT) )
			{
				datAfterByNV[0]->Fill(ze2->reco.mz);
				//if(ze2->reco.nverts < nvLow) continue;
				if(ze2->reco.nverts > nvHigh)
				{
					datAfterByNV[nvBins+1]->Fill(ze2->reco.mz);
					//break;
				}
				for(int i=1; i<=nvBins; i++)
				{
					if( ( ze2->reco.nverts > nvLow+nvStep*(i-1) ) && ( ze2->reco.nverts < nvLow+nvStep*i ) )
					{
						datAfterByNV[i]->Fill(ze2->reco.mz);
						break;
					}
				}
				
				datAfterByNV[0]->Fill(ze2->reco.mz);
				for(int j=1; j<=ptBins; j++)
				{
					if(ze2->reco.pt[1] < ptLow) break;
					if(ze2->reco.pt[1] > ptHigh)
					{
						datAfterByPt[ptBins+1]->Fill(ze2->reco.mz);
						break;
					}
					if( ( ze2->reco.pt[1] > ptLow+ptStep*(j-1) ) && ( ze2->reco.pt[1] < ptLow+ptStep*j ) )
					{
						datAfterByPt[j]->Fill(ze2->reco.mz);
						break;
					}
				}
				
				eIsoAfter->Fill(ze2->reco.Eiso[1]);
				hIsoAfter->Fill(ze2->reco.Hiso[1]);
				rNineAfter->Fill(ze2->reco.RNine[1]);
				HoEMAfter->Fill(ze2->reco.HoEM[1]);
				sIeIeAfter->Fill(ze2->reco.Sieie[1]);
			}
		}
		
		run = ze2->GetNextEvent();
	}
	
	test->Draw();
	c1->Print("Efficiencies/test.png");
	c1->Clear();
	
	TLorentzVector elec1, elec2, theZ;
	float pt, eta, phi, E;
	int randSeed=123456; //using constant seed for reproducibility
	TRandom3* rand = new TRandom3(randSeed);
	
	int nEntries = 0;
	run = true;
	while (run) //fill MC hist and do smearing on the events
	{    	
		ze1->Entries();
		nEntries++;
		if( (fabs(ze1->reco.eta[1])>2.5) && (fabs(ze1->reco.eta[1])<3.0) && (ze1->reco.isSelected(1,NT)) 
			)//range not limited for MC, since it does not know about crystal badness
		{
			pt = ze1->reco.pt[0];    	
			eta = ze1->reco.eta[0];
			phi = ze1->reco.phi[0];
			E = pt*cosh(eta);
			elec1.SetPtEtaPhiE(pt,eta,phi,E);
			
			pt = ze1->reco.pt[1]; //believe it or not, HERE is the smearing part:
			pt *= ( 1 + rand->Gaus(mean,sigma) );//note: using PROPORTIONAL correction now    	
			eta = ze1->reco.eta[1];
			phi = ze1->reco.phi[1];    
			E = pt*cosh(eta);	
			elec2.SetPtEtaPhiE(pt,eta,phi,E);
			
			theZ = elec1+elec2;
			
			unSmearedMbyPt[0]->Fill(ze1->reco.mz);
			smearedMbyPt[0]->Fill(theZ.M());
			for(int j=1; j<=ptBins; j++)
			{
				if(pt < ptLow) break;
				if(pt > ptHigh)
				{
					unSmearedMbyPt[ptBins+1]->Fill(ze1->reco.mz);
					smearedMbyPt[ptBins+1]->Fill(theZ.M());
					break;
				}
				//std::cout<<"Ping!"<<std::endl;
				if( ( pt > ptLow+ptStep*(j-1) ) && ( pt < ptLow+ptStep*j ) )
				{
					unSmearedMbyPt[j]->Fill(ze1->reco.mz);
					smearedMbyPt[j]->Fill(theZ.M());
					break;
				}
			}
		}
		
		run = ze1->GetNextEvent();
	}
	
	std::cout<<"Hists made, starting the fitting process."<<std::endl;
	
	//arrays to hold efficiency values and the points they correspond to:
	double effVsNV1[nvBins+2], effVsPt1[ptBins+2], effVsNV2[nvBins+2], effVsPt2[ptBins+2];
	double ptBinArr[ptBins+1], nvBinArr[nvBins+1];
	double effVsNVComb[nvBins+1], effVsPtComb[ptBins+1];
	double effVsNVUpStErr[nvBins+1], effVsPtUpStErr[ptBins+1], effVsNVDownStErr[nvBins+1], effVsPtDownStErr[ptBins+1];
	double effVsNVUpErr[nvBins+1], effVsPtUpErr[ptBins+1], effVsNVDownErr[nvBins+1], effVsPtDownErr[ptBins+1];
	
	//make a hist to show efficiency vs. nvert and Pt by fitting the peak
	TH1D* effVsNVert1 = new TH1D("effVsNVert1","NT Cut Efficiency vs. nVert By Fit",nvBins+1,nvLow,nvHigh+nvStep);
	TH1D* effVsPtrans1 = new TH1D("effVsPtrans1","NT Cut Efficiency vs. Pt By Fit",ptBins+1,ptLow,ptHigh+ptStep);
	
	//and by background subtraction method
	TH1D* effVsNVert2 = new TH1D("effVsNVert2","NT Cut Efficiency vs. nVert By Bkgd Subtraction",nvBins+1,nvLow,nvHigh+nvStep);
	TH1D* effVsPtrans2 = new TH1D("effVsPtrans2","NT Cut Efficiency vs. Pt By Bkgd Subtraction",ptBins+1,ptLow,ptHigh+ptStep);
	
	//hists for BayesDivide; they will be filled with the AVERAGE signal values from the two methods
	TH1D* signalBeforeByNv = new TH1D("signalBeforeByNv","",nvBins+1,nvLow,nvHigh+nvStep);
	TH1D* signalBeforeByPt = new TH1D("signalBeforeByPt","",ptBins+1,ptLow,ptHigh+ptStep);
	TH1D* signalAfterByNv = new TH1D("signalAfterByNv","",nvBins+1,nvLow,nvHigh+nvStep);
	TH1D* signalAfterByPt = new TH1D("signalAfterByPt","",ptBins+1,ptLow,ptHigh+ptStep);

	double sigEvtBefore1, sigEvtAfter1, sigEvtBefore2, sigEvtAfter2;
	//do the NV efficiencies
	for(int i=0; i<=nvBins+1; i++)
	{
		char pngNameBefore[128];
		char pngNameAfter[128];
		sprintf(pngNameBefore, "Efficiencies/EffByNvertex/EffNvertBeforeBin%d.png",i);
		sprintf(pngNameAfter, "Efficiencies/EffByNvertex/EffNvertAfterBin%d.png",i);
		
		fitStruct before(smearedMbyPt[0],datBeforeByNV[i]);
		before.doFit();
		sigEvtBefore1 = before.getSignalEvt1();
		sigEvtBefore2 = before.getSignalEvt2();
		c1->Print(pngNameBefore);
		c1->Clear();
		
		fitStruct after(smearedMbyPt[0],datAfterByNV[i]);
		after.doFit();
		sigEvtAfter1 = after.getSignalEvt1();
		sigEvtAfter2 = after.getSignalEvt2();
		c1->Print(pngNameAfter);
		c1->Clear();
		
		effVsNV1[i] = sigEvtAfter1/sigEvtBefore1;		
		effVsNVert1->SetBinContent(i,effVsNV1[i]);
		
		effVsNV2[i] = sigEvtAfter2/sigEvtBefore2;		
		effVsNVert2->SetBinContent(i,effVsNV2[i]);
		
		signalBeforeByNv->SetBinContent(i, 0.5*(sigEvtBefore1+sigEvtBefore2) );
		signalAfterByNv->SetBinContent(i, 0.5*(sigEvtAfter1+sigEvtAfter2) );
		
		//Efficiency Stats struct:
		EfficiencyStatistics esNv( effVsNV1[i], (int)sigEvtBefore1 );
		//effVsNVert1->SetBinError(i,esNv.sigma(1)-effVsNV1[i]);
		
		if(i>0)//necessary to avoid the "everything" entries
		{
			nvBinArr[i-1] = nvLow+nvStep*(i-0.5);
			effVsNVComb[i-1] = 0.5*(effVsNV1[i]+effVsNV2[i]);
			//stat errors
			effVsNVUpStErr[i-1] = esNv.sigma(1) - effVsNV1[i];
			effVsNVDownStErr[i-1] = effVsNV1[i]-esNv.sigma(-1);
			//full errors, add in quadrature!
			effVsNVUpErr[i-1] = sqrt( pow(esNv.sigma(1) - effVsNV1[i],2) + pow(0.5*(effVsNV1[i]-effVsNV2[i]),2) );
			effVsNVDownErr[i-1] = sqrt( pow(effVsNV1[i]-esNv.sigma(-1),2) + pow(0.5*(effVsNV1[i]-effVsNV2[i]),2) );
		}
	}
	
	//now do the pt ones
	for(int j=0; j<=ptBins+1; j++)
	{
		char pngNameBefore[128];
		char pngNameAfter[128];
		char unSmearedPng[128];
		sprintf(pngNameBefore, "Efficiencies/EffByPt/EffPtBeforeBin%d.png",j);
		sprintf(pngNameAfter, "Efficiencies/EffByPt/EffPtAfterBin%d.png",j);
		sprintf(unSmearedPng, "Efficiencies/UnsmearedM/UnSmearedMbyPt%d.png",j);
		
		fitStruct before(smearedMbyPt[j],datBeforeByPt[j]);
		before.doFit();
		sigEvtBefore1 = before.getSignalEvt1();
		sigEvtBefore2 = before.getSignalEvt2();
		c1->Print(pngNameBefore);
		c1->Clear();
		
		fitStruct after(smearedMbyPt[j],datAfterByPt[j]);
		after.doFit();
		sigEvtAfter1 = after.getSignalEvt1();
		sigEvtAfter2 = after.getSignalEvt2();
		c1->Print(pngNameAfter);
		c1->Clear();
		
		effVsPt1[j] = sigEvtAfter1/sigEvtBefore1;
		effVsPtrans1->SetBinContent(j,effVsPt1[j]);
		
		effVsPt2[j] = sigEvtAfter2/sigEvtBefore2;
		effVsPtrans2->SetBinContent(j,effVsPt2[j]);
		
		signalBeforeByPt->SetBinContent(j, 0.5*(sigEvtBefore1+sigEvtBefore2) );
		signalAfterByPt->SetBinContent(j, 0.5*(sigEvtAfter1+sigEvtAfter2) );
		
		//Efficiency Stats struct:
		EfficiencyStatistics esPt( effVsPt1[j], (int)sigEvtBefore1 );
		//effVsPtrans1->SetBinError(j,esPt.sigma(1)-effVsPt1[j]);
		
		if(j>0)//necessary ot avoid the "everything" entries
		{
			ptBinArr[j-1] = ptLow+ptStep*(j-0.5);
			effVsPtComb[j-1] = 0.5*(effVsPt1[j]+effVsPt2[j]);
			//stat errors
			effVsPtUpStErr[j-1] = esPt.sigma(1)-effVsPt1[j];
			effVsPtDownStErr[j-1] = effVsPt1[j]-esPt.sigma(-1);
			//full errors
			effVsPtUpErr[j-1] = sqrt( pow(esPt.sigma(1)-effVsPt1[j],2) + pow(0.5*(effVsPt1[j]-effVsPt2[j]),2) );
			effVsPtDownErr[j-1] = sqrt( pow(effVsPt1[j]-esPt.sigma(-1),2) + pow(0.5*(effVsPt1[j]-effVsPt2[j]),2) );
		}
		
		//draw diagnostic unsmeared M hists:
		
		unSmearedMbyPt[j]->SetLineColor(2);
		unSmearedMbyPt[j]->SetLineWidth(2);
		unSmearedMbyPt[j]->GetXaxis()->SetRangeUser(60,130);
		unSmearedMbyPt[j]->GetXaxis()->SetTitle("M_{ee} (GeV)");
		unSmearedMbyPt[j]->GetYaxis()->SetTitle("Events");
		unSmearedMbyPt[j]->Draw("HIST");
		c1->Print(unSmearedPng);
		c1->Clear();
	}
	
	//shower-shape plots now done in separate script
	/*eIsoBefore->SetTitle(";EcalIso;Events");
	eIsoBefore->SetLineWidth(3);
	eIsoBefore->SetLineColor(4);
	eIsoBefore->Draw("hist");
	eIsoAfter->SetTitle("");
	eIsoAfter->SetLineWidth(3);
	eIsoAfter->SetLineColor(2);
	eIsoAfter->Draw("same hist");
	c1->Print("Efficiencies/EcalIso.png");
	c1->Clear();
	
	hIsoBefore->SetTitle(";HcalIso;Events");
	hIsoBefore->SetLineWidth(3);
	hIsoBefore->SetLineColor(4);
	hIsoBefore->Draw("hist");
	hIsoAfter->SetTitle("");
	hIsoAfter->SetLineWidth(3);
	hIsoAfter->SetLineColor(2);
	hIsoAfter->Draw("same hist");
	c1->Print("Efficiencies/HcalIso.png");
	c1->Clear();
	
	HoEMBefore->SetTitle(";E_{H}/E_{EM};Events");
	HoEMBefore->SetLineWidth(3);
	HoEMBefore->SetLineColor(4);
	HoEMBefore->Draw("hist");
	HoEMAfter->SetTitle("");
	HoEMAfter->SetLineWidth(3);
	HoEMAfter->SetLineColor(2);
	HoEMAfter->Draw("same hist");
	c1->Print("Efficiencies/HadrOverEM.png");
	c1->Clear();
	
	rNineBefore->SetTitle(";R9;Events");
	rNineBefore->SetLineWidth(3);
	rNineBefore->SetLineColor(4);
	rNineBefore->Draw("hist");
	rNineAfter->SetTitle("");
	rNineAfter->SetLineWidth(3);
	rNineAfter->SetLineColor(2);
	rNineAfter->Draw("same hist");
	c1->Print("Efficiencies/R9.png");
	c1->Clear();
	
	sIeIeBefore->SetTitle(";#sigma_{i#eta,i#eta};Events");
	sIeIeBefore->SetLineWidth(3);
	sIeIeBefore->SetLineColor(4);
	sIeIeBefore->Draw("hist");
	sIeIeAfter->SetTitle("");
	sIeIeAfter->SetLineWidth(3);
	sIeIeAfter->SetLineColor(2);
	sIeIeAfter->Draw("same hist");
	c1->Print("Efficiencies/Sieie.png");
	c1->Clear();*/
	
	//Efficiency plots
	effVsNVert1->SetMarkerStyle(20);
	effVsNVert1->Draw("P");
	c1->Print("Efficiencies/EffPlots/EfficiencyVsNvertFit.png");
	
	effVsPtrans1->SetMarkerStyle(20);
	effVsPtrans1->Draw("P");
	c1->Print("Efficiencies/EffPlots/EfficiencyVsPtFit.png");
	
	effVsNVert2->SetMarkerStyle(20);
	effVsNVert2->Draw("P");
	c1->Print("Efficiencies/EffPlots/EfficiencyVsNvertBkgdSub.png");
	
	effVsPtrans2->SetMarkerStyle(20);
	effVsPtrans2->Draw("P");
	c1->Print("Efficiencies/EffPlots/EfficiencyVsPtBkgdSub.png");
	
	for(int k=0; k<(int)highPtHitmap.size(); k++)
	{
		char pngName[128];
		sprintf(pngName, "Efficiencies/Hitmaps/highPtHitmap%d.png", k);
		highPtHitmap[k]->Draw("BOX");
		c1->Print(pngName);
		c1->Clear();
	}
	
	//comparison plots
	effVsNVert2->SetTitle(";N(vertices);Efficiency");
	effVsNVert2->SetMarkerStyle(20);
	effVsNVert2->SetMarkerColor(2);
	effVsNVert2->GetYaxis()->SetRangeUser(0.25,1);
	effVsNVert2->Draw("P");
	effVsNVert1->SetMarkerStyle(21);
	effVsNVert1->SetMarkerColor(4);
	effVsNVert1->SetLineColor(4);
	effVsNVert1->Draw("SAME P");
	l1->AddEntry(effVsNVert2,"Bkgd Subtrac");
	l1->AddEntry(effVsNVert1,"Peak Fit");
	l1->Draw();
	c1->Print("Efficiencies/EffPlots/EfficiencyVsNvertComparison.png");
	//l1->Clear();
	
	effVsPtrans2->SetTitle(";e p_{T} (GeV);Efficiency");
	effVsPtrans2->SetMarkerStyle(20);
	effVsPtrans2->SetMarkerColor(2);
	effVsPtrans2->GetYaxis()->SetRangeUser(0.25,1);
	effVsPtrans2->Draw("P");
	effVsPtrans1->SetMarkerStyle(21);
	effVsPtrans1->SetMarkerColor(4);
	effVsPtrans1->SetLineColor(4);
	effVsPtrans1->Draw("SAME P");
	l2->AddEntry(effVsPtrans2,"Bkgd Subtrac");
	l2->AddEntry(effVsPtrans1,"Peak Fit");
	l2->Draw();
	c1->Print("Efficiencies/EffPlots/EfficiencyVsPtComparison.png");
	c1->Clear();
	
	//combination plots
	TGraphAsymmErrors* EffVsNvCombined1 = new TGraphAsymmErrors(nvBins+1, nvBinArr, effVsNVComb, 0,0,effVsNVDownErr, effVsNVUpErr);
	EffVsNvCombined1->SetTitle(";N(vertices);Efficiency");
	EffVsNvCombined1->GetYaxis()->SetRangeUser(0.25,1.05);
	EffVsNvCombined1->SetMarkerStyle(20);
	EffVsNvCombined1->SetMarkerColor(kBlue+3);
	EffVsNvCombined1->SetLineColor(4);
	EffVsNvCombined1->Draw("AP");
	TGraphAsymmErrors* EffVsNvCombined2 = new TGraphAsymmErrors(nvBins+1, nvBinArr, effVsNVComb, 0,0,effVsNVDownStErr, effVsNVUpStErr);	
	EffVsNvCombined2->SetTitle(";;");
	EffVsNvCombined2->SetMarkerStyle(20);
	EffVsNvCombined2->SetMarkerColor(kBlue+3);
	EffVsNvCombined2->SetLineColor(4);
	EffVsNvCombined2->Draw("SAME P");
	c1->Update();
	c1->Print("Efficiencies/EffPlots/EfficiencyVsNvertCombined.png");
	
	TGraphAsymmErrors* EffVsPtCombined1 = new TGraphAsymmErrors(ptBins+1, ptBinArr, effVsPtComb, 0,0,effVsPtDownErr, effVsPtUpErr);
	EffVsPtCombined1->SetTitle(";e p_{T} (GeV);Efficiency");
	EffVsPtCombined1->GetYaxis()->SetRangeUser(0.25,1.05);
	EffVsPtCombined1->SetMarkerStyle(20);
	EffVsPtCombined1->SetMarkerColor(kBlue+3);
	EffVsPtCombined1->SetLineColor(4);
	EffVsPtCombined1->Draw("AP");
	TGraphAsymmErrors* EffVsPtCombined2 = new TGraphAsymmErrors(ptBins+1, ptBinArr, effVsPtComb, 0,0,effVsPtDownStErr, effVsPtUpStErr);
	EffVsPtCombined2->SetTitle(";;");
	EffVsPtCombined2->SetMarkerStyle(20);
	EffVsPtCombined2->SetMarkerColor(kBlue+3);
	EffVsPtCombined2->SetLineColor(4);
	EffVsPtCombined2->Draw("SAME P");
	c1->Update();
	c1->Print("Efficiencies/EffPlots/EfficiencyVsPtCombined.png");
	
	TGraphAsymmErrors *EffVsNvWithRoot = new TGraphAsymmErrors();
	TGraphAsymmErrors *EffVsPtWithRoot = new TGraphAsymmErrors();
	
	//stuff for making double error-bars using Root's stat uncertainties		
	double NvEff[nvBins+1], NvErrUp[nvBins+1], NvErrDown[nvBins+1];
	double pointX(0), pointY(0);
	double PtEff[ptBins+1], PtErrUp[ptBins+1], PtErrDown[ptBins+1];
	
	EffVsNvWithRoot->BayesDivide(signalAfterByNv,signalBeforeByNv);
	for(int i=0; i<=nvBins+1; i++)
	{
		EffVsNvWithRoot->SetPointEXlow(i,0);
		EffVsNvWithRoot->SetPointEXhigh(i,0);
		NvErrUp[i] = sqrt( pow(0.5*(effVsNV1[i+1]-effVsNV2[i+1]),2) + pow(EffVsNvWithRoot->GetErrorYhigh(i) ,2) );
		NvErrDown[i] = sqrt( pow(0.5*(effVsNV1[i+1]-effVsNV2[i+1]),2) + pow(EffVsNvWithRoot->GetErrorYlow(i) ,2) );
		EffVsNvWithRoot->GetPoint(i,pointX,pointY);
		NvEff[i] = pointY;
	}
	EffVsNvWithRoot->SetTitle(";N(vertices);Efficiency");
	EffVsNvWithRoot->GetYaxis()->SetRangeUser(0.25,1.05);
	EffVsNvWithRoot->GetXaxis()->SetRangeUser(nvLow,nvHigh+nvStep);
	EffVsNvWithRoot->SetMarkerStyle(20);
	EffVsNvWithRoot->SetMarkerColor(2);
	EffVsNvWithRoot->SetLineColor(2);
	EffVsNvWithRoot->Draw("AP");
	TGraphAsymmErrors *EffVsNvWithRoot2 = new TGraphAsymmErrors(nvBins+1, nvBinArr, NvEff, 0,0,NvErrDown, NvErrUp);
	EffVsNvWithRoot2->SetMarkerStyle(20);
	EffVsNvWithRoot2->SetMarkerColor(4);
	EffVsNvWithRoot2->SetLineColor(4);
	EffVsNvWithRoot2->Draw("SAME P");
	c1->Update();
	c1->Print("Efficiencies/EffPlots/EfficiencyVsNvertRoot.png");
	
	EffVsPtWithRoot->BayesDivide(signalAfterByPt,signalBeforeByPt);
	for(int j=0; j<=ptBins+1; j++)
	{
		EffVsPtWithRoot->SetPointEXlow(j,0);
		EffVsPtWithRoot->SetPointEXhigh(j,0);
		PtErrUp[j] = sqrt( pow(0.5*(effVsPt1[j+1]-effVsPt2[j+1]),2) + pow(EffVsPtWithRoot->GetErrorYhigh(j) ,2) );
		PtErrDown[j] = sqrt( pow(0.5*(effVsPt1[j+1]-effVsPt2[j+1]),2) + pow(EffVsPtWithRoot->GetErrorYlow(j) ,2) );
		EffVsPtWithRoot->GetPoint(j,pointX,pointY);
		PtEff[j] = pointY;
	}
	EffVsPtWithRoot->SetTitle(";e p_{T} (GeV);Efficiency");
	EffVsPtWithRoot->GetYaxis()->SetRangeUser(0.25,1.05);
	EffVsPtWithRoot->GetXaxis()->SetRangeUser(ptLow,ptHigh+ptStep);
	EffVsPtWithRoot->SetMarkerStyle(20);
	EffVsPtWithRoot->SetMarkerColor(2);
	EffVsPtWithRoot->SetLineColor(2);
	EffVsPtWithRoot->Draw("AP");
	TGraphAsymmErrors *EffVsPtWithRoot2 = new TGraphAsymmErrors(ptBins+1, ptBinArr, PtEff, 0,0,PtErrDown, PtErrUp);
	EffVsPtWithRoot2->SetMarkerStyle(20);
	EffVsPtWithRoot2->SetMarkerColor(4);
	EffVsPtWithRoot2->SetLineColor(4);
	EffVsPtWithRoot2->Draw("SAME P");
	c1->Update();
	c1->Print("Efficiencies/EffPlots/EfficiencyVsPtRoot.png");
	
	
	
	c1->Close();
	
	file->Write();
	
	std::cout<<"There are "<<rNineBefore->Integral()<<" entries in R9 Before cuta and "<<rNineAfter->Integral()<<" after."<<std::endl;
	
	return 0;
}

























