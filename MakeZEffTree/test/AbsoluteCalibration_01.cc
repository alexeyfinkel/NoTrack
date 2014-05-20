//this is a copy of the MC projected cal consts script, only ran on data.

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

int calibrateData()
{
    gROOT->SetStyle("Plain");	
	gStyle->SetErrorX(0);
	gStyle->SetOptStat(0);
	
	TCanvas* c1 = new TCanvas("C1", "c1", 1000, 800);
	c1->cd();
	//c1->SetLeftMargin(0.08);
	//c1->SetRightMargin(0.12);
    TLegend* l1 = new TLegend(0.6,0.6,0.9,0.9);
	l1->SetFillColor(10);

	//iteration number (should be about 10-15)
	const int nIterations = 15;
    
    //let's be fancy and use a map!
	std::map< std::pair<int,int>, TH1D* > ratioMapP, ratioMapN, constMapN, constMapP;
    char title[256],name[128];//, filename[128];
    
    for(int i=25;i<76;i++)
	{
		for(int j=25;j<76;j++)
		{
            //if( (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5))  < 132.25) || (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5)) > 342.25) ) continue;
            sprintf(name,"RatioIx%diy%dzP",i,j);
            sprintf(title,"E_{exp}/E_{obs}, NT+ i_{x} = %d, i_{y} = %d;E_{exp}/E_{obs};Events",i,j);
            ratioMapP[std::make_pair(i,j)]= new TH1D(name,title,24,0.4,1.6);
            sprintf(name,"RatioIx%diy%dzN",i,j);
            sprintf(title,"E_{exp}/E_{obs}, NT- i_{x} = %d, i_{y} = %d;E_{exp}/E_{obs};Events",i,j);			
            ratioMapN[std::make_pair(i,j)]= new TH1D(name,title,24,0.4,1.6);
            sprintf(name,"CalConstIx%diy%dzN",i,j);
            sprintf(title,"Calib. Const, NT- i_{x} = %d, i_{y} = %d;Iteration;CalConst",i,j);	
            constMapN[std::make_pair(i,j)]= new TH1D(name,title,nIterations,0,nIterations);  
            sprintf(name,"CalConstIx%diy%dzP",i,j);
            sprintf(title,"Calib. Const, NT+ i_{x} = %d, i_{y} = %d;Iteration;CalConst",i,j);	
            constMapP[std::make_pair(i,j)]= new TH1D(name,title,nIterations,0,nIterations);	

            //initialize consts	to zero, just in case (not actually used):
            constMapN[std::make_pair(i,j)]->SetBinContent(0,1);
            constMapP[std::make_pair(i,j)]->SetBinContent(0,1);
		}
	}
    //cal.consts. hists:
    TH2D *decalConstsN = new TH2D("decalConstsN","Decalibration Constants, NT-;i_{x};i_{y}",40,30,70,40,30,70);
	TH2D *decalConstsP = new TH2D("decalConstsp","Decalibration Constants, NT+;i_{x};i_{y}",40,30,70,40,30,70);
	TH2D *absConstsN = new TH2D("absConstsN","Absolute Calibration Constants, NT-;i_{x};i_{y}",40,30,70,40,30,70);
	TH2D *absConstsP = new TH2D("absConstsP","Absolute Calibration Constants, NT+;i_{x};i_{y}",40,30,70,40,30,70);
	TH1D *absConsts1D = new TH1D("absConsts1D","Absolute Calibration Constant Distribution",50,0,2.5);
	
	//Z mass hists:
	TH1D *mZbefore = new TH1D("mZbefore","Z Mass with Old Calibration;M_{ee} (GeV)",30,60,120);
	TH1D *mZoriginal = new TH1D("mZoriginal","Decalibrated Z Mass;M_{ee} (GeV)",30,60,120);
	TH1D *mZafter = new TH1D("mZafter","Z Mass Ater Absolute Calibration;M_{ee} (GeV)",30,60,120);
    
    //file of decalibration constants
	ofstream calFile, slopeFile;
	calFile.open("Calibration/Absolute_Consts_Cumulative_Data.txt");
	if(!calFile.is_open())
	{
		std::cout<<"Failed to create cal. constants file. Existing"<<std::endl;
		return 1;
	}
	calFile<<std::setprecision(7)<<std::fixed;
    
    char line[1000];
    int cix, ciy, iz;
    float calConst, sigma;
    ifstream oldConstsFile;
    //read in cumulative consts file
    oldConstsFile.open("/afs/cern.ch/work/a/afinkel/public/NoTrack/IC_FILES/EcalIntercalibConstants_V20120620_piZPhiSEtaScale2012_IOV2_AlphaStudies.txt");
	if(!oldConstsFile.is_open())
	{
		std::cout<<"Failed to open cumul. cal. constants file. Existing"<<std::endl;
		return 1;
	}
    do
    {
        oldConstsFile.getline(line, 1000)>>cix>>ciy>>iz>>calConst>>sigma;
        if(iz==0)continue;
        if(  (((cix-50.5)*(cix-50.5)+(ciy-50.5)*(ciy-50.5))  < 132.25) || (((cix-50.5)*(cix-50.5)+(ciy-50.5)*(ciy-50.5)) > 342.25)  ) continue;
        if(iz==-1) 
        {
            decalConstsN->SetBinContent(cix-30,ciy-30,calConst);
        }
        else if(iz==1)
        {
            decalConstsP->SetBinContent(cix-30,ciy-30,calConst);
        }
        else 
        {
            std::cout<<"Error reading in a cal. const. Exiting"<<std::endl;
            return 1;
        }
        
    } while(!oldConstsFile.eof());
    
    //and file to store absolute consts:
    ofstream absConstsFile;
    absConstsFile.open("Calibration/History/AbsoluteCumulativeConsts.txt");
    if(!absConstsFile.is_open())
    {
    	std::cout<<"Failed to create cal. constants file. Existing"<<std::endl;
		return 1;
    }
    absConstsFile<<std::setprecision(7)<<std::fixed;
			
	//grab a data ntuple
	TFile* f2 = new TFile("/afs/cern.ch/work/a/afinkel/public/NoTrack/Ntuples/Data2012/DE_ReReco_2012Full_WithClusters_001.root");	
	if(f2 == NULL)
	{
		std::cout<<"Failed to open Data file. Exiting."<<std::endl;
		return 1;
	}
    
    //gaussian fits
    TF1 *ratioFit = new TF1("fit","gaus",-1,2);
    
    //-------------------------------------------------start iterationloop here
	//std::cout<<"Ping2!"<<std::endl;
	for( int iter=1; iter<=nIterations; iter++ )
	{
		//make data ntuple:
		ZEffTree* ze2 = new ZEffTree(*f2,false);
		//Fill the data histograms:
		int ix,iy;
		int hix, hiy;
		const double Mz = 91.1876; //nominal Z mass
		double expectedPt;
		double observedPt;
		double correction;
		
		TLorentzVector elec1, elec2, theZ;
		float pt, eta, phi, E;
        
		for(int event=0; event<ze2->Entries(); event++ ) //fill the data hists
		{
			ix=ze2->reco.ix[1];//seed ix
			iy=ze2->reco.iy[1];//seed iy
            
			if( (fabs(ze2->reco.eta[1])>2.5) && (fabs(ze2->reco.eta[1])<3.0) 
				 &&(ix>29) &&(ix<71) &&(iy>29) &&(iy<71)
				 && ze2->reco.isSelected(1,"NTLooseElectronId-EtaDet") //using only events that pass selection
			  ) 
			{	
				if( (((ix-50.5)*(ix-50.5)+(iy-50.5)*(iy-50.5))  < 132.25) || (((ix-50.5)*(ix-50.5)+(iy-50.5)*(iy-50.5)) > 342.25) )
				{
					ze2->GetNextEvent();
					continue;
				}

				//compute the correction factor for observed energy
				else correction=0;
                for(unsigned int k=0; k<(ze2->ixs->size()); k++ )
                {	
                    hix = ze2->ixs->at(k);
                    hiy = ze2->iys->at(k);
                    if( (((hix-50.5)*(hix-50.5)+(hiy-50.5)*(hiy-50.5))  < 132.25) || (((hix-50.5)*(hix-50.5)+(hiy-50.5)*(hiy-50.5)) > 342.25) )
                    {
                        correction+=ze2->hitEnergyFractions->at(k);
                        continue;
                    }
                    if(ze2->reco.eta[1]>0) correction+=(constMapP[std::make_pair(hix,hiy)]->GetBinContent(iter-1))*(ze2->hitEnergyFractions->at(k))/decalConstsP->GetBinContent(hix-30,hiy-30);//positive endcap 
                    else correction+=(constMapN[std::make_pair(hix,hiy)]->GetBinContent(iter-1))*(ze2->hitEnergyFractions->at(k))/decalConstsN->GetBinContent(hix-30,hiy-30);//negative endcap
                }
				//recalculate observed Pt using correction:
				observedPt = ze2->reco.pt[1]*correction;				
				
				//recalculate Z mass using found correction
				pt = ze2->reco.pt[0];    	//tracked electron
				eta = ze2->reco.eta[0];
				phi = ze2->reco.phi[0];
				E = pt*cosh(eta);
				elec1.SetPtEtaPhiE(pt,eta,phi,E);				
				eta = ze2->reco.eta[1];
				phi = ze2->reco.phi[1];    
				E = observedPt*cosh(eta);	
				elec2.SetPtEtaPhiE(observedPt,eta,phi,E);
				theZ = elec1+elec2;		//new Z vector
				
				//fill mass hists:
				if(iter==1)
				{
					mZbefore->Fill( ze2->reco.mz );
					mZoriginal->Fill( theZ.M() );
				}
				if(iter==nIterations) mZafter->Fill( theZ.M() );
				
				//"expected energy
				expectedPt = Mz*Mz / ( 2*ze2->reco.pt[0]*( cosh(ze2->reco.eta[1]-ze2->reco.eta[0]) - cos(ze2->reco.phi[1]-ze2->reco.phi[0]) ) );
				//NOTE: this is the energy of the entire CLUSTER, not just SEED!
				//Though most of it still comes from the seed crystal... so good enough for now.				
				
				if(ze2->reco.eta[1]>0)//positive endcap
				{				
					if( (ze2->reco.mz>70) && (ze2->reco.mz<110) )
					{
						for(unsigned int k=0; k<(ze2->ixs->size()); k++ )
						{
                            hix = ze2->ixs->at(k);
							hiy = ze2->iys->at(k);
							if((hix<25)||(hix>75)||(hiy<25)||(hiy>75))
							{
								std::cout<<"Out of bounds: ("<<hix<<","<<hiy<<") Z+; frac = "<<ze2->hitEnergyFractions->at(k)<<std::endl;
								continue;
							}
							ratioMapP[std::make_pair(hix,hiy)]->Fill(expectedPt/observedPt,ze2->hitEnergyFractions->at(k));  //NOTE: Using Expected/Observed!
						}
					}
				}
				else//negative endcap
				{
					//histMapN[std::make_pair(ix,iy)]->Fill(ze2->reco.mz);
					if( (ze2->reco.mz>80) && (ze2->reco.mz<110) )
					{
						for(unsigned int k=0; k<(ze2->ixs->size()); k++ )
						{
							hix = ze2->ixs->at(k);
							hiy = ze2->iys->at(k);
							if((hix<25)||(hix>75)||(hiy<25)||(hiy>75))
							{
								std::cout<<"Out of bounds: ("<<hix<<","<<hiy<<") Z-; frac = "<<ze2->hitEnergyFractions->at(k)<<std::endl;
								continue;
							}
							ratioMapN[std::make_pair(hix,hiy)]->Fill(expectedPt/observedPt,ze2->hitEnergyFractions->at(k));  //NOTE: Using Expected/Observed!
						}
					}
				}                 
            }
			ze2->GetNextEvent();
		}//end loop over events in ntuple
		
        //do ratio and Z-peak fitting and fill 2-D hists:
        
		for(int i=25;i<76;i++)
		{
			for(int j=30;j<71;j++)
			{
				if( (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5))  < 132.25)  || (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5)) > 342.25) ) 
                {
                    constMapN[std::make_pair(i,j)]->SetBinContent(iter,1);//probably redundant, since corrections not applied in this range anyway...
                    constMapP[std::make_pair(i,j)]->SetBinContent(iter,1);
                    continue;
                }

				//do ratio fits and hists
                //negative side
                if(ratioMapN[std::make_pair(i,j)]->GetEntries() > 10)
                {
                    ratioFit->SetParameter(1,1);
                    ratioFit->SetParameter(2,0.2);
                    ratioFit->SetParLimits(1,0.5,1.5);
                    ratioFit->SetParLimits(2,0.1,0.6);
                    ratioMapN[std::make_pair(i,j)]->Fit(ratioFit,"QLLMN","",0.4,1.6);					
                    constMapN[std::make_pair(i,j)]->SetBinContent(iter,(constMapN[std::make_pair(i,j)]->GetBinContent(iter-1) )*(ratioFit->GetParameter(1) ) );
                    constMapN[std::make_pair(i,j)]->SetBinError(iter, ratioFit->GetParError(1));
                    if(iter==nIterations) 
                    {
                    	absConsts1D->Fill( constMapN[std::make_pair(i,j)]->GetBinContent(iter) );
                    	absConstsN->SetBinContent(i-30,j-30,constMapN[std::make_pair(i,j)]->GetBinContent(iter) );
                    }
                }
                else
                {
                    constMapN[std::make_pair(i,j)]->SetBinContent(iter,1 );
                    constMapN[std::make_pair(i,j)]->SetBinError(iter, 999);
                }
                //now positive side
                if(ratioMapP[std::make_pair(i,j)]->GetEntries() > 10)
                {
                    ratioFit->SetParameter(1,1);
                    ratioFit->SetParameter(2,0.2);
                    ratioFit->SetParLimits(1,0.5,1.5);
                    ratioFit->SetParLimits(2,0.1,0.6);
                    ratioMapP[std::make_pair(i,j)]->Fit(ratioFit,"QLLMN","",0.4,1.6);					
                    constMapP[std::make_pair(i,j)]->SetBinContent(iter, (constMapP[std::make_pair(i,j)]->GetBinContent(iter-1) )*(ratioFit->GetParameter(1) ));
                    constMapP[std::make_pair(i,j)]->SetBinError(iter, ratioFit->GetParError(1));      
                    if(iter==nIterations) 
                    {
                    	absConsts1D->Fill( constMapP[std::make_pair(i,j)]->GetBinContent(iter) );                  
                    	absConstsP->SetBinContent(i-30,j-30,constMapP[std::make_pair(i,j)]->GetBinContent(iter) );
                    }
                }
                else
                {
                    constMapP[std::make_pair(i,j)]->SetBinContent(iter,1 );
                    constMapP[std::make_pair(i,j)]->SetBinError(iter, 999);
                }
                if(iter != nIterations)
                {
                    ratioMapN[std::make_pair(i,j)]->Reset();
                    ratioMapP[std::make_pair(i,j)]->Reset();
                }
                else
                {
                    if(constMapP[std::make_pair(i,j)]->GetBinContent(iter)!=1)
                        absConstsFile<<i<<"\t"<<j<<"\t1\t"<<(float)(constMapP[std::make_pair(i,j)]->GetBinContent(iter))<<"\t"<<(float)ratioFit->GetParError(1)<<"\n";
                    else
                        absConstsFile<<i<<"\t"<<j<<"\t1\t"<<-1<<"\t"<<999<<"\n";
                    if(constMapN[std::make_pair(i,j)]->GetBinContent(iter)!=1)
                        absConstsFile<<i<<"\t"<<j<<"\t-1\t"<<(float)(constMapN[std::make_pair(i,j)]->GetBinContent(iter))<<"\t"<<(float)ratioFit->GetParError(1)<<"\n";
                    else
                        absConstsFile<<i<<"\t"<<j<<"\t-1\t"<<-1<<"\t"<<999<<"\n";
                }
			}
		}
        //delete the ntuple:
		delete[] ze2;
		std::cout<<iter<<" iterations complete!"<<std::endl;

	//----------------------------------------------------end iteration loop here
    }
    
    std::cout<<"Iterations done!"<<std::endl;
    //PLOTS FOR DAYS!!
    absConstsN->GetZaxis()->SetRangeUser(0.4,2);
	//absConstsN->GetZaxis()->SetLabelSize(0.03);
	absConstsN->Draw("colz");
	c1->Print("Calibration/History/absoluteCumulativeConstsN.png");
	c1->Clear();
    absConstsP->GetZaxis()->SetRangeUser(0.4,2);
	//absConstsP->GetZaxis()->SetLabelSize(0.03);
	absConstsP->Draw("colz");
	c1->Print("Calibration/History/absoluteCumulativeConstsP.png");
	c1->Clear();
	absConsts1D->SetMarkerStyle(20);
	absConsts1D->SetMarkerSize(1.4);
	absConsts1D->Draw("P");
	c1->Print("Calibration/History/absoluteCumulativeConstsDist.png");
	c1->Clear();

    TF1 *voigtFit1 = new TF1("fit1",voigt,70,110,6);
    //TF1 *bkgd = new TF1("bkgd","[0]+[1]*(x-80)",80,105 );
    
    voigtFit1->SetParameter(1,90);
    voigtFit1->SetParLimits(1,80,100);
    voigtFit1->SetParameter(2,5);
    voigtFit1->SetParLimits(2,1,20);
    voigtFit1->FixParameter(3,4.9904);
    voigtFit1->SetParameter(4,1000);
    voigtFit1->SetParameter(5,-10);
    voigtFit1->SetLineWidth(5);
    voigtFit1->SetLineColor(kBlue+1);
    mZbefore->GetXaxis()->SetRangeUser(60,120);
    mZbefore->GetYaxis()->SetRangeUser(0,mZafter->GetMaximum()*1.2);
    mZbefore->SetMarkerStyle(20);
    mZbefore->SetMarkerColor(kBlue+1);
    mZbefore->SetMarkerSize(1.4);
    mZbefore->SetLineColor(kBlue+1);
    mZbefore->SetLineWidth(5);
    mZbefore->Fit(voigtFit1,"QLMN","",70,110);
    mZbefore->Draw("P");
    voigtFit1->DrawClone("same l");
    sprintf(title,"Old Calib: #mu = %.3g, #sigma = %.3g",(float)voigtFit1->GetParameter(1),(float)voigtFit1->GetParameter(2) );
    l1->AddEntry(mZbefore,title,"l");
    
	voigtFit1->SetParameter(1,90);
    voigtFit1->SetParLimits(1,80,100);
    voigtFit1->SetParameter(2,5);
    voigtFit1->SetParLimits(2,1,20);
    voigtFit1->FixParameter(3,4.9904);
    voigtFit1->SetParameter(4,1000);
    voigtFit1->SetParameter(5,-10);
    voigtFit1->SetLineWidth(5);
    voigtFit1->SetLineColor(kRed+2);
    mZoriginal->SetMarkerStyle(20);
    mZoriginal->SetMarkerColor(kRed+2);
    mZoriginal->SetMarkerSize(1.4);
    mZoriginal->SetLineColor(kRed+2);
    mZoriginal->SetLineWidth(5);
    mZoriginal->Fit(voigtFit1,"QLMN","",70,110);
    mZoriginal->Draw("same p");
    voigtFit1->DrawClone("same l");
    sprintf(title,"Decalibrated: #mu = %.3g, #sigma = %.3g",(float)voigtFit1->GetParameter(1),(float)voigtFit1->GetParameter(2) );
    l1->AddEntry(mZoriginal,title,"l");
    
    voigtFit1->SetParameter(1,90);
    voigtFit1->SetParLimits(1,80,100);
    voigtFit1->SetParameter(2,5);
    voigtFit1->SetParLimits(2,1,20);
    voigtFit1->FixParameter(3,4.9904);
    voigtFit1->SetParameter(4,1000);
    voigtFit1->SetParameter(5,-10);
    voigtFit1->SetLineWidth(5);
    voigtFit1->SetLineColor(kGreen+3);
    mZafter->SetMarkerStyle(20);
    mZafter->SetMarkerColor(kGreen+3);
    mZafter->SetMarkerSize(1.4);
    mZafter->SetLineColor(kGreen+3);
    mZafter->SetLineWidth(5);
    mZafter->Fit(voigtFit1,"QLMN","",70,110);
    mZafter->Draw("same p");
    voigtFit1->DrawClone("same l");
    sprintf(title,"Abs Calib: #mu = %.3g, #sigma = %.3g",(float)voigtFit1->GetParameter(1),(float)voigtFit1->GetParameter(2) );
    l1->AddEntry(mZafter,title,"l");
    
    l1->Draw();
    c1->Print("Calibration/History/mZold-Decal-New.png");
    
    c1->Close();
    return 0;
}
