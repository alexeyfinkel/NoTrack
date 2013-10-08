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
#include "TLorentzVector.h"


#include "../src/ZEffTree.h"


bool isPassing(ZEffTree*);

int crystalCalibAllEE()
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
	
	//iteration number (should be about 10-15)
	const int nIterations = 15;
	
	//let's be fancy and use a map!
	std::map< pair<int,int>, TH1D* > ratioMapP, ratioMapN, constMapN, constMapP, zMassMapN, zMassMapP; 
	char title[256],name[128], filename[128];
	//and some hists
	std::vector<TH1D*> meansN, meansP, meansAll;
	std::vector<TH1D*> widthsN, widthsP, widthsAll;
	std::vector<TH1D*> calConstsN, calConstsP, calConstsAll, calConstsAllLowEta;
	std::vector<TH1D*> zMassN, zMassP, zMassAll;
	
	for(int i=0;i<101;i++)
	{
		for(int j=0;j<101;j++)
		{
			//if( (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5))  < 132.25) || (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5)) > 342.25) ) continue;
			sprintf(name,"RatioIx%diy%dzP",i,j);
			sprintf(title,"E_{exp}/E_{obs}, EE+ i_{x} = %d, i_{y} = %d;E_{exp}/E_{obs};Events",i,j);
			ratioMapP[std::make_pair(i,j)]= new TH1D(name,title,12,0.4,1.6);
			sprintf(name,"RatioIx%diy%dzN",i,j);
			sprintf(title,"E_{exp}/E_{obs}, EE- i_{x} = %d, i_{y} = %d;E_{exp}/E_{obs};Events",i,j);			
			ratioMapN[std::make_pair(i,j)]= new TH1D(name,title,12,0.4,1.6);
			sprintf(name,"CalConstIx%diy%dzN",i,j);
			sprintf(title,"Calib. Const, EE- i_{x} = %d, i_{y} = %d;Iteration;CalConst",i,j);	
			constMapN[std::make_pair(i,j)]= new TH1D(name,title,nIterations,0,nIterations);  
			sprintf(name,"CalConstIx%diy%dzP",i,j);
			sprintf(title,"Calib. Const, EE+ i_{x} = %d, i_{y} = %d;Iteration;CalConst",i,j);	
			constMapP[std::make_pair(i,j)]= new TH1D(name,title,nIterations,0,nIterations);	
			sprintf(name,"zMassIx%diy%dzN",i,j);
			sprintf(title,"M_{ee}, EE- i_{x} = %d, i_{y} = %d;M_{ee} (GeV);Events/2GeV",i,j);	
			zMassMapN[std::make_pair(i,j)]= new TH1D(name,title,30,60,120);
			sprintf(name,"zMassIx%diy%dzP",i,j);
			sprintf(title,"M_{ee}, EE+ i_{x} = %d, i_{y} = %d;M_{ee} (GeV);Events/2GeV",i,j);	
			zMassMapP[std::make_pair(i,j)]= new TH1D(name,title,30,60,120);
			
			//initialize consts	to zero, just in case (not actually used):
			constMapN[std::make_pair(i,j)]->SetBinContent(0,1);
			constMapP[std::make_pair(i,j)]->SetBinContent(0,1);
		}
	}
	//std::cout<<"Maps made; proceeding to make ntuple."<<std::endl;
	
	//initialize hists
	for( int i=1; i<=nIterations; i++ )
	{
		//mean ratios
		sprintf(name,"meansP_%d",i);
		sprintf(title,"E_{exp}/E_{obs}, EE+, %d Iterations;E_{exp}/E_{obs};Events",i);
		meansP.push_back(new TH1D(name,title,40,0.8,1.2));
		sprintf(name,"meansN_%d",i);
		sprintf(title,"E_{exp}/E_{obs}, EE-, %d Iterations;E_{exp}/E_{obs};Events",i);
		meansN.push_back(new TH1D(name,title,40,0.8,1.2));
		sprintf(name,"meansAll_%d",i);
		sprintf(title,"E_{exp}/E_{obs}, All NT, %d Iterations;E_{exp}/E_{obs};Events",i);
		meansAll.push_back(new TH1D(name,title,40,0.8,1.2));
		//widths or ratio distributions per xtl
		sprintf(name,"widthsP_%d",i);
		sprintf(title,"#sigma(E_{exp}/E_{obs}), EE+, %d Iterations;#sigma(E_{exp}/E_{obs});Events",i);
		widthsP.push_back(new TH1D(name,title,20,0.08,0.18));
		sprintf(name,"widthsN_%d",i);
		sprintf(title,"#sigma(E_{exp}/E_{obs}), EE-, %d Iterations;#sigma(E_{exp}/E_{obs});Events",i);
		widthsN.push_back(new TH1D(name,title,20,0.08,0.18));
		sprintf(name,"widthsAll_%d",i);
		sprintf(title,"#sigma(E_{exp}/E_{obs}), All NT, %d Iterations;#sigma(E_{exp}/E_{obs});Events",i);
		widthsAll.push_back(new TH1D(name,title,20,0.08,0.18));
		//calib. const. distributions
		sprintf(name,"calConstsP_%d",i);
		sprintf(title,"Cal. Consts, EE+, %d Iterations;CalConst;Events",i);
		calConstsP.push_back(new TH1D(name,title,24,0.4,1.6));
		sprintf(name,"calConstsN_%d",i);
		sprintf(title,"Cal. Consts, EE-, %d Iterations;CalConst;Events",i);
		calConstsN.push_back(new TH1D(name,title,24,0.4,1.6));
		sprintf(name,"calConstsAll_%d",i);
		sprintf(title,"Cal. Consts, All NT, %d Iterations;CalConst;Events",i);
		calConstsAll.push_back(new TH1D(name,title,24,0.4,1.6));
        sprintf(name,"calConstsAllLowEta_%d",i);
		sprintf(title,"Cal. Consts, All Low #eta, %d Iterations;CalConst;Events",i);
		calConstsAllLowEta.push_back(new TH1D(name,title,24,0.6,1.4));
		//Z mass hists
		sprintf(name,"zMassN_%d",i);
		sprintf(title,"M_{ee}, EE-, %d Iterations;M_{ee} (GeV);Events",i);
		zMassN.push_back(new TH1D(name,title,30,60,120));
		sprintf(name,"zMassP_%d",i);
		sprintf(title,"M_{ee}, EE+, %d Iterations;M_{ee} (GeV);Events",i);
		zMassP.push_back(new TH1D(name,title,30,60,120));
		sprintf(name,"zMassAll_%d",i);
		sprintf(title,"M_{ee}, All NT, %d Iterations;M_{ee} (GeV);Events",i);
		zMassAll.push_back(new TH1D(name,title,30,60,120));
	}
	//2-D hists for final Mean Ratios, Z masses before and after:
	TH2D *meansMapBeforeN = new TH2D("meansMapBeforeN","Mean Before, EE-;i_{x};i_{y}",100,0,100,100,0,100);
	TH2D *meansMapBeforeP = new TH2D("meansMapBeforeP","Mean Before, EE+;i_{x};i_{y}",100,0,100,100,0,100);
	TH2D *meansMapAfterN = new TH2D("meansMapAfterN","Mean After, EE-;i_{x};i_{y}",100,0,100,100,0,100);
	TH2D *meansMapAfterP = new TH2D("meansMapAfterP","Mean After, EE+;i_{x};i_{y}",100,0,100,100,0,100);
	TH2D *zMassBeforeN = new TH2D("zMassBeforeN","Z Mass Before Calibration, EE-;i_{x};i_{y}",100,0,100,100,0,100);
	TH2D *zMassBeforeP = new TH2D("zMassBeforeP","Z Mass Before Calibration, EE+;i_{x};i_{y}",100,0,100,100,0,100);
	TH2D *zMassAfterN = new TH2D("zMassAfterN","Z Mass After Calibration, EE-;i_{x};i_{y}",100,0,100,100,0,100);
	TH2D *zMassAfterP = new TH2D("zMassAfterP","Z Mass After Calibration, EE+;i_{x};i_{y}",100,0,100,100,0,100);
	TH2D *finalCalConstsN = new TH2D("finalCalConstsN","Final Calibration Constants, EE-;i_{x};i_{y}",100,0,100,100,0,100);
	TH2D *finalCalConstsP = new TH2D("finalCalConstsp","Final Calibration Constants, EE+;i_{x};i_{y}",100,0,100,100,0,100);
	
	
		
	//file to store calibration constants
	ofstream calFile;
	calFile.open("Calibration/AllEE/CalConsts_allEE.txt",ios::trunc);
	if(!calFile.is_open())
	{
		std::cout<<"Failed to create cal. constants file. Existing"<<std::endl;
		return 1;
	}
	calFile<<std::setprecision(7)<<std::fixed;
			
	//grab a data ntuple
	TFile* f2 = new TFile("/local/cms/user/finkel/NoTrack/Ntuple/DE_ReReco_2012Full_allEE_001.root");	
	if(f2 == NULL)
	{
		std::cout<<"Failed to open Data file. Exiting."<<std::endl;
		return 1;
	}
	
	
	//-------------------------------------------------start iterationloop here
	std::cout<<"Starting first iteration."<<std::endl;
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
			if(  isPassing(ze2)     //attempt at my combination EE+EE criterion
                 //(fabs(ze2->reco.eta[1])>2.5) && (fabs(ze2->reco.eta[1])<3.0) 
				 //&&(ix>29) &&(ix<71) &&(iy>29) &&(iy<71)
				 //&& ze2->reco.isSelected(1,"NTLooseElectronId-EtaDet") //using only events that pass selection now!
			  ) 
			{	
				if( (((ix-50.5)*(ix-50.5)+(iy-50.5)*(iy-50.5))  < 132.25) || (((ix-50.5)*(ix-50.5)+(iy-50.5)*(iy-50.5)) > 2550.25) )
				{
					ze2->GetNextEvent();
					continue;
				}

                //if(fabs(ze2->reco.eta[1])<2.5) std::cout<<"Low-eta found: Eta = "<<ze2->reco.eta[1]<<std::endl;
				//compute the correction factor for observed energy
				if(iter==1) correction = 1;
				correction=0;
                for(unsigned int k=0; k<(ze2->ixs->size()); k++ )
                {	
                    hix = ze2->ixs->at(k);
                    hiy = ze2->iys->at(k);
                    if(ze2->reco.eta[1]>0) correction+=(constMapP[std::make_pair(hix,hiy)]->GetBinContent(iter-1))*(ze2->hitEnergyFractions->at(k));//positive endcap 
                    else correction+=(constMapN[std::make_pair(hix,hiy)]->GetBinContent(iter-1))*(ze2->hitEnergyFractions->at(k));//negative endcap
                }
				//recalculate observed Pt using correction:
				observedPt = ze2->reco.pt[1]*correction;
				
				//recalculate Z mass using found correction
				pt = ze2->reco.pt[0];    	//tracked electron
				eta = ze2->reco.eta[0];
				phi = ze2->reco.phi[0];
				E = pt*cosh(eta);
				elec1.SetPtEtaPhiE(pt,eta,phi,E);				
				pt = ze2->reco.pt[1]; 		//now untracked
				pt *= correction;	
				eta = ze2->reco.eta[1];
				phi = ze2->reco.phi[1];    
				E = pt*cosh(eta);	
				elec2.SetPtEtaPhiE(pt,eta,phi,E);
				theZ = elec1+elec2;		//new Z vector
				zMassAll[iter-1]->Fill(theZ.M());
				
				//"expected energy
				expectedPt = Mz*Mz / ( 2*ze2->reco.pt[0]*( cosh(ze2->reco.eta[1]-ze2->reco.eta[0]) - cos(ze2->reco.phi[1]-ze2->reco.phi[0]) ) );
				//NOTE: this is the energy of the entire CLUSTER, not just SEED!
				//Though most of it still comes from the seed crystal... so good enough for now.				
				
				//std::cout<<"ExpectedPt = "<<expectedPt<<"; ObservedPt = "<<observedPt<<std::endl;
				if(ze2->reco.eta[1]>0)//positive endcap
				{				
					//histMapP[std::make_pair(ix,iy)]->Fill(ze2->reco.mz);
					if( (ze2->reco.mz>80) && (ze2->reco.mz<110) )
					{
						zMassP[iter-1]->Fill(theZ.M());
						zMassMapP[std::make_pair(ix,iy)]->Fill(theZ.M());
						//std::cout<<"\n\nNEXT HIT ("<<ze2->ixs->size()<<" xtl hits):"<<std::endl;
						for(unsigned int k=0; k<(ze2->ixs->size()); k++ )
						{
							hix = ze2->ixs->at(k);
							hiy = ze2->iys->at(k);
							if((hix<0)||(hix>101)||(hiy<0)||(hiy>101))
							{
								std::cout<<"Out of bounds: ("<<hix<<","<<hiy<<") Z+; frac = "<<ze2->hitEnergyFractions->at(k)<<std::endl;
								continue;
							}
							//std::cout<<"Xtal ("<<(ze2->ixs->at(k))<<","<<(ze2->iys->at(k))<<") Z+, fraction = "<<(ze2->hitEnergyFractions->at(k))<<std::endl;
							ratioMapP[std::make_pair(hix,hiy)]->Fill(expectedPt/observedPt,ze2->hitEnergyFractions->at(k));  //NOTE: Using Expected/Observed!
							
						}
					}
				}
				else//negative endcap
				{
					//histMapN[std::make_pair(ix,iy)]->Fill(ze2->reco.mz);
					if( (ze2->reco.mz>80) && (ze2->reco.mz<110) )
					{
						zMassN[iter-1]->Fill(theZ.M());
						zMassMapN[std::make_pair(ix,iy)]->Fill(theZ.M());
						//std::cout<<"\n\nNEXT HIT ("<<ze2->ixs->size()<<" xtl hits):"<<std::endl;
						for(unsigned int k=0; k<(ze2->ixs->size()); k++ )
						{
							hix = ze2->ixs->at(k);
							hiy = ze2->iys->at(k);
							if((hix<0)||(hix>101)||(hiy<0)||(hiy>101))
							{
								std::cout<<"Out of bounds: ("<<hix<<","<<hiy<<") Z-; frac = "<<ze2->hitEnergyFractions->at(k)<<std::endl;
								continue;
							}
							//if(ze2->hitEnergyFractions->at(k)<0) continue;
							//std::cout<<"Xtal ("<<(ze2->ixs->at(k))<<","<<(ze2->iys->at(k))<<") Z-, fraction = "<<(ze2->hitEnergyFractions->at(k))<<std::endl;
							ratioMapN[std::make_pair(hix,hiy)]->Fill(expectedPt/observedPt,ze2->hitEnergyFractions->at(k));  //NOTE: Using Expected/Observed!
						}
					}
				}
			}
			ze2->GetNextEvent();
		}//end loop over events in ntuple
		
		//compute next set of constants
		TF1 *ratioFit = new TF1("fit","gaus",0.4,1.6);
		TF1 *massFit = new TF1("fit","gaus",70,110);
		
        //do ratio and Z-peak fitting and fill 2-D hists:
        
		for(int i=0;i<101;i++)
		{
			for(int j=0;j<101;j++)
			{
				if( (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5))  < 132.25)  || (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5)) > 2550.25) ) continue;
				
				//do ratio fits and hists
                //negative side
				if(ratioMapN[std::make_pair(i,j)]->GetEntries() > 20)
				{
					ratioMapN[std::make_pair(i,j)]->Fit(ratioFit,"QLMN","",0.4,1.6);					
					constMapN[std::make_pair(i,j)]->SetBinContent(iter,(constMapN[std::make_pair(i,j)]->GetBinContent(iter-1) )*(ratioFit->GetParameter(1) ) );
					//constMapN[std::make_pair(i,j)]->SetBinError(iter, ratioFit->GetParError(1));
					calConstsN[iter-1]->Fill( constMapN[std::make_pair(i,j)]->GetBinContent(iter) );
					calConstsAll[iter-1]->Fill( constMapN[std::make_pair(i,j)]->GetBinContent(iter) );
                    if( ((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5)) > 342.25 )
                    {
                        calConstsAllLowEta[iter-1]->Fill( constMapN[std::make_pair(i,j)]->GetBinContent(iter) );
                    }
					meansN[iter-1]->Fill( ratioFit->GetParameter(1) );
					meansAll[iter-1]->Fill( ratioFit->GetParameter(1) );
					if(iter==1) 
                    {
                        meansMapBeforeN->SetBinContent(i,j,ratioFit->GetParameter(1));
                        zMassMapN[std::make_pair(i,j)]->Fit(massFit,"QLMN","",80,110);
                        zMassBeforeN->SetBinContent(i,j,massFit->GetParameter(1));
                    }
					if(iter==nIterations) 
					{
						meansMapAfterN->SetBinContent(i,j,ratioFit->GetParameter(1) );
						finalCalConstsN->SetBinContent(i,j,constMapN[std::make_pair(i,j)]->GetBinContent(iter) );
                        zMassMapN[std::make_pair(i,j)]->Fit(massFit,"QLMN","",75,105);
                        zMassAfterN->SetBinContent(i,j,massFit->GetParameter(1));
                        calFile<<i<<"\t"<<j<<"\t-1\t"<<(float)(constMapN[std::make_pair(i,j)]->GetBinContent(iter))<<"\t"<<(float)ratioFit->GetParError(1)<<"\n";
					}
				
					widthsN[iter-1]->Fill(ratioFit->GetParameter(2));
					widthsAll[iter-1]->Fill(ratioFit->GetParameter(2));
				}
				//sample xtl
				if(i==5 && j==50) std::cout<<"Mean = "<<ratioFit->GetParameter(1)<<", width = "<<ratioFit->GetParameter(2)<<", entries = "<<zMassMapN[std::make_pair(i,j)]->GetEntries()<<std::endl;
				
                //now positive side
				if(ratioMapP[std::make_pair(i,j)]->GetEntries() > 20)
				{
					ratioMapP[std::make_pair(i,j)]->Fit(ratioFit,"QLMN","",0.4,1.6);					
					constMapP[std::make_pair(i,j)]->SetBinContent(iter, (constMapP[std::make_pair(i,j)]->GetBinContent(iter-1) )*(ratioFit->GetParameter(1) ));
					//constMapP[std::make_pair(i,j)]->SetBinError(iter, ratioFit->GetParError(1));
					calConstsP[iter-1]->Fill( constMapP[std::make_pair(i,j)]->GetBinContent(iter) );
					calConstsAll[iter-1]->Fill( constMapP[std::make_pair(i,j)]->GetBinContent(iter) );	
                    if( ((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5)) > 342.25 )
                    {
                        calConstsAllLowEta[iter-1]->Fill( constMapP[std::make_pair(i,j)]->GetBinContent(iter) );
                    }
					meansP[iter-1]->Fill( ratioFit->GetParameter(1) );
					meansAll[iter-1]->Fill( ratioFit->GetParameter(1) );			
					if(iter==1) 
                    {
                        meansMapBeforeP->SetBinContent(i,j,ratioFit->GetParameter(1));
                        zMassMapP[std::make_pair(i,j)]->Fit(massFit,"QLMN","",80,110);
                        zMassBeforeP->SetBinContent(i,j,massFit->GetParameter(1));
                        
                    }
					if(iter==nIterations) 
					{
						meansMapAfterP->SetBinContent(i,j,ratioFit->GetParameter(1) );
						finalCalConstsP->SetBinContent(i,j,constMapP[std::make_pair(i,j)]->GetBinContent(iter) );
                        zMassMapP[std::make_pair(i,j)]->Fit(massFit,"QLMN","",75,105);
                        zMassAfterP->SetBinContent(i,j,massFit->GetParameter(1));
                        calFile<<i<<"\t"<<j<<"\t1\t"<<(float)(constMapP[std::make_pair(i,j)]->GetBinContent(iter))<<"\t"<<(float)ratioFit->GetParError(1)<<"\n";
                        
					}
				
					widthsP[iter-1]->Fill(ratioFit->GetParameter(2));
					widthsAll[iter-1]->Fill(ratioFit->GetParameter(2));
				}
				if(iter != nIterations)
				{
					ratioMapN[std::make_pair(i,j)]->Reset();
					ratioMapP[std::make_pair(i,j)]->Reset();
					zMassMapN[std::make_pair(i,j)]->Reset();
					zMassMapP[std::make_pair(i,j)]->Reset();
				}
			}
		}
		
		//draw cal. consts distributions
		sprintf(filename,"Calibration/AllEE/CalConstsByIter/CalConstsN_after%d.png",iter);
		calConstsN[iter-1]->SetMarkerStyle(20);
		calConstsN[iter-1]->Draw("P");
		c1->Print(filename);
		c1->Clear();
		sprintf(filename,"Calibration/AllEE/CalConstsByIter/CalConstsP_after%d.png",iter);
		calConstsP[iter-1]->SetMarkerStyle(20);
		calConstsP[iter-1]->Draw("P");
		c1->Print(filename);
		c1->Clear();
		sprintf(filename,"Calibration/AllEE/CalConstsByIter/CalConstsAll_after%d.png",iter);
		calConstsAll[iter-1]->SetMarkerStyle(20);
		calConstsAll[iter-1]->Draw("P");
		c1->Print(filename);
		c1->Clear();
        sprintf(filename,"Calibration/AllEE/CalConstsByIter/CalConstsAllLowEta_after%d.png",iter);
		calConstsAllLowEta[iter-1]->SetMarkerStyle(20);
		calConstsAllLowEta[iter-1]->Draw("P");
		c1->Print(filename);
		c1->Clear();
		//draw widths distributions
		sprintf(filename,"Calibration/AllEE/WidthsByIter/WidthsN_after%d.png",iter);
		widthsN[iter-1]->SetMarkerStyle(20);
		widthsN[iter-1]->Draw("P");
		c1->Print(filename);
		c1->Clear();
		sprintf(filename,"Calibration/AllEE/WidthsByIter/WidthsP_after%d.png",iter);
		widthsP[iter-1]->SetMarkerStyle(20);
		widthsP[iter-1]->Draw("P");
		c1->Print(filename);
		c1->Clear();
		sprintf(filename,"Calibration/AllEE/WidthsByIter/WidthsAll_after%d.png",iter);
		widthsAll[iter-1]->SetMarkerStyle(20);
		widthsAll[iter-1]->Draw("P");
		c1->Print(filename);
		c1->Clear();
		//draw Z mass plots
		sprintf(filename,"Calibration/AllEE/ZMassByIter/ZMassN_after%d.png",iter);
		zMassN[iter-1]->SetMarkerStyle(20);
		zMassN[iter-1]->Draw("P");
		c1->Print(filename);
		c1->Clear();
		sprintf(filename,"Calibration/AllEE/ZMassByIter/ZMassP_after%d.png",iter);
		zMassP[iter-1]->SetMarkerStyle(20);
		zMassP[iter-1]->Draw("P");
		c1->Print(filename);
		c1->Clear();
		sprintf(filename,"Calibration/AllEE/ZMassByIter/ZMassAll_after%d.png",iter);
		zMassAll[iter-1]->SetMarkerStyle(20);
		zMassAll[iter-1]->Draw("P");
		c1->Print(filename);
		c1->Clear();
		
		//delete the ntuple:
		delete[] ze2;
		std::cout<<iter<<" iterations complete!"<<std::endl;
	}
	
	//----------------------------------------------------end iteration loop here

	//draw sample plot
	ratioMapN[std::make_pair(43,36)]->SetMarkerStyle(20);
	ratioMapN[std::make_pair(43,36)]->Draw("P");
	c1->Print("Calibration/AllEE/RatioPlots/Ratio_N_ix48_iy38.png");
	c1->Clear();
	constMapN[std::make_pair(48,38)]->SetMarkerStyle(20);
	constMapN[std::make_pair(48,38)]->Draw("P");
	c1->Print("Calibration/AllEE/ConvergencePlots/CalConst_N_ix48_iy38.png");
	c1->Clear();
	//initial means
	meansMapBeforeN->GetZaxis()->SetRangeUser(0.8,1.2);
	meansMapBeforeN->GetZaxis()->SetLabelSize(0.03);
	meansMapBeforeN->Draw("colz");
	c1->Print("Calibration/AllEE/InitialMeansN.png");
	c1->Clear();
	meansMapBeforeP->GetZaxis()->SetRangeUser(0.8,1.2);
	meansMapBeforeP->GetZaxis()->SetLabelSize(0.03);
	meansMapBeforeP->Draw("colz");
	c1->Print("Calibration/AllEE/InitialMeansP.png");
	c1->Clear();
	//final means
	meansMapAfterN->GetZaxis()->SetRangeUser(0.8,1.2);
	meansMapAfterN->GetZaxis()->SetLabelSize(0.03);
	meansMapAfterN->Draw("colz");
	c1->Print("Calibration/AllEE/FinalMeansN.png");
	c1->Clear();
	meansMapAfterP->GetZaxis()->SetRangeUser(0.8,1.2);
	meansMapAfterP->GetZaxis()->SetLabelSize(0.03);
	meansMapAfterP->Draw("colz");
	c1->Print("Calibration/AllEE/FinalMeansP.png");
	c1->Clear();
	//Z mass before
	zMassBeforeN->GetZaxis()->SetRangeUser(82,102);
	zMassBeforeN->GetZaxis()->SetLabelSize(0.03);
	zMassBeforeN->Draw("colz");
	c1->Print("Calibration/AllEE/ZMassBeforeN.png");
	c1->Clear();
	zMassBeforeP->GetZaxis()->SetRangeUser(82,102);
	zMassBeforeP->GetZaxis()->SetLabelSize(0.03);
	zMassBeforeP->Draw("colz");
	c1->Print("Calibration/AllEE/ZMassBeforeP.png");
	c1->Clear();
	//Z mass after
	zMassAfterN->GetZaxis()->SetRangeUser(82,102);
	zMassAfterN->GetZaxis()->SetLabelSize(0.03);
	zMassAfterN->Draw("colz");
	c1->Print("Calibration/AllEE/ZMassAfterN.png");
	c1->Clear();
	zMassAfterP->GetZaxis()->SetRangeUser(82,102);
	zMassAfterP->GetZaxis()->SetLabelSize(0.03);
	zMassAfterP->Draw("colz");
	c1->Print("Calibration/AllEE/ZMassAfterP.png");
	c1->Clear();
	//final calibration constants
	finalCalConstsN->GetZaxis()->SetRangeUser(0.8,1.2);
	finalCalConstsN->GetZaxis()->SetLabelSize(0.03);
	finalCalConstsN->Draw("colz");
	c1->Print("Calibration/AllEE/FinalCalConstsN.png");
	c1->Clear();
	finalCalConstsP->GetZaxis()->SetRangeUser(0.8,1.2);
	finalCalConstsP->GetZaxis()->SetLabelSize(0.03);
	finalCalConstsP->Draw("colz");
	c1->Print("Calibration/AllEE/FinalCalConstsP.png");
	c1->Clear();
		
	c1->Close();
	return 0;
}


bool isPassing(ZEffTree* ze)
{
    if(  ze->reco.isSelected(1,"NTLooseElectronId-EtaDet")) return true;
    if(  ze->reco.Sieie[1]<0.029
      && ze->reco.RNine[1]<1.02 && ze->reco.RNine[1]>0.89
      && ze->reco.Eiso[1]<0.035
      && ze->reco.Hiso[1]<0.11
      && ze->reco.HoEM[1]<0.05
      && ze->reco.mz>60 && ze->reco.mz<120 ) return true;
    return false;
}