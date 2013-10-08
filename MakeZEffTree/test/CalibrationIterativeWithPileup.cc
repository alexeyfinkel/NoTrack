//this script is supposed to run iterative process in PU bind and project back to 0 PU.

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

Double_t voigt(Double_t *x, Double_t *par)
{
    return par[0]*(TMath::Voigt((x[0]-par[1]),par[2],par[3]));
}

int calibrate()
{
    gROOT->SetStyle("Plain");	
	gStyle->SetErrorX(0);
	gStyle->SetOptStat(0);
	
	TCanvas* c1 = new TCanvas("C1", "c1", 1000, 800);
	c1->cd();
	//c1->SetLeftMargin(0.08);
	//c1->SetRightMargin(0.12);

	//iteration number (should be about 10-15)
	const int nIterations = 15;
    
    //let's be fancy and use a map!
	std::map< triplet, TH1D* > ratioMapP, ratioMapN, constMapN, constMapP, zMassMapN, zMassMapP;
    std::map< std::pair<int,int>, TH1D* > constByPuN, constByPuP;
    char title[256],name[128], filename[128];
    
    for(int i=25;i<76;i++)
	{
		for(int j=25;j<76;j++)
		{
            sprintf(name,"constsByPU_ix%diy%dzP",i,j);
            sprintf(title,"Const by PU, NT+ i_{x} = %d, i_{y} = %d;PU;Const",i,j);
            constByPuP[std::make_pair(i,j)]=new TH1D(name,title,4,5,25);
            sprintf(name,"constsByPU_ix%diy%dzN",i,j);
            sprintf(title,"Const by PU, NT- i_{x} = %d, i_{y} = %d;PU;Const",i,j);
            constByPuN[std::make_pair(i,j)]=new TH1D(name,title,4,5,25);
            
            for(int k=0;k<4;k++)
			{
                //if( (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5))  < 132.25) || (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5)) > 342.25) ) continue;
                sprintf(name,"RatioIx%diy%dzP_pu%d",i,j,k);
                sprintf(title,"E_{exp}/E_{obs}, NT+ i_{x} = %d, i_{y} = %d, PU %d-%d;E_{exp}/E_{obs};Events",i,j,5*(k+1),5*(k+2));
                ratioMapP[triplet(i,j,k)]= new TH1D(name,title,24,0.4,1.6);
                sprintf(name,"RatioIx%diy%dzN_pu%d",i,j,k);
                sprintf(title,"E_{exp}/E_{obs}, NT- i_{x} = %d, i_{y} = %d, PU %d-%d;E_{exp}/E_{obs};Events",i,j,5*(k+1),5*(k+2));			
                ratioMapN[triplet(i,j,k)]= new TH1D(name,title,24,0.4,1.6);
                sprintf(name,"CalConstIx%diy%dzN_pu%d",i,j,k);
                sprintf(title,"Calib. Const, NT- i_{x} = %d, i_{y} = %d, PU %d-%d;Iteration;CalConst",i,j,5*(k+1),5*(k+2));	
                constMapN[triplet(i,j,k)]= new TH1D(name,title,nIterations,0,nIterations);  
                sprintf(name,"CalConstIx%diy%dzP_pu%d",i,j,k);
                sprintf(title,"Calib. Const, NT+ i_{x} = %d, i_{y} = %d, PU %d-%d;Iteration;CalConst",i,j,5*(k+1),5*(k+2));	
                constMapP[triplet(i,j,k)]= new TH1D(name,title,nIterations,0,nIterations);	
                sprintf(name,"zMassIx%diy%dzN_pu%d",i,j,k);
                sprintf(title,"M_{ee}, NT- i_{x} = %d, i_{y} = %d, PU %d-%d;M_{ee} (GeV);Events/2GeV",i,j,5*(k+1),5*(k+2));	
                zMassMapN[triplet(i,j,k)]= new TH1D(name,title,30,60,120);
                sprintf(name,"zMassIx%diy%dzP_pu%d",i,j,k);
                sprintf(title,"M_{ee}, NT+ i_{x} = %d, i_{y} = %d, PU %d-%d;M_{ee} (GeV);Events/2GeV",i,j,5*(k+1),5*(k+2));	
                zMassMapP[triplet(i,j,k)]= new TH1D(name,title,30,60,120);

                //initialize consts	to zero, just in case (not actually used):
                constMapN[triplet(i,j,k)]->SetBinContent(0,1);
                constMapP[triplet(i,j,k)]->SetBinContent(0,1);
            }
		}
	}
    //Z mass by pileup
    std::vector<TH1D*> ZbeforeByPU, ZafterByPU;
    for(int bin=0;bin<4;bin++)
    {
        sprintf(name,"Zbefore_pu%d",bin);
		sprintf(title,"M_{Z} Before, PU %d-%d;M_{ee};Events",5*(bin+1),5*(bin+2));
        ZbeforeByPU.push_back(new TH1D(name,title,30,60,120));
        sprintf(name,"Zafter_pu%d",bin);
		sprintf(title,"M_{Z} After, PU %d-%d;M_{ee};Events",5*(bin+1),5*(bin+2));
        ZafterByPU.push_back(new TH1D(name,title,30,60,120));
    }
    TH1D *ZmassVsPU = new TH1D("zMassVsPU","Z Mass vs PU,;nVert;M_{peak}",4,5,25);
    TH1D *Slopes1D = new TH1D("Slopes1D","PU-dependence slopes,;slope;crystals",20,-0.02,0.02);
    TH1D *Ratios1D = new TH1D("Ratios1D","Projected/Cumulative Distribution;C_{0}/C_{tot};crystals",20,0.6,1.5);
    TH1D *slopesVsEtaN = new TH1D("slopesVsEtaN","Ring-avg. slope vs. #eta, NT-;#eta-ring;slope",7,0,7);
    TH1D *slopesVsEtaP = new TH1D("slopesVsEtaP","Ring-avg. slope vs. #eta, NT+;#eta-ring;slope",7,0,7);
    TH1D *ratiosVsEtaN = new TH1D("ratiosVsEtaN","Ring-avg. ratios vs. #eta, NT-;#eta-ring;slope",7,0,7);
    TH1D *ratiosVsEtaP = new TH1D("ratiosVsEtaP","Ring-avg. ratios vs. #eta, NT+;#eta-ring;slope",7,0,7);
    
    //ring hists for plots and ratios
    std::vector<TH1D*> slopesByRingN, slopesByRingP, ratiosByRingN, ratiosByRingP;
    for(int i=0;i<7;i++)
    {
        sprintf(name,"slopesRing%dN",i);
        sprintf(title,"Slopes in NT- #eta-ring %d;#eta-ring;slope",i+1);
        slopesByRingN.push_back(new TH1D(name,title,10,-0.05,0.05));
        sprintf(name,"slopesRing%dP",i);
        sprintf(title,"Slopes in NT+ #eta-ring %d;#eta-ring;slope",i+1);
        slopesByRingP.push_back(new TH1D(name,title,10,-0.05,0.05));
        sprintf(name,"ratiosRing%dN",i);
        sprintf(title,"Ratios in NT- #eta-ring %d;#eta-ring;slope",i+1);
        ratiosByRingN.push_back(new TH1D(name,title,10,0.7,1.4));
        sprintf(name,"ratiosRing%dP",i);
        sprintf(title,"Ratios in NT+ #eta-ring %d;#eta-ring;slope",i+1);
        ratiosByRingP.push_back(new TH1D(name,title,10,0.7,1.4));
    }
    
    
    //2-D hists for final Mean Ratios, Z masses before and after:
	//TH2D *meansMapBeforeN = new TH2D("meansMapBeforeN","Mean Before, NT-;i_{x};i_{y}",40,30,70,40,30,70);
	//TH2D *meansMapBeforeP = new TH2D("meansMapBeforeP","Mean Before, NT+;i_{x};i_{y}",40,30,70,40,30,70);
	//TH2D *meansMapAfterN = new TH2D("meansMapAfterN","Mean After, NT-;i_{x};i_{y}",40,30,70,40,30,70);
	//TH2D *meansMapAfterP = new TH2D("meansMapAfterP","Mean After, NT+;i_{x};i_{y}",40,30,70,40,30,70);
	TH2D *zMassBeforeN = new TH2D("zMassBeforeN","Z Mass Before Calibration, NT-;i_{x};i_{y}",40,30,70,40,30,70);
	TH2D *zMassBeforeP = new TH2D("zMassBeforeP","Z Mass Before Calibration, NT+;i_{x};i_{y}",40,30,70,40,30,70);
	TH2D *zMassAfterN = new TH2D("zMassAfterN","Z Mass After Calibration, NT-;i_{x};i_{y}",40,30,70,40,30,70);
	TH2D *zMassAfterP = new TH2D("zMassAfterP","Z Mass After Calibration, NT+;i_{x};i_{y}",40,30,70,40,30,70);
	TH2D *finalCalConstsN = new TH2D("finalCalConstsN","Final Calibration Constants, NT-;i_{x};i_{y}",40,30,70,40,30,70);
	TH2D *finalCalConstsP = new TH2D("finalCalConstsp","Final Calibration Constants, NT+;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *projectedCalConstsN = new TH2D("projectedCalConstsN","Projected Calibration Constants, NT-;i_{x};i_{y}",40,30,70,40,30,70);
	TH2D *projectedCalConstsP = new TH2D("projectedCalConstsp","Projected Calibration Constants, NT+;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *Slopes2DN = new TH2D("Slopes2DN","PU-dependence slopes, NT-,;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *Slopes2DP = new TH2D("Slopes2DP","PU-dependence slopes, NT+,;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *projectionRatiosN = new TH2D("projectionRatiosN","Projected/Cumulative, NT-;i_{x};i_{y}",40,30,70,40,30,70);
	TH2D *projectionRatiosP = new TH2D("projectionRatiosP","Projected/Cumulative, NT+;i_{x};i_{y}",40,30,70,40,30,70);

    
    //file to store calibration constants
	ofstream calFile;
	calFile.open("Calibration/PU/CalConsts_PUregress_Data.txt");
	if(!calFile.is_open())
	{
		std::cout<<"Failed to create cal. constants file. Existing"<<std::endl;
		return 1;
	}
	calFile<<std::setprecision(7)<<std::fixed;
			
	//grab a data ntuple
	TFile* f2 = new TFile("/local/cms/user/finkel/NoTrack/Ntuple/DE_ReReco_2012Full_WithClusters_001.root");	
	if(f2 == NULL)
	{
		std::cout<<"Failed to open Data file. Exiting."<<std::endl;
		return 1;
	}
    
    //gaussian fits
    TF1 *ratioFit = new TF1("fit","gaus",0.4,1.6);
	TF1 *massFit = new TF1("fit","gaus",70,110);
    
    //tracking averages:
    std::vector<int> evtsInBin(4,0), nVertInBin(4,0);
    std::vector<double> binAveragePU(4,0), zPeaks(4,0);
    //-------------------------------------------------start iterationloop here
	//std::cout<<"Ping2!"<<std::endl;
	for( int iter=1; iter<=nIterations; iter++ )
	{
		//make data ntuple:
		ZEffTree* ze2 = new ZEffTree(*f2,false);
		//Fill the data histograms:
		int ix,iy;
		int hix, hiy;
        int nPu, puBin;
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
            nPu=ze2->reco.nverts-1;
            if(nPu>=25) puBin=3;
            else if(nPu<=5 ) puBin=0;
            else puBin=(int)((nPu-5)/5);
            
            if(iter==1)
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
				if(iter==1) correction = 1;
				correction=0;
                for(unsigned int k=0; k<(ze2->ixs->size()); k++ )
                {	
                    hix = ze2->ixs->at(k);
                    hiy = ze2->iys->at(k);
                    if(ze2->reco.eta[1]>0) correction+=(constMapP[triplet(hix,hiy,puBin)]->GetBinContent(iter-1))*(ze2->hitEnergyFractions->at(k));//positive endcap 
                    else correction+=(constMapN[triplet(hix,hiy,puBin)]->GetBinContent(iter-1))*(ze2->hitEnergyFractions->at(k));//negative endcap
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
				//zMassAll[iter-1]->Fill(theZ.M());
				
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
						//zMassP[iter-1]->Fill(theZ.M());
						//zMassMapP[std::make_pair(ix,iy)]->Fill(theZ.M());
						//std::cout<<"\n\nNEXT HIT ("<<ze2->ixs->size()<<" xtl hits):"<<std::endl;
						for(unsigned int k=0; k<(ze2->ixs->size()); k++ )
						{
                            hix = ze2->ixs->at(k);
							hiy = ze2->iys->at(k);
							if((hix<25)||(hix>75)||(hiy<25)||(hiy>75))
							{
								std::cout<<"Out of bounds: ("<<hix<<","<<hiy<<") Z+; frac = "<<ze2->hitEnergyFractions->at(k)<<std::endl;
								continue;
							}
							//std::cout<<"Xtal ("<<(ze2->ixs->at(k))<<","<<(ze2->iys->at(k))<<") Z+, fraction = "<<(ze2->hitEnergyFractions->at(k))<<std::endl;
							ratioMapP[triplet(hix,hiy,puBin)]->Fill(expectedPt/observedPt,ze2->hitEnergyFractions->at(k));  //NOTE: Using Expected/Observed!
							
						}
					}
				}
				else//negative endcap
				{
					//histMapN[std::make_pair(ix,iy)]->Fill(ze2->reco.mz);
					if( (ze2->reco.mz>80) && (ze2->reco.mz<110) )
					{
						//zMassN[iter-1]->Fill(theZ.M());
						//zMassMapN[triplet(ix,iy,puBin)]->Fill(theZ.M());
						//std::cout<<"\n\nNEXT HIT ("<<ze2->ixs->size()<<" xtl hits):"<<std::endl;
						for(unsigned int k=0; k<(ze2->ixs->size()); k++ )
						{
							hix = ze2->ixs->at(k);
							hiy = ze2->iys->at(k);
							if((hix<25)||(hix>75)||(hiy<25)||(hiy>75))
							{
								std::cout<<"Out of bounds: ("<<hix<<","<<hiy<<") Z-; frac = "<<ze2->hitEnergyFractions->at(k)<<std::endl;
								continue;
							}
							//if(ze2->hitEnergyFractions->at(k)<0) continue;
							//std::cout<<"Xtal ("<<(ze2->ixs->at(k))<<","<<(ze2->iys->at(k))<<") Z-, fraction = "<<(ze2->hitEnergyFractions->at(k))<<std::endl;
							ratioMapN[triplet(hix,hiy,puBin)]->Fill(expectedPt/observedPt,ze2->hitEnergyFractions->at(k));  //NOTE: Using Expected/Observed!
						}
					}
				}
                
                if(iter==1) ZbeforeByPU[puBin]->Fill(ze2->reco.mz);
                if(iter==nIterations) ZafterByPU[puBin]->Fill(theZ.M());                    
            }
			ze2->GetNextEvent();
		}//end loop over events in ntuple
		
        //do ratio and Z-peak fitting and fill 2-D hists:
        
		for(int i=25;i<76;i++)
		{
			for(int j=30;j<71;j++)
			{
				if( (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5))  < 132.25)  || (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5)) > 342.25) ) continue;
				for(int k=0;k<4;k++)
                {
                    //do ratio fits and hists
                    //negative side
                    if(ratioMapN[triplet(i,j,k)]->GetEntries() > 20)
                    {
                        ratioMapN[triplet(i,j,k)]->Fit(ratioFit,"QLMN","",0.4,1.6);					
                        constMapN[triplet(i,j,k)]->SetBinContent(iter,(constMapN[triplet(i,j,k)]->GetBinContent(iter-1) )*(ratioFit->GetParameter(1) ) );
                        constMapN[triplet(i,j,k)]->SetBinError(iter, ratioFit->GetParError(1));
                        //calConstsN[iter-1]->Fill( constMapN[triplet(i,j,k)]->GetBinContent(iter) );
                        //calConstsAll[iter-1]->Fill( constMapN[triplet(i,j,k)]->GetBinContent(iter) );
                        //meansN[iter-1]->Fill( ratioFit->GetParameter(1) );
                        //meansAll[iter-1]->Fill( ratioFit->GetParameter(1) );
                        if(iter==1) 
                        {
                            zMassMapN[triplet(i,j,k)]->Fit(massFit,"QLMN","",80,110);
                            zMassBeforeN->SetBinContent(i-30,j-30,massFit->GetParameter(1));
                        }
                        if(iter==nIterations) 
                        {
                            //meansMapAfterN->SetBinContent(i-30,j-30,ratioFit->GetParameter(1) );
                            finalCalConstsN->SetBinContent(i-30,j-30,constMapN[triplet(i,j,k)]->GetBinContent(iter) );
                            zMassMapN[triplet(i,j,k)]->Fit(massFit,"QLMN","",75,105);
                            zMassAfterN->SetBinContent(i-30,j-30,massFit->GetParameter(1));
                            calFile<<i<<"\t"<<j<<"\t-1\t"<<(float)(constMapN[triplet(i,j,k)]->GetBinContent(iter))<<"\t"<<(float)ratioFit->GetParError(1)<<"\n";
                        }
                    }
                    //now positive side
                    if(ratioMapP[triplet(i,j,k)]->GetEntries() > 20)
                    {
                        ratioMapP[triplet(i,j,k)]->Fit(ratioFit,"QLMN","",0.4,1.6);					
                        constMapP[triplet(i,j,k)]->SetBinContent(iter, (constMapP[triplet(i,j,k)]->GetBinContent(iter-1) )*(ratioFit->GetParameter(1) ));
                        constMapP[triplet(i,j,k)]->SetBinError(iter, ratioFit->GetParError(1));
                        //calConstsP[iter-1]->Fill( constMapP[triplet(i,j,k)]->GetBinContent(iter) );
                        //calConstsAll[iter-1]->Fill( constMapP[triplet(i,j,k)]->GetBinContent(iter) );	
                        //meansP[iter-1]->Fill( ratioFit->GetParameter(1) );
                        //meansAll[iter-1]->Fill( ratioFit->GetParameter(1) );			
                        if(iter==1) 
                        {
                            //meansMapBeforeP->SetBinContent(i-30,j-30,ratioFit->GetParameter(1));
                            zMassMapP[triplet(i,j,k)]->Fit(massFit,"QLMN","",80,110);
                            zMassBeforeP->SetBinContent(i-30,j-30,massFit->GetParameter(1));

                        }
                        if(iter==nIterations) 
                        {
                            //meansMapAfterP->SetBinContent(i-30,j-30,ratioFit->GetParameter(1) );
                            finalCalConstsP->SetBinContent(i-30,j-30,constMapP[triplet(i,j,k)]->GetBinContent(iter) );
                            zMassMapP[triplet(i,j,k)]->Fit(massFit,"QLMN","",75,105);
                            zMassAfterP->SetBinContent(i-30,j-30,massFit->GetParameter(1));
                            calFile<<i<<"\t"<<j<<"\t1\t"<<(float)(constMapP[triplet(i,j,k)]->GetBinContent(iter))<<"\t"<<(float)ratioFit->GetParError(1)<<"\n";

                        }
                    }
                    if(iter != nIterations)
                    {
                        ratioMapN[triplet(i,j,k)]->Reset();
                        ratioMapP[triplet(i,j,k)]->Reset();
                        zMassMapN[triplet(i,j,k)]->Reset();
                        zMassMapP[triplet(i,j,k)]->Reset();
                    }
                }
			}
		}
        //delete the ntuple:
		delete[] ze2;
		std::cout<<iter<<" iterations complete!"<<std::endl;

	//----------------------------------------------------end iteration loop here
    }
    
    std::cout<<"Iterations done!"<<std::endl;

    TF1 *puFit = new TF1("fit","pol1",5,25);
    //get final consts by PU, cram them into constByPU 
    int ring;
    
    for(int i=25;i<76;i++)
	{
		for(int j=25;j<76;j++)
		{
            projectedCalConstsP->SetBinContent(i-30,j-30, -1 );
            projectedCalConstsN->SetBinContent(i-30,j-30, -1 );
            Slopes2DN->SetBinContent(i-30,j-30, -1 );
            Slopes2DP->SetBinContent(i-30,j-30, -1 );
            if( (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5))  < 132.25)  || (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5)) > 342.25) ) continue;
            ring = (int)(sqrt( (i-50.5)*(i-50.5)+(j-50.5)*(j-50.5) )-11.5);
            if(ring<0 || ring >6) std::cout<<"Invalid Ring: "<<ring<<std::endl;
            for(int k=0;k<4;k++)
            {
                constByPuP[std::make_pair(i,j)]->SetBinContent(k+1,constMapP[triplet(i,j,k)]->GetBinContent(nIterations));
                constByPuP[std::make_pair(i,j)]->SetBinError(k+1,constMapP[triplet(i,j,k)]->GetBinError(nIterations));
                constByPuN[std::make_pair(i,j)]->SetBinContent(k+1,constMapN[triplet(i,j,k)]->GetBinContent(nIterations));
                constByPuN[std::make_pair(i,j)]->SetBinError(k+1,constMapN[triplet(i,j,k)]->GetBinError(nIterations));
            }
            if(constByPuP[std::make_pair(i,j)]->GetEntries()==4)
            {
                puFit->SetParameter(0,1);
                puFit->SetParameter(1,0);
                //puFit->SetParameter(2,0);
                //constByPuP[std::make_pair(i,j)]->Draw("E");
                //constByPuP[std::make_pair(i,j)]->GetYaxis()->SetRangeUser(0.6,1.4);
                constByPuP[std::make_pair(i,j)]->Fit(puFit,"QLMN","",5,25);
                //sprintf(title,"Calibration/PU/CrystalFits/ix%d_iy%d_P.png",i,j);
                //c1->Print(title);
                projectedCalConstsP->SetBinContent(i-30,j-30, puFit->GetParameter(0) );
                Slopes2DP->SetBinContent(i-30,j-30, puFit->GetParameter(1) );
                projectionRatiosP->SetBinContent(i-30,j-30, puFit->GetParameter(0)/finalCalConstsP->GetBinContent(i-30,j-30) );
                Slopes1D->Fill(puFit->GetParameter(1));
                Ratios1D->Fill(projectionRatiosP->GetBinContent(i-30,j-30));
                slopesByRingP[ring]->Fill(puFit->GetParameter(1));
                ratiosByRingP[ring]->Fill(puFit->GetParameter(0)/finalCalConstsP->GetBinContent(i-30,j-30));
            }
            else std::cout<<"Bad Fit: ix = "<<i<<", iy = "<<j<<", iz = 1"<<std::endl;
            if(constByPuN[std::make_pair(i,j)]->GetEntries()==4)
            {
                puFit->SetParameter(0,1);
                puFit->SetParameter(1,0);
                //puFit->SetParameter(2,0);
                //constByPuN[std::make_pair(i,j)]->Draw("E");
                //constByPuN[std::make_pair(i,j)]->GetYaxis()->SetRangeUser(0.6,1.4);
                constByPuN[std::make_pair(i,j)]->Fit(puFit,"QLMN","",5,25);
                //sprintf(title,"Calibration/PU/CrystalFits/ix%d_iy%d_N.png",i,j);
                //c1->Print(title);
                projectedCalConstsN->SetBinContent(i-30,j-30, puFit->GetParameter(0) );
                Slopes2DN->SetBinContent(i-30,j-30, puFit->GetParameter(1) );
                projectionRatiosN->SetBinContent(i-30,j-30, puFit->GetParameter(0)/finalCalConstsN->GetBinContent(i-30,j-30) );
                Slopes1D->Fill(puFit->GetParameter(1));
                Ratios1D->Fill(projectionRatiosN->GetBinContent(i-30,j-30));
                slopesByRingN[ring]->Fill(puFit->GetParameter(1));
                ratiosByRingN[ring]->Fill(puFit->GetParameter(0)/finalCalConstsN->GetBinContent(i-30,j-30));
            }
            else std::cout<<"Bad Fit: ix = "<<i<<", iy = "<<j<<", iz = -1"<<std::endl;
        }
    }

    //std::cout<<"Ping!"<<std::endl;
    for(int r=1;r<8;r++)
    {
        slopesVsEtaN->SetBinContent(r,slopesByRingN[r-1]->GetMean());
        slopesVsEtaP->SetBinContent(r,slopesByRingP[r-1]->GetMean());
        ratiosVsEtaN->SetBinContent(r,ratiosByRingN[r-1]->GetMean());
        ratiosVsEtaP->SetBinContent(r,ratiosByRingP[r-1]->GetMean());
    }
    
     //projected calibration constants
	projectedCalConstsN->GetZaxis()->SetRangeUser(0.75,1.25);
	projectedCalConstsN->GetZaxis()->SetLabelSize(0.03);
	projectedCalConstsN->Draw("colz");
	c1->Print("Calibration/PU/ProjectedConstsN.png");
	c1->Clear();
	projectedCalConstsP->GetZaxis()->SetRangeUser(0.75,1.25);
	projectedCalConstsP->GetZaxis()->SetLabelSize(0.03);
	projectedCalConstsP->Draw("colz");
	c1->Print("Calibration/PU/ProjectedCalConstsP.png");
	c1->Clear();
    Slopes2DN->GetZaxis()->SetRangeUser(-0.01,0.01);
	Slopes2DN->GetZaxis()->SetLabelSize(0.03);
	Slopes2DN->Draw("colz");
	c1->Print("Calibration/PU/Slopes2DN.png");
	c1->Clear();
    Slopes2DP->GetZaxis()->SetRangeUser(-0.01,0.01);
	Slopes2DP->GetZaxis()->SetLabelSize(0.03);
	Slopes2DP->Draw("colz");
	c1->Print("Calibration/PU/Slopes2DP.png");
	c1->Clear();
    projectionRatiosP->GetZaxis()->SetRangeUser(0.85,1.2);
	projectionRatiosP->GetZaxis()->SetLabelSize(0.03);
	projectionRatiosP->Draw("colz");
	c1->Print("Calibration/PU/ProjectionRatiosP.png");
	c1->Clear();
    projectionRatiosN->GetZaxis()->SetRangeUser(0.85,1.2);
	projectionRatiosN->GetZaxis()->SetLabelSize(0.03);
	projectionRatiosN->Draw("colz");
	c1->Print("Calibration/PU/ProjectionRatiosN.png");
	c1->Clear();
   
    TF1 *voigtFit1 = new TF1("fit1",voigt,80,105,4);
    
    //c1->SetMargin(0.08,0.03,0.09,0.05);
    TLegend* l1 = new TLegend(0.7,0.7,0.9,0.9);
	l1->SetFillColor(10);
    
    
    
    for(int k=0;k<4;k++)
    {
        sprintf(title,"Before Calibration, PU %d-%d;M_{ee} (GeV);Events/GeV",5*(k+1),5*(k+2) );
        ZbeforeByPU[k]->SetTitle(title);
        //ZbeforeByPU[k]->GetYaxis()->SetTitleOffset(1);
        //ZbeforeByPU[k]->GetYaxis()->SetTitleSize(0.04);
        ZbeforeByPU[k]->GetXaxis()->SetRangeUser(70,110);
        ZbeforeByPU[k]->SetMarkerStyle(21);
        ZbeforeByPU[k]->SetLineColor(kBlue);
        ZbeforeByPU[k]->SetMarkerColor(kBlue+3);
        ZbeforeByPU[k]->SetMarkerSize(1.3);
        voigtFit1->SetParameter(1,90);
        voigtFit1->SetParLimits(1,80,100);
        voigtFit1->SetParameter(2,5);
        voigtFit1->SetParLimits(2,1,20);
        voigtFit1->FixParameter(3,4.9904);
        voigtFit1->SetLineWidth(4);
        voigtFit1->SetLineColor(kBlue);
        ZbeforeByPU[k]->Fit(voigtFit1,"QLM","",80,105);
        sprintf(title,"#mu = %.3g, #sigma = %.3g",(float)voigtFit1->GetParameter(1),(float)voigtFit1->GetParameter(2) );
        l1->AddEntry(voigtFit1,title,"");
        sprintf(filename,"Calibration/PU/ZmassBefore_pu%d.png",k);
        ZbeforeByPU[k]->Draw("P");
        l1->Draw();
        c1->Print(filename);
        l1->Clear();
        c1->Clear();
        
        ZmassVsPU->SetBinContent(k+1,voigtFit1->GetParameter(1));
        ZmassVsPU->SetBinError(k+1,voigtFit1->GetParError(1));
        
        binAveragePU[k]=(double)nVertInBin[k]/(double)evtsInBin[k];
        zPeaks[k]=voigtFit1->GetParameter(1);
        std::cout<<"bin "<<k<<" avPU "<<binAveragePU[k]<<" Z-mass "<<zPeaks[k]<<std::endl;
    }
    
    //for make stupid TGraph work
    TH2D* dummy = new TH2D("dummy","Z-peak location Vs. Average PU;nVert+1;M_{Z}",1,0,25,1,90,95);    
    dummy->Draw();
    //dummy->GetXaxis()->SetRangeUser(0,25);
    //dummy->GetYaxis()->SetRangeUser(90,95);
    puFit->SetLineColor(kBlue);
    TGraph* zPeakByPu = new TGraph(4,&binAveragePU[0],&zPeaks[0]);
    zPeakByPu->SetMarkerStyle(20);
    zPeakByPu->SetMarkerSize(1.3);
    zPeakByPu->SetMarkerColor(kBlue+2);
    zPeakByPu->Fit(puFit,"QM","",0,25);
    zPeakByPu->GetXaxis()->SetRangeUser(0,25);
    zPeakByPu->GetYaxis()->SetRangeUser(90,95);
    zPeakByPu->Draw("SAME P");
    c1->Print("Calibration/PU/ZmassVsAveragePU.png");    
    
    Slopes1D->SetMarkerStyle(20);
    Slopes1D->Draw("P");
    c1->Print("Calibration/PU/Slopes1D.png");
    Ratios1D->SetMarkerStyle(20);
    Ratios1D->Draw("P");
    c1->Print("Calibration/PU/Ratios1D.png");
    slopesVsEtaN->SetMarkerStyle(20);
    slopesVsEtaN->SetMarkerStyle(20);
    slopesVsEtaN->SetMarkerSize(2);
    slopesVsEtaN->SetMarkerColor(kRed+2);
    slopesVsEtaN->GetYaxis()->SetRangeUser(-0.01,0);
    slopesVsEtaN->Draw("P");
    c1->Print("Calibration/PU/slopesVsEtaN.png");
    slopesVsEtaP->SetMarkerStyle(20);
    slopesVsEtaP->SetMarkerStyle(20);
    slopesVsEtaP->SetMarkerSize(2);
    slopesVsEtaP->SetMarkerColor(kRed+2);
    slopesVsEtaP->GetYaxis()->SetRangeUser(-0.01,0);
    slopesVsEtaP->Draw("P");
    c1->Print("Calibration/PU/slopesVsEtaP.png");
    ratiosVsEtaN->SetMarkerStyle(20);
    ratiosVsEtaN->SetMarkerStyle(20);
    ratiosVsEtaN->SetMarkerSize(2);
    ratiosVsEtaN->SetMarkerColor(kGreen+2);
    ratiosVsEtaN->GetYaxis()->SetRangeUser(1,1.2);
    ratiosVsEtaN->Draw("P");
    c1->Print("Calibration/PU/ratiosVsEtaN.png");
    ratiosVsEtaP->SetMarkerStyle(20);
    ratiosVsEtaP->SetMarkerSize(2);
    ratiosVsEtaP->SetMarkerColor(kGreen+2);
    ratiosVsEtaP->GetYaxis()->SetRangeUser(1,1.2);
    ratiosVsEtaP->Draw("P");
    c1->Print("Calibration/PU/ratiosVsEtaP.png");
    
    /*ZmassVsPU->SetMarkerStyle(20);
    ZmassVsPU->GetYaxis()->SetRangeUser(87,97);
    ZmassVsPU->Draw("E");
    puFit->SetParameter(0,92);
    puFit->SetParameter(1,0);
    //puFit->SetParameter(2,0);
    ZmassVsPU->Fit(puFit,"QLM","",5,25);
    c1->Print("Calibration/PU/ZmassVaPU.png");*/
    
    c1->Close();
    return 0;
}