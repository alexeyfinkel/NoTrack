//This is the script for validating the converged calibration constants.
//It should read in the cal. consts files from the "ODDS" and "EVENS" data subsets
//and apply them to the complimentary sets to see the effectiveness on restoring the mean and Z-mass uniformity

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
#include "TMath.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TColor.h"
#include "TProfile.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TColor.h"


#include "../src/ZEffTree.h"

Double_t voigt(Double_t *x, Double_t *par)
{
    return par[0]*(TMath::Voigt((x[0]-par[1]),par[2],par[3]));
}

int calConstValidation()
{
	gROOT->SetStyle("Plain");	
	gStyle->SetErrorX(0);
	gStyle->SetOptStat(0);
	
	TCanvas* c1 = new TCanvas("C1", "c1", 1200, 900);
	c1->cd();
	//c1->SetLeftMargin(0.08);
	//c1->SetRightMargin(0.12);
	
    TLegend* l1 = new TLegend(0.68,0.65,0.97,0.97);
	l1->SetFillColor(10);
    l1->SetTextSize(0.025);
    l1->SetBorderSize(0);
    //maps:
	std::map< pair<int,int>, TH1D* > ratioMapBeforeP, ratioMapBeforeN, zMassMapBeforeN, zMassMapBeforeP, ratioMapAfterP, ratioMapAfterN, zMassMapAfterN, zMassMapAfterP; 
	char title[256],name[128];// filename[128];
    for(int i=25;i<76;i++)
	{
		for(int j=25;j<76;j++)
		{
			//if( (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5))  < 132.25) || (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5)) > 342.25) ) continue;
			sprintf(name,"RatioIx%diy%dzP_Before",i,j);
			sprintf(title,"E_{exp}/E_{obs}, NT+ i_{x} = %d, i_{y} = %d;E_{exp}/E_{obs};Events",i,j);
			ratioMapBeforeP[std::make_pair(i,j)]= new TH1D(name,title,12,0.4,1.6);
			sprintf(name,"RatioIx%diy%dzN_Before",i,j);
			sprintf(title,"E_{exp}/E_{obs}, NT- i_{x} = %d, i_{y} = %d;E_{exp}/E_{obs};Events",i,j);			
			ratioMapBeforeN[std::make_pair(i,j)]= new TH1D(name,title,12,0.4,1.6);	
			sprintf(name,"zMassIx%diy%dzN_Before",i,j);
			sprintf(title,"M_{ee}, NT- i_{x} = %d, i_{y} = %d;M_{ee} (GeV);Events/2GeV",i,j);	
			zMassMapBeforeN[std::make_pair(i,j)]= new TH1D(name,title,30,60,120);
			sprintf(name,"zMassIx%diy%dzP_Before",i,j);
			sprintf(title,"M_{ee}, NT+ i_{x} = %d, i_{y} = %d;M_{ee} (GeV);Events/2GeV",i,j);	
			zMassMapBeforeP[std::make_pair(i,j)]= new TH1D(name,title,30,60,120);
            
            sprintf(name,"RatioIx%diy%dzP_After",i,j);
			sprintf(title,"E_{exp}/E_{obs}, NT+ i_{x} = %d, i_{y} = %d;E_{exp}/E_{obs};Events",i,j);
			ratioMapAfterP[std::make_pair(i,j)]= new TH1D(name,title,12,0.4,1.6);
			sprintf(name,"RatioIx%diy%dzN_After",i,j);
			sprintf(title,"E_{exp}/E_{obs}, NT- i_{x} = %d, i_{y} = %d;E_{exp}/E_{obs};Events",i,j);			
			ratioMapAfterN[std::make_pair(i,j)]= new TH1D(name,title,12,0.4,1.6);	
			sprintf(name,"zMassIx%diy%dzN_After",i,j);
			sprintf(title,"M_{ee}, NT- i_{x} = %d, i_{y} = %d;M_{ee} (GeV);Events/2GeV",i,j);	
			zMassMapAfterN[std::make_pair(i,j)]= new TH1D(name,title,30,60,120);
			sprintf(name,"zMassIx%diy%dzP_After",i,j);
			sprintf(title,"M_{ee}, NT+ i_{x} = %d, i_{y} = %d;M_{ee} (GeV);Events/2GeV",i,j);	
			zMassMapAfterP[std::make_pair(i,j)]= new TH1D(name,title,30,60,120);

		}
	}
    //odds 1-D hists:
    TH1D *ODDSmeansDistBeforeP = new TH1D("ODDSmeansDistBeforeP","Odd Evts, #mu_{E_{exp}/E_{obs}}, NT+, Before;#mu_{E_{exp}/E_{obs}};Events",40,0.8,1.2);
    TH1D *ODDSmeansDistBeforeN = new TH1D("ODDSmeansDistBeforeN","Odd Evts, #mu_{E_{exp}/E_{obs}}, NT-, Before;#mu_{E_{exp}/E_{obs}};Events",40,0.8,1.2);
    TH1D *ODDSmeansDistAfterP = new TH1D("ODDSmeansDistAfterP","Odd Evts, #mu_{E_{exp}/E_{obs}}, NT+, After;#mu_{E_{exp}/E_{obs}};Events",40,0.8,1.2);
    TH1D *ODDSmeansDistAfterN = new TH1D("ODDSmeansDistAfterN","Odd Evts, #mu_{E_{exp}/E_{obs}}, NT-, After;#mu_{E_{exp}/E_{obs}};Events",40,0.8,1.2);
    //TH1D *ODDSzMassDistBeforeP = new TH1D("ODDSzMassDistBeforeP","Odd Evts, M_{ee}, NT+, Before;M_{ee}",20,70,110);
    //TH1D *ODDSzMassDistAfterP = new TH1D("ODDSzMassDistAfterP","Odd Evts, M_{ee}, NT+, After;M_{ee};Events",20,70,110);
    TH1D *ODDSzMassDistAfterAll = new TH1D("ODDSzMassDistAfterAll","Odd Evts, M_{ee}, After;M_{ee}",40,70,110);
    TH1D *ODDSzMassDistBeforeAll = new TH1D("ODDSzMassDistBeforeAll","Odd Evts, M_{ee}, Before;M_{ee}",40,70,110);
    //and evens:
    TH1D *EVENSmeansDistBeforeP = new TH1D("EVENSmeansDistBeforeP","Odd Evts, #mu_{E_{exp}/E_{obs}}, NT+, Before;#mu_{E_{exp}/E_{obs}};Events",40,0.8,1.2);
    TH1D *EVENSmeansDistBeforeN = new TH1D("EVENSmeansDistBeforeN","Odd Evts, #mu_{E_{exp}/E_{obs}}, NT-, Before;#mu_{E_{exp}/E_{obs}};Events",40,0.8,1.2);
    TH1D *EVENSmeansDistAfterP = new TH1D("EVENSmeansDistAfterP","Odd Evts, #mu_{E_{exp}/E_{obs}}, NT+, After;#mu_{E_{exp}/E_{obs}};Events",40,0.8,1.2);
    TH1D *EVENSmeansDistAfterN = new TH1D("EVENSmeansDistAfterN","Odd Evts, #mu_{E_{exp}/E_{obs}}, NT-, After;#mu_{E_{exp}/E_{obs}};Events",40,0.8,1.2);
    //TH1D *EVENSzMassDistBeforeP = new TH1D("EVENSzMassDistBeforeP","Odd Evts, M_{ee}, NT+, Before;M_{ee}",20,70,110);
    //TH1D *EVENSzMassDistAfterP = new TH1D("EVENSzMassDistAfterP","Odd Evts, M_{ee}, NT+, After;M_{ee};Events",20,70,110);
    TH1D *EVENSzMassDistAfterAll = new TH1D("EVENSzMassDistAfterAll","Odd Evts, M_{ee}, After;M_{ee}",40,70,110);
    TH1D *EVENSzMassDistBeforeAll = new TH1D("EVENSzMassDistBeforeAll","Odd Evts, M_{ee}, Before;M_{ee}",40,70,110);
    //ratios:
    TH1D *ratiosP = new TH1D("ratiosP","Odd/Even Cal. Consts NT+;C_{odd}/C_{even};Events",20,0.9,1.1);
    TH1D *ratiosN = new TH1D("ratiosN","Odd/Even Cal. Consts NT-;C_{odd}/C_{even};Events",20,0.9,1.1);
    TH1D *ratiosAll = new TH1D("ratiosAll","Odd/Even Cal. Consts All NT;C_{odd}/C_{even};Events",20,0.9,1.1);
    
    TH1D *constsDiff = new TH1D("constsDiff","C_{odd}-C_{even};C_{odd}-C_{even};Events",20,-0.2,0.2);
    TH1D *constsDiffPull = new TH1D("constsDiffPull","Odd/Even Pull Distribution;P;Events",20,-5,5);
    
    //2-D ratios: (coming soon))
    TH2D *calConstRatiosN = new TH2D("calConstRatiosN","C_{odd}/C_{even}, NT-;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *calConstRatiosP = new TH2D("calConstRatiosP","C_{odd}/C_{even}, NT+;i_{x};i_{y}",40,30,70,40,30,70);
    
    //odds 2D hists:
    TH2D *ODDSmeansBeforeN = new TH2D("ODDSmeansBeforeN","Odd Evts. Mean Before, NT-;i_{x};i_{y}",40,30,70,40,30,70);
	TH2D *ODDSmeansBeforeP = new TH2D("ODDSmeansBeforeP","Odd Evts. Mean Before, NT+;i_{x};i_{y}",40,30,70,40,30,70);  
    TH2D *ODDSmeansAfterN = new TH2D("ODDSmeansAfterN","Odd Evts. Mean After, NT-;i_{x};i_{y}",40,30,70,40,30,70);
	TH2D *ODDSmeansAfterP = new TH2D("ODDSmeansAfterP","Odd Evts. Mean After, NT+;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *ODDSzMassBeforeN = new TH2D("ODDSzMassBeforeN","Odd Evts. M_{ee} Before, NT-;i_{x};i_{y}",40,30,70,40,30,70);
	TH2D *ODDSzMassBeforeP = new TH2D("ODDSzMassBeforeP","Odd Evts. M_{ee} Before, NT+;i_{x};i_{y}",40,30,70,40,30,70);  
    TH2D *ODDSzMassAfterN = new TH2D("ODDSzMassAfterN","Odd Evts. M_{ee} After, NT-;i_{x};i_{y}",40,30,70,40,30,70);
	TH2D *ODDSzMassAfterP = new TH2D("ODDSzMassAfterP","Odd Evts. M_{ee} After, NT+;i_{x};i_{y}",40,30,70,40,30,70);   
    //and evens:
    TH2D *EVENSmeansBeforeN = new TH2D("EVENSmeansBeforeN","Even Evts. Mean Before, NT-;i_{x};i_{y}",40,30,70,40,30,70);
	TH2D *EVENSmeansBeforeP = new TH2D("EVENSmeansBeforeP","Even Evts. Mean Before, NT+;i_{x};i_{y}",40,30,70,40,30,70);  
    TH2D *EVENSmeansAfterN = new TH2D("EVENSmeansAfterN","Even Evts. Mean After, NT-;i_{x};i_{y}",40,30,70,40,30,70);
	TH2D *EVENSmeansAfterP = new TH2D("EVENSmeansAfterP","Even Evts. Mean After, NT+;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *EVENSzMassBeforeN = new TH2D("EVENSzMassBeforeN","Even Evts. M_{ee} Before, NT-;i_{x};i_{y}",40,30,70,40,30,70);
	TH2D *EVENSzMassBeforeP = new TH2D("EVENSzMassBeforeP","Even Evts. M_{ee} Before, NT+;i_{x};i_{y}",40,30,70,40,30,70);  
    TH2D *EVENSzMassAfterN = new TH2D("EVENSzMassAfterN","Even Evts. M_{ee} After, NT-;i_{x};i_{y}",40,30,70,40,30,70);
	TH2D *EVENSzMassAfterP = new TH2D("EVENSzMassAfterP","Even Evts. M_{ee} After, NT+;i_{x};i_{y}",40,30,70,40,30,70);
    
    //cal. consts histograms:
    TH2D* ODDScalConstsN = new TH2D("ODDScalConstsN","Calib. Consts for Odd Evts. (found from Even), NT-;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D* ODDScalConstsP = new TH2D("ODDScalConstsP","Calib. Consts for Odd Evts. (found from Even), NT+;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D* EVENScalConstsN = new TH2D("EVENScalConstsN","Calib. Consts for Even Evts. (found from Odd), NT-;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D* EVENScalConstsP = new TH2D("EVENScalConstsP","Calib. Consts for Even Evts. (found from Odd), NT+;i_{x};i_{y}",40,30,70,40,30,70);
    //and their uncertainties (for pull)
    TH2D* ODDSsigmaN = new TH2D("ODDSsigmaN","Cal. Const. Uncertainty for Odd Evts. (found from ODDS), NT-;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D* ODDSsigmaP = new TH2D("ODDSsigmaP","Cal. Const. Uncertainty for Odd Evts. (found from ODDS), NT+;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D* EVENSsigmaN = new TH2D("EVENSsigmaN","Cal. Const. Uncertainty for Even Evts. (found from EVENS), NT-;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D* EVENSsigmaP = new TH2D("EVENSsigmaP","Cal. Const. Uncertainty for Even Evts. (found from EVENS), NT+;i_{x};i_{y}",40,30,70,40,30,70);
    
    char line[1000];
    int cix, ciy, iz;
    float calConst, sigma;
    ifstream oddsFile;
    //read in ODDS cal. consts file (read into "EVENS" hists!):
	oddsFile.open("Calibration/ODDS/CalConsts_Data_ODDS.txt");
	if(!oddsFile.is_open())
	{
		std::cout<<"Failed to open cal. constants file. Existing"<<std::endl;
		return 1;
	}
    do
    {
        oddsFile.getline(line, 1000)>>cix>>ciy>>iz>>calConst>>sigma;
        if(iz==-1) 
        {
            EVENScalConstsN->SetBinContent(cix-30,ciy-30,calConst);
            EVENSsigmaN->SetBinContent(cix-30,ciy-30,sigma);
        }
        else if(iz==1)
        {
            EVENScalConstsP->SetBinContent(cix-30,ciy-30,calConst);
            EVENSsigmaP->SetBinContent(cix-30,ciy-30,sigma);
        }
        else 
        {
            std::cout<<"Error reading in a cal. const. Exiting"<<std::endl;
            return 1;
        }
        
    } while(!oddsFile.eof());
    
    //now read in EVENS (use for ODDS!):
    ifstream evensFile;
	evensFile.open("Calibration/EVENS/CalConsts_Data_EVENS.txt");
	if(!evensFile.is_open())
	{
		std::cout<<"Failed to open cal. constants file. Existing"<<std::endl;
		return 1;
	}
    do
    {
        evensFile.getline(line, 1000)>>cix>>ciy>>iz>>calConst>>sigma;
        if(iz==-1) 
        {
            ODDScalConstsN->SetBinContent(cix-30,ciy-30,calConst);
            ODDSsigmaN->SetBinContent(cix-30,ciy-30,sigma);
        }
        else if(iz==1)
        {
            ODDScalConstsP->SetBinContent(cix-30,ciy-30,calConst);
            ODDSsigmaP->SetBinContent(cix-30,ciy-30,sigma);
        }
        else 
        {
            std::cout<<"Error reading in a cal. const. Exiting"<<std::endl;
            return 1;
        }
    } while(!evensFile.eof());
    
    //grab a data ntuple
    TFile* f1 = new TFile("/local/cms/user/finkel/NoTrack/Ntuple/DE_ReReco_2012Full_WithClusters_001_ODDS.root");	
	if(f1 == NULL)
	{
		std::cout<<"Failed to open Odds file. Exiting."<<std::endl;
		return 1;
	}
	TFile* f2 = new TFile("/local/cms/user/finkel/NoTrack/Ntuple/DE_ReReco_2012Full_WithClusters_001_EVENS.root");	
	if(f2 == NULL)
	{
		std::cout<<"Failed to open Evens file. Exiting."<<std::endl;
		return 1;
	}
    //make data ntuple:
    ZEffTree* ze1 = new ZEffTree(*f1,false);
	ZEffTree* ze2 = new ZEffTree(*f2,false);
    
    //kinematic vars:
    int ix,iy;
    int hix, hiy;
    const double Mz = 91.1876; //nominal Z mass
    double expectedPt;
    double observedPt;
    double correction;

    TLorentzVector elec1, elec2, theZ;
    float pt, eta, phi, E;
    
    //loop over ODDS dataset-----------------------------------------------------------------------------------------------------
    for(int event=0; event<ze1->Entries(); event++ ) 
    {
        ix=ze1->reco.ix[1];
        iy=ze1->reco.iy[1];
        if( (fabs(ze1->reco.eta[1])>2.5) && (fabs(ze1->reco.eta[1])<3.0) 
             &&(ix>29) &&(ix<71) &&(iy>29) &&(iy<71)
             && ze1->reco.isSelected(1,"NTLooseElectronId-EtaDet") //using only events that pass selection now!
          ) 
        {	
            if( (((ix-50.5)*(ix-50.5)+(iy-50.5)*(iy-50.5))  < 132.25) || (((ix-50.5)*(ix-50.5)+(iy-50.5)*(iy-50.5)) > 342.25) )
            {
                ze1->GetNextEvent();
                continue;
            }
            correction=0;
            for(unsigned int k=0; k<(ze1->ixs->size()); k++ )
            {	
                hix = ze1->ixs->at(k);
                hiy = ze1->iys->at(k);
                //if( (((hix-50.5)*(hix-50.5)+(hiy-50.5)*(hiy-50.5))  < 132.25) || (((hix-50.5)*(hix-50.5)+(hiy-50.5)*(hiy-50.5)) > 342.25) ) correction+=(ze1->hitEnergyFractions->at(k));
                if(ze1->reco.eta[1]>0) correction+= ( (ODDScalConstsP->GetBinContent(hix-30,hiy-30,calConst) )*(ze1->hitEnergyFractions->at(k)) );//positive endcap 
                else correction+=( (ODDScalConstsN->GetBinContent(hix-30,hiy-30,calConst))*(ze1->hitEnergyFractions->at(k)) );//negative endcap
            }					
            //recalculate observed Pt using correction:
            observedPt = ze1->reco.pt[1]*correction;

            //recalculate Z mass using found correction
            pt = ze1->reco.pt[0];    	//tracked electron
            eta = ze1->reco.eta[0];
            phi = ze1->reco.phi[0];
            E = pt*cosh(eta);
            elec1.SetPtEtaPhiE(pt,eta,phi,E);				
            pt = ze1->reco.pt[1]; 		//now untracked
            pt *= correction;	
            eta = ze1->reco.eta[1];
            phi = ze1->reco.phi[1];    
            E = pt*cosh(eta);	
            elec2.SetPtEtaPhiE(pt,eta,phi,E);
            theZ = elec1+elec2;		//new Z vector
            ODDSzMassDistBeforeAll->Fill(ze1->reco.mz);
            ODDSzMassDistAfterAll->Fill(theZ.M());
            //"expected energy
            expectedPt = Mz*Mz / ( 2*ze1->reco.pt[0]*( cosh(ze1->reco.eta[1]-ze1->reco.eta[0]) - cos(ze1->reco.phi[1]-ze1->reco.phi[0]) ) );            
            if(ze1->reco.eta[1]>0)//positive endcap
            {				
                if( (ze1->reco.mz>75) && (ze1->reco.mz<110) )
                {
                    zMassMapBeforeP[std::make_pair(ix,iy)]->Fill(ze1->reco.mz);
                    zMassMapAfterP[std::make_pair(ix,iy)]->Fill(theZ.M());
                    for(unsigned int k=0; k<(ze1->ixs->size()); k++ )
                    {
                        hix = ze1->ixs->at(k);
                        hiy = ze1->iys->at(k);
                        ratioMapBeforeP[std::make_pair(hix,hiy)]->Fill(expectedPt/ze1->reco.pt[1],ze1->hitEnergyFractions->at(k));
                        ratioMapAfterP[std::make_pair(hix,hiy)]->Fill(expectedPt/observedPt,ze1->hitEnergyFractions->at(k));  //NOTE: Using Expected/Observed!

                    }
                }
            }

            if(ze1->reco.eta[1]<0)//negative endcap
            {				
                if( (ze1->reco.mz>75) && (ze1->reco.mz<110) )
                {
                    zMassMapBeforeN[std::make_pair(ix,iy)]->Fill(ze1->reco.mz);
                    zMassMapAfterN[std::make_pair(ix,iy)]->Fill(theZ.M());
                    for(unsigned int k=0; k<(ze1->ixs->size()); k++ )
                    {
                        hix = ze1->ixs->at(k);
                        hiy = ze1->iys->at(k);
                        ratioMapBeforeN[std::make_pair(hix,hiy)]->Fill(expectedPt/ze1->reco.pt[1],ze1->hitEnergyFractions->at(k));
                        ratioMapAfterN[std::make_pair(hix,hiy)]->Fill(expectedPt/observedPt,ze1->hitEnergyFractions->at(k));  //NOTE: Using Expected/Observed!

                    }
                }
            }
        }
        ze1->GetNextEvent();
    }
    //end loop over events in ODD ntuple-----------------------------------------------------------------------------------------
    //do fits;
    TF1 *ratioFit = new TF1("fit","gaus",0.4,1.6);
	TF1 *massFit = new TF1("fit","gaus",70,110);
    TF1 *pullFit = new TF1("fit","gaus",-5,5);
    TF1 *voigtFit1 = new TF1("fit1",voigt,80,105,4);
    TF1 *voigtFit2 = new TF1("fit2",voigt,80,105,4);
    
    for(int i=30;i<71;i++)
    {
        for(int j=30;j<71;j++)
        {
            if( (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5))  < 132.25) || (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5)) > 342.25) ) continue;
            if(ratioMapBeforeN[std::make_pair(i,j)]->GetEntries() > 30)
            {
                ratioMapBeforeN[std::make_pair(i,j)]->Fit(ratioFit,"QLMN","",0.4,1.6);
                ODDSmeansBeforeN->SetBinContent(i-30,j-30,ratioFit->GetParameter(1));
                ODDSmeansDistBeforeN->Fill(ratioFit->GetParameter(1));
                ratioMapAfterN[std::make_pair(i,j)]->Fit(ratioFit,"QLMN","",0.4,1.6);
                ODDSmeansAfterN->SetBinContent(i-30,j-30,ratioFit->GetParameter(1));
                ODDSmeansDistAfterN->Fill(ratioFit->GetParameter(1));
                zMassMapBeforeN[std::make_pair(i,j)]->Fit(massFit,"QLMN","",80,110);
                ODDSzMassBeforeN->SetBinContent(i-30,j-30,massFit->GetParameter(1));
                zMassMapAfterN[std::make_pair(i,j)]->Fit(massFit,"QLMN","",75,105);
                ODDSzMassAfterN->SetBinContent(i-30,j-30,massFit->GetParameter(1));
            }
            if(ratioMapBeforeP[std::make_pair(i,j)]->GetEntries() > 30)
            {
                ratioMapBeforeP[std::make_pair(i,j)]->Fit(ratioFit,"QLMN","",0.4,1.6);
                ODDSmeansBeforeP->SetBinContent(i-30,j-30,ratioFit->GetParameter(1));
                ODDSmeansDistBeforeP->Fill(ratioFit->GetParameter(1));
                ratioMapAfterP[std::make_pair(i,j)]->Fit(ratioFit,"QLMN","",0.4,1.6);
                ODDSmeansAfterP->SetBinContent(i-30,j-30,ratioFit->GetParameter(1));
                ODDSmeansDistAfterP->Fill(ratioFit->GetParameter(1));
                zMassMapBeforeP[std::make_pair(i,j)]->Fit(massFit,"QLMN","",80,110);
                ODDSzMassBeforeP->SetBinContent(i-30,j-30,massFit->GetParameter(1));
                //ODDSzMassDistBeforeA->Fill(massFit->GetParameter(1));
                zMassMapAfterP[std::make_pair(i,j)]->Fit(massFit,"QLMN","",75,105);
                ODDSzMassAfterP->SetBinContent(i-30,j-30,massFit->GetParameter(1));
                //ODDSzMassDistAfterP->Fill(massFit->GetParameter(1));
            }
        }
    }
    // and make 2-D plots:
    ODDSmeansBeforeN->GetZaxis()->SetRangeUser(0.75,1.25);
	ODDSmeansBeforeN->GetZaxis()->SetLabelSize(0.03);
	ODDSmeansBeforeN->Draw("colz");
	c1->Print("Calibration/ODDS/ODDSMeansBeforeN.png");
	c1->Clear();
    ODDSmeansBeforeP->GetZaxis()->SetRangeUser(0.75,1.25);
	ODDSmeansBeforeP->GetZaxis()->SetLabelSize(0.03);
	ODDSmeansBeforeP->Draw("colz");
	c1->Print("Calibration/ODDS/ODDSMeansBeforeP.png");
	c1->Clear();
    ODDSmeansAfterN->GetZaxis()->SetRangeUser(0.75,1.25);
	ODDSmeansAfterN->GetZaxis()->SetLabelSize(0.03);
	ODDSmeansAfterN->Draw("colz");
	c1->Print("Calibration/ODDS/ODDSMeansAfterN.png");
	c1->Clear();
    ODDSmeansAfterP->GetZaxis()->SetRangeUser(0.75,1.25);
	ODDSmeansAfterP->GetZaxis()->SetLabelSize(0.03);
	ODDSmeansAfterP->Draw("colz");
	c1->Print("Calibration/ODDS/ODDSMeansAfterP.png");
	c1->Clear();
    ODDSzMassBeforeN->GetZaxis()->SetRangeUser(82,102);
	ODDSzMassBeforeN->GetZaxis()->SetLabelSize(0.03);
	ODDSzMassBeforeN->Draw("colz");
	c1->Print("Calibration/ODDS/ODDSZMassBeforeN.png");
	c1->Clear();
    ODDSzMassBeforeP->GetZaxis()->SetRangeUser(82,102);
	ODDSzMassBeforeP->GetZaxis()->SetLabelSize(0.03);
	ODDSzMassBeforeP->Draw("colz");
	c1->Print("Calibration/ODDS/ODDSZMassBeforeP.png");
	c1->Clear();
    ODDSzMassAfterN->GetZaxis()->SetRangeUser(82,102);
	ODDSzMassAfterN->GetZaxis()->SetLabelSize(0.03);
	ODDSzMassAfterN->Draw("colz");
	c1->Print("Calibration/ODDS/ODDSZMassAfterN.png");
	c1->Clear();
    ODDSzMassAfterP->GetZaxis()->SetRangeUser(82,102);
	ODDSzMassAfterP->GetZaxis()->SetLabelSize(0.03);
	ODDSzMassAfterP->Draw("colz");
	c1->Print("Calibration/ODDS/ODDSZMassAfterP.png");
	c1->Clear();
    //draw 1-D hists
    ODDSmeansDistBeforeN->SetMarkerStyle(20);
	ODDSmeansDistBeforeN->Draw("P");
	c1->Print("Calibration/ODDS/ODDSMeansDistBeforeN.png");
	c1->Clear();
    ODDSmeansDistBeforeP->SetMarkerStyle(20);
	ODDSmeansDistBeforeP->Draw("P");
	c1->Print("Calibration/ODDS/ODDSMeansDistBeforeP.png");
	c1->Clear();
    ODDSmeansDistAfterN->SetMarkerStyle(20);
	ODDSmeansDistAfterN->Draw("P");
	c1->Print("Calibration/ODDS/ODDSMeansDistAfterN.png");
	c1->Clear();
    ODDSmeansDistAfterP->SetMarkerStyle(20);
	ODDSmeansDistAfterP->Draw("P");
	c1->Print("Calibration/ODDS/ODDSMeansDistAfterP.png");
	c1->Clear();
    //ODDSzMassDistAfterP->SetMarkerStyle(20);
	//ODDSzMassDistAfterP->Draw("P");
	//c1->Print("Calibration/ODDS/ODDSZMassDistAfterP.png");
	//c1->Clear();
        
    voigtFit1->SetParameter(0,1e5);
    voigtFit1->SetParameter(1,90);
    voigtFit1->SetParLimits(1,80,100);
    voigtFit1->SetParameter(2,5);
    voigtFit1->SetParLimits(2,1,20);
    voigtFit1->SetParameter(3,5);
    voigtFit1->SetParLimits(3,1,20);
    ODDSzMassDistBeforeAll->SetMarkerStyle(20);
    ODDSzMassDistBeforeAll->Fit(voigtFit1,"LMN","",80,105);
    //sprintf(title,"#mu = %f\n#sigma = %f\nlg = %f",(float)voigtFit1->GetParameter(1),(float)voigtFit1->GetParameter(2),(float)voigtFit1->GetParameter(3)/2);
    ODDSzMassDistBeforeAll->Draw("P");
    sprintf(title,"#mu = %f",(float)voigtFit1->GetParameter(1) );
    l1->AddEntry(ODDSzMassDistBeforeAll,title,"l");
    sprintf(title,"#sigma = %f",(float)voigtFit1->GetParameter(2) );
    l1->AddEntry(ODDSzMassDistBeforeAll,title,"");
    sprintf(title,"lg = %f",(float)voigtFit1->GetParameter(3)/2 );
    l1->AddEntry(ODDSzMassDistBeforeAll,title,"");
    l1->Draw();
	c1->Print("Calibration/ODDS/ODDSZMassDistBeforeAll.png");
	c1->Clear();
    l1->Clear();
    voigtFit1->SetParameter(1,90);
    voigtFit1->SetParLimits(1,80,100);
    voigtFit1->SetParameter(2,5);
    voigtFit1->SetParLimits(2,1,20);
    voigtFit1->SetParameter(3,5);
    voigtFit1->SetParLimits(3,1,20);
    ODDSzMassDistAfterAll->SetMarkerStyle(20);
    ODDSzMassDistAfterAll->Fit(voigtFit1,"LMN","",80,105);
    //sprintf(title,"#mu = %f\n#sigma = %f\nlg = %f",(float)voigtFit1->GetParameter(1),(float)voigtFit1->GetParameter(2),(float)voigtFit1->GetParameter(3)/2);
    ODDSzMassDistAfterAll->Draw("P");
    sprintf(title,"#mu = %f",(float)voigtFit1->GetParameter(1) );
    l1->AddEntry(ODDSzMassDistBeforeAll,title,"l");
    sprintf(title,"#sigma = %f",(float)voigtFit1->GetParameter(2) );
    l1->AddEntry(ODDSzMassDistBeforeAll,title,"");
    sprintf(title,"lg = %f",(float)voigtFit1->GetParameter(3)/2 );
    l1->AddEntry(ODDSzMassDistBeforeAll,title,"");
    l1->Draw();
	c1->Print("Calibration/ODDS/ODDSZMassDistAfterAll.png");
	c1->Clear();
    l1->Clear();
    
    //nice combination plot
    c1->SetMargin(0.12,0.02,0.09,0.025);
    
    ODDSzMassDistAfterAll->SetTitle(";M_{ee} (GeV);Events/GeV");
    ODDSzMassDistAfterAll->GetYaxis()->SetTitleOffset(1.4);
    ODDSzMassDistAfterAll->GetYaxis()->SetTitleSize(0.045);
    ODDSzMassDistAfterAll->SetMarkerStyle(21);
    ODDSzMassDistAfterAll->SetLineColor(kBlue);
    ODDSzMassDistAfterAll->SetMarkerColor(kBlue+3);
    ODDSzMassDistAfterAll->SetMarkerSize(1.3);
    voigtFit1->SetParameter(1,90);
    voigtFit1->SetParLimits(1,80,100);
    voigtFit1->SetParameter(2,5);
    voigtFit1->SetParLimits(2,1,20);
    voigtFit1->FixParameter(3,4.9904);
    voigtFit1->SetLineWidth(4);
    voigtFit1->SetLineColor(kBlue);
    ODDSzMassDistAfterAll->Fit(voigtFit1,"LMN","",80,105);
    sprintf(title,"After Cross-Calibration" );
    l1->AddEntry(voigtFit1,title,"l");
    sprintf(title,"#mu = %.3g, #sigma = %.3g",(float)voigtFit1->GetParameter(1),(float)voigtFit1->GetParameter(2) );
    l1->AddEntry(voigtFit1,title,"");
    
    ODDSzMassDistBeforeAll->SetTitle(";;");
    ODDSzMassDistBeforeAll->SetMarkerStyle(20);
    ODDSzMassDistBeforeAll->SetLineColor(kRed);
    ODDSzMassDistBeforeAll->SetMarkerColor(kRed+2);
    ODDSzMassDistBeforeAll->SetMarkerSize(1.3);
    voigtFit2->SetParameter(1,90);
    voigtFit2->SetParLimits(1,80,100);
    voigtFit2->SetParameter(2,5);
    voigtFit2->SetParLimits(2,1,20);
    voigtFit2->FixParameter(3,4.9904);
    voigtFit2->SetLineWidth(4);
    voigtFit2->SetLineColor(kRed);
    ODDSzMassDistBeforeAll->Fit(voigtFit2,"LMN","",80,105);
    sprintf(title,"Before Cross-Calibration" );
    l1->AddEntry(voigtFit2,title,"l");
    sprintf(title,"#mu = %.3g, #sigma = %.3g",(float)voigtFit2->GetParameter(1),(float)voigtFit2->GetParameter(2) );
    l1->AddEntry(voigtFit2,title,"");
    
    c1->Clear();//just in case    
    ODDSzMassDistAfterAll->DrawCopy("E");
    voigtFit1->Draw("SAME");
    ODDSzMassDistAfterAll->DrawCopy(" SAME E");
    //ODDSzMassDistBeforeAll->DrawCopy("SAME E");
    voigtFit2->Draw("SAME");
    ODDSzMassDistBeforeAll->DrawCopy("SAME E");
    l1->Draw();
    
    TLatex theTitle;
	theTitle.SetTextSize(0.03);
	theTitle.SetTextFont(42);
	theTitle.SetNDC(true);    
    theTitle.DrawLatex(0.16,0.9,"CMS preliminary");
    theTitle.DrawLatex(0.16,0.85,"|#eta_{1}| < 2.5,  2.5 < |#eta_{2}| < 3.0");
    theTitle.DrawLatex(0.16,0.8,"2012 -- #sqrt{s} = 8 TeV");
    theTitle.DrawLatex(0.16,0.75,"#intL = 19 fb^{-1}, split sample calibration");
    
    //c1->Update();
    c1->Print("Calibration/ODDSzMassCompare.eps");
	c1->Clear();
    l1->Clear();
    
    std::cout<<"ODDS done."<<std::endl;
    
    c1->SetMargin(0.08,0.04,0.1,0.07);
    
    //NOW loop over the EVENS ntuple---------------------------------------------------------------------------------------------
    for(int event=0; event<ze2->Entries(); event++ ) 
    {
        ix=ze2->reco.ix[1];
        iy=ze2->reco.iy[1];
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
            correction=0;
            for(unsigned int k=0; k<(ze2->ixs->size()); k++ )
            {	
                hix = ze2->ixs->at(k);
                hiy = ze2->iys->at(k);
                //if( (((hix-50.5)*(hix-50.5)+(hiy-50.5)*(hiy-50.5))  < 132.25) || (((hix-50.5)*(hix-50.5)+(hiy-50.5)*(hiy-50.5)) > 342.25) ) correction+=(ze2->hitEnergyFractions->at(k));
                if(ze2->reco.eta[1]>0) correction+= ( (EVENScalConstsP->GetBinContent(hix-30,hiy-30,calConst) )*(ze2->hitEnergyFractions->at(k)) );//positive endcap 
                else correction+=( (EVENScalConstsN->GetBinContent(hix-30,hiy-30,calConst))*(ze2->hitEnergyFractions->at(k)) );//negative endcap
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
            EVENSzMassDistBeforeAll->Fill(ze2->reco.mz);
            EVENSzMassDistAfterAll->Fill(theZ.M());
            //"expected energy
            expectedPt = Mz*Mz / ( 2*ze2->reco.pt[0]*( cosh(ze2->reco.eta[1]-ze2->reco.eta[0]) - cos(ze2->reco.phi[1]-ze2->reco.phi[0]) ) );            
            if(ze2->reco.eta[1]>0)//positive endcap
            {				
                if( (ze2->reco.mz>75) && (ze2->reco.mz<110) )
                {
                    zMassMapBeforeP[std::make_pair(ix,iy)]->Fill(ze2->reco.mz);
                    zMassMapAfterP[std::make_pair(ix,iy)]->Fill(theZ.M());
                    for(unsigned int k=0; k<(ze2->ixs->size()); k++ )
                    {
                        hix = ze2->ixs->at(k);
                        hiy = ze2->iys->at(k);
                        ratioMapBeforeP[std::make_pair(hix,hiy)]->Fill(expectedPt/ze2->reco.pt[1],ze2->hitEnergyFractions->at(k));
                        ratioMapAfterP[std::make_pair(hix,hiy)]->Fill(expectedPt/observedPt,ze2->hitEnergyFractions->at(k));  //NOTE: Using Expected/Observed!

                    }
                }
            }

            if(ze2->reco.eta[1]<0)//negative endcap
            {				
                if( (ze2->reco.mz>75) && (ze2->reco.mz<110) )
                {
                    zMassMapBeforeN[std::make_pair(hix,hiy)]->Fill(ze2->reco.mz);
                    zMassMapAfterN[std::make_pair(hix,hiy)]->Fill(theZ.M());
                    for(unsigned int k=0; k<(ze2->ixs->size()); k++ )
                    {
                        hix = ze2->ixs->at(k);
                        hiy = ze2->iys->at(k);
                        ratioMapBeforeN[std::make_pair(hix,hiy)]->Fill(expectedPt/ze2->reco.pt[1],ze2->hitEnergyFractions->at(k));
                        ratioMapAfterN[std::make_pair(hix,hiy)]->Fill(expectedPt/observedPt,ze2->hitEnergyFractions->at(k));  //NOTE: Using Expected/Observed!

                    }
                }
            }
        }
        ze2->GetNextEvent();
    }
    //end loop over events in EVEN ntuple-----------------------------------------------------------------------------------------
    
    for(int i=30;i<71;i++)
    {
        for(int j=30;j<71;j++)
        {
            if( (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5))  < 132.25) || (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5)) > 342.25) ) continue;
            if(ratioMapBeforeN[std::make_pair(i,j)]->GetEntries() > 30)
            {
                ratioMapBeforeN[std::make_pair(i,j)]->Fit(ratioFit,"QLMN","",0.4,1.6);
                EVENSmeansBeforeN->SetBinContent(i-30,j-30,ratioFit->GetParameter(1));
                EVENSmeansDistBeforeN->Fill(ratioFit->GetParameter(1));
                ratioMapAfterN[std::make_pair(i,j)]->Fit(ratioFit,"QLMN","",0.4,1.6);
                EVENSmeansAfterN->SetBinContent(i-30,j-30,ratioFit->GetParameter(1));
                EVENSmeansDistAfterN->Fill(ratioFit->GetParameter(1));
                zMassMapBeforeN[std::make_pair(i,j)]->Fit(massFit,"QLMN","",80,110);
                EVENSzMassBeforeN->SetBinContent(i-30,j-30,massFit->GetParameter(1));
                //EVENSzMassDistBeforeN->Fill(massFit->GetParameter(1));
                zMassMapAfterN[std::make_pair(i,j)]->Fit(massFit,"QLMN","",75,105);
                EVENSzMassAfterN->SetBinContent(i-30,j-30,massFit->GetParameter(1));
                //EVENSzMassDistAfterN->Fill(massFit->GetParameter(1));
            }
            if(ratioMapBeforeP[std::make_pair(i,j)]->GetEntries() > 30)
            {
                ratioMapBeforeP[std::make_pair(i,j)]->Fit(ratioFit,"QLMN","",0.4,1.6);
                EVENSmeansBeforeP->SetBinContent(i-30,j-30,ratioFit->GetParameter(1));
                EVENSmeansDistBeforeP->Fill(ratioFit->GetParameter(1));
                ratioMapAfterP[std::make_pair(i,j)]->Fit(ratioFit,"QLMN","",0.4,1.6);
                EVENSmeansAfterP->SetBinContent(i-30,j-30,ratioFit->GetParameter(1));
                EVENSmeansDistAfterP->Fill(ratioFit->GetParameter(1));
                zMassMapBeforeP[std::make_pair(i,j)]->Fit(massFit,"QLMN","",80,110);
                EVENSzMassBeforeP->SetBinContent(i-30,j-30,massFit->GetParameter(1));
                //EVENSzMassDistBeforeP->Fill(massFit->GetParameter(1));
                zMassMapAfterP[std::make_pair(i,j)]->Fit(massFit,"QLMN","",75,105);
                EVENSzMassAfterP->SetBinContent(i-30,j-30,massFit->GetParameter(1));
                //EVENSzMassDistAfterP->Fill(massFit->GetParameter(1));
            }
        }
    }
    // and make 2-D plots:
    EVENSmeansBeforeN->GetZaxis()->SetRangeUser(0.75,1.25);
	EVENSmeansBeforeN->GetZaxis()->SetLabelSize(0.03);
	EVENSmeansBeforeN->Draw("colz");
	c1->Print("Calibration/EVENS/EVENSMeansBeforeN.png");
	c1->Clear();
    EVENSmeansBeforeP->GetZaxis()->SetRangeUser(0.75,1.25);
	EVENSmeansBeforeP->GetZaxis()->SetLabelSize(0.03);
	EVENSmeansBeforeP->Draw("colz");
	c1->Print("Calibration/EVENS/EVENSMeansBeforeP.png");
	c1->Clear();
    EVENSmeansAfterN->GetZaxis()->SetRangeUser(0.75,1.25);
	EVENSmeansAfterN->GetZaxis()->SetLabelSize(0.03);
	EVENSmeansAfterN->Draw("colz");
	c1->Print("Calibration/EVENS/EVENSMeansAfterN.png");
	c1->Clear();
    EVENSmeansAfterP->GetZaxis()->SetRangeUser(0.75,1.25);
	EVENSmeansAfterP->GetZaxis()->SetLabelSize(0.03);
	EVENSmeansAfterP->Draw("colz");
	c1->Print("Calibration/EVENS/EVENSMeansAfterP.png");
	c1->Clear();
    EVENSzMassBeforeN->GetZaxis()->SetRangeUser(82,102);
	EVENSzMassBeforeN->GetZaxis()->SetLabelSize(0.03);
	EVENSzMassBeforeN->Draw("colz");
	c1->Print("Calibration/EVENS/EVENSZMassBeforeN.png");
	c1->Clear();
    EVENSzMassBeforeP->GetZaxis()->SetRangeUser(82,102);
	EVENSzMassBeforeP->GetZaxis()->SetLabelSize(0.03);
	EVENSzMassBeforeP->Draw("colz");
	c1->Print("Calibration/EVENS/EVENSZMassBeforeP.png");
	c1->Clear();
    EVENSzMassAfterN->GetZaxis()->SetRangeUser(82,102);
	EVENSzMassAfterN->GetZaxis()->SetLabelSize(0.03);
	EVENSzMassAfterN->Draw("colz");
	c1->Print("Calibration/EVENS/EVENSZMassAfterN.png");
	c1->Clear();
    EVENSzMassAfterP->GetZaxis()->SetRangeUser(82,102);
	EVENSzMassAfterP->GetZaxis()->SetLabelSize(0.03);
	EVENSzMassAfterP->Draw("colz");
	c1->Print("Calibration/EVENS/EVENSZMassAfterP.png");
	c1->Clear();
    //draw 1-D hists
    EVENSmeansDistBeforeN->SetMarkerStyle(20);
	EVENSmeansDistBeforeN->Draw("P");
	c1->Print("Calibration/EVENS/EVENSMeansDistBeforeN.png");
	c1->Clear();
    EVENSmeansDistBeforeP->SetMarkerStyle(20);
	EVENSmeansDistBeforeP->Draw("P");
	c1->Print("Calibration/EVENS/EVENSMeansDistBeforeP.png");
	c1->Clear();
    EVENSmeansDistAfterN->SetMarkerStyle(20);
	EVENSmeansDistAfterN->Draw("P");
	c1->Print("Calibration/EVENS/EVENSMeansDistAfterN.png");
	c1->Clear();
    EVENSmeansDistAfterP->SetMarkerStyle(20);
	EVENSmeansDistAfterP->Draw("P");
	c1->Print("Calibration/EVENS/EVENSMeansDistAfterP.png");
	c1->Clear();
    //EVENSzMassDistBeforeN->SetMarkerStyle(20);
	//EVENSzMassDistBeforeN->Draw("P");
	//c1->Print("Calibration/EVENS/EVENSZMassDistBeforeN.png");
	//c1->Clear();
    //EVENSzMassDistBeforeP->SetMarkerStyle(20);
	//EVENSzMassDistBeforeP->Draw("P");
	//c1->Print("Calibration/EVENS/EVENSZMassDistBeforeP.png");
	//c1->Clear();
    //EVENSzMassDistAfterN->SetMarkerStyle(20);
	//EVENSzMassDistAfterN->Draw("P");
	//c1->Print("Calibration/EVENS/EVENSZMassDistAfterN.png");
	//c1->Clear();
    //EVENSzMassDistAfterP->SetMarkerStyle(20);
	//EVENSzMassDistAfterP->Draw("P");
	//c1->Print("Calibration/EVENS/EVENSZMassDistAfterP.png");
	//c1->Clear();
    
    voigtFit1->SetParameter(0,1e4);
    voigtFit1->SetParameter(1,90);
    voigtFit1->SetParLimits(1,80,100);
    voigtFit1->SetParameter(2,5);
    voigtFit1->SetParLimits(2,1,20);
    voigtFit1->SetParameter(3,5);
    voigtFit1->SetParLimits(3,1,20);
    EVENSzMassDistBeforeAll->SetMarkerStyle(20);
    EVENSzMassDistBeforeAll->Fit(voigtFit1,"LM","",80,105);
    //sprintf(title,"#mu = %f\n#sigma = %f\nlg = %f",(float)voigtFit1->GetParameter(1),(float)voigtFit1->GetParameter(2),(float)voigtFit1->GetParameter(3)/2);
    EVENSzMassDistBeforeAll->Draw("P");
    sprintf(title,"#mu = %f",(float)voigtFit1->GetParameter(1) );
    l1->AddEntry(EVENSzMassDistBeforeAll,title,"l");
    sprintf(title,"#sigma = %f",(float)voigtFit1->GetParameter(2) );
    l1->AddEntry(EVENSzMassDistBeforeAll,title,"");
    sprintf(title,"lg = %f",(float)voigtFit1->GetParameter(3)/2 );
    l1->AddEntry(EVENSzMassDistBeforeAll,title,"");
    l1->Draw();
	c1->Print("Calibration/EVENS/EVENSZMassDistBeforeAll.png");
	c1->Clear();
    l1->Clear();
    voigtFit1->SetParameter(1,90);
    voigtFit1->SetParLimits(1,80,100);
    voigtFit1->SetParameter(2,5);
    voigtFit1->SetParLimits(2,1,20);
    voigtFit1->SetParameter(3,5);
    voigtFit1->SetParLimits(3,1,20);
    EVENSzMassDistAfterAll->SetMarkerStyle(20);
    EVENSzMassDistAfterAll->Fit(voigtFit1,"LM","",80,105);
    //sprintf(title,"#mu = %f\n#sigma = %f\nlg = %f",(float)voigtFit1->GetParameter(1),(float)voigtFit1->GetParameter(2),(float)voigtFit1->GetParameter(3)/2);
    EVENSzMassDistAfterAll->Draw("P");
    sprintf(title,"#mu = %f",(float)voigtFit1->GetParameter(1) );
    l1->AddEntry(EVENSzMassDistBeforeAll,title,"l");
    sprintf(title,"#sigma = %f",(float)voigtFit1->GetParameter(2) );
    l1->AddEntry(EVENSzMassDistBeforeAll,title,"");
    sprintf(title,"lg = %f",(float)voigtFit1->GetParameter(3)/2 );
    l1->AddEntry(EVENSzMassDistBeforeAll,title,"");
    l1->Draw();
	c1->Print("Calibration/EVENS/EVENSZMassDistAfterAll.png");
	c1->Clear();
    l1->Clear();
        
    //produce ratio hists:-------------------------------------------------------------------------------------------------------
    //double rn, rp;
    double cOdd, cEven;
    double sigmaOdd, sigmaEven;
    double sigmaOld, sigmaAdd;
    for(int i=30;i<71;i++)
    {
        for(int j=30;j<71;j++)
        {
            if( (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5))  < 132.25) || (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5)) > 342.25) ) continue;
            //negative endcap
            cOdd = ODDScalConstsN->GetBinContent(i-30,j-30);
            cEven = EVENScalConstsN->GetBinContent(i-30,j-30);
            //rn = ODDScalConstsN->GetBinContent(i-30,j-30)/EVENScalConstsN->GetBinContent(i-30,j-30);
            //rp = ODDScalConstsP->GetBinContent(i-30,j-30)/EVENScalConstsP->GetBinContent(i-30,j-30);
            ratiosN->Fill( cOdd/cEven );
            ratiosAll->Fill( cOdd/cEven );
            calConstRatiosN->SetBinContent(i-30,j-30,cOdd/cEven);
            constsDiff->Fill( cOdd-cEven );
            sigmaOld = constsDiff->GetBinError(constsDiff->FindBin(cOdd-cEven));
            sigmaOdd = ODDSsigmaN->GetBinContent(i-30,j-30);
            sigmaEven = EVENSsigmaN->GetBinContent(i-30,j-30);
            sigmaAdd = sqrt(sigmaOdd*sigmaOdd+sigmaEven*sigmaEven);
            constsDiff->SetBinError(constsDiff->FindBin(cOdd-cEven), sigmaOld+sigmaAdd);
            constsDiffPull->Fill( (cOdd-cEven)/sigmaAdd );
            
            //positive endcap
            cOdd = ODDScalConstsP->GetBinContent(i-30,j-30);
            cEven = EVENScalConstsP->GetBinContent(i-30,j-30);
            ratiosP->Fill( cOdd/cEven );
            ratiosAll->Fill( cOdd/cEven );
            calConstRatiosP->SetBinContent(i-30,j-30,cOdd/cEven);
            constsDiff->Fill( cOdd-cEven );
            sigmaOld = constsDiff->GetBinError(constsDiff->FindBin(cOdd-cEven));
            sigmaOdd = ODDSsigmaP->GetBinContent(i-30,j-30);
            sigmaEven = EVENSsigmaP->GetBinContent(i-30,j-30);
            sigmaAdd = sqrt(sigmaOdd*sigmaOdd+sigmaEven*sigmaEven);
            constsDiff->SetBinError(constsDiff->FindBin(cOdd-cEven), sigmaOld+sigmaAdd);
            constsDiffPull->Fill( (cOdd-cEven)/sigmaAdd );
        }
    }
    //now make a pull distribution: 
    for(unsigned int bin=1;bin<(constsDiff->GetNbinsX());bin++)
    {
        constsDiff->SetBinError(bin, sqrt(constsDiff->GetBinError(bin) ) );
        //constsDiffPull->Fill( constsDiff->GetBinCenter(bin)/constsDiff->GetBinError(bin) );
        //(this one seems like I am doing it wrong...)
    }
    //draw them:
    TLegend* l2 = new TLegend(0.70,0.70,0.96,0.93);
    l2->SetFillColor(0);
    l2->SetTextSize(0.05);
    ratiosN->SetMarkerStyle(20);
    ratiosN->SetMarkerSize(1.2);
    ratioFit->SetParameter(2,0.01);
    ratiosN->Fit(ratioFit,"QLM","",0.95,1.05);
    sprintf(title,"#sigma = %.3g",(float)ratioFit->GetParameter(2) );
    l2->AddEntry(ratioFit,title,"");
	ratiosN->Draw("E");
    l2->Draw();
	c1->Print("Calibration/OddsEvensRatiosDistN.png");
	c1->Clear();
    l2->Clear();
    ratiosP->SetMarkerStyle(20);
    ratiosP->SetMarkerSize(1.2);
    ratioFit->SetParameter(2,0.01);
    ratiosP->Fit(ratioFit,"QLM","",0.95,1.05);
	ratiosP->Draw("E");
    sprintf(title,"#sigma = %.3g",(float)ratioFit->GetParameter(2) );
    l2->AddEntry(ratioFit,title,"");
    l2->Draw();
	c1->Print("Calibration/OddsEvensRatiosDistP.png");
	c1->Clear();
    l2->Clear();
    ratiosAll->SetMarkerStyle(20);
    ratiosAll->SetMarkerSize(1.4);
    ratioFit->SetParameter(2,0.01);
    ratiosAll->Fit(ratioFit,"QLM","",0.95,1.05);
    sprintf(title,"#sigma = %.3g",(float)ratioFit->GetParameter(2) );
    l2->AddEntry(ratioFit,title,"");
	ratiosAll->Draw("E");
    l2->Draw();
	c1->Print("Calibration/OddsEvensRatiosDistAll.png");
	c1->Clear();
    l2->Clear();
    constsDiff->SetMarkerStyle(20);
	constsDiff->Draw("P");
	c1->Print("Calibration/OddsEvensDifferenceAll.png");
	c1->Clear();
    
    constsDiffPull->SetMarkerStyle(20);
    constsDiffPull->Fit(pullFit,"QLM","",-5,5);
    sprintf(title,"#mu = %f, #sigma = %f",(float)pullFit->GetParameter(1),(float)pullFit->GetParameter(2));
    constsDiffPull->Draw("P");
    l1->AddEntry(constsDiffPull,title,"");
    l1->Draw();
	c1->Print("Calibration/OddsEvensDifferencePull.png");
	c1->Clear();
    
    calConstRatiosN->GetZaxis()->SetRangeUser(0.9,1.1);
	calConstRatiosN->GetZaxis()->SetLabelSize(0.03);
	calConstRatiosN->Draw("colz");
	c1->Print("Calibration/calConstRatiosN.png");
	c1->Clear();
    calConstRatiosP->GetZaxis()->SetRangeUser(0.9,1.1);
	calConstRatiosP->GetZaxis()->SetLabelSize(0.03);
	calConstRatiosP->Draw("colz");
	c1->Print("Calibration/calConstRatiosP.png");
	c1->Clear();
    
    c1->Close();
    
    return 0;
}