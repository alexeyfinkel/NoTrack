//This is the script for validating the converged calibration constants.
//It should read in the cal. consts files from the "" and "EVENS" data subsets
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


#include "../../../src/ZEffTree.h"

Double_t voigt(Double_t *x, Double_t *par)
{
    return (par[0]*(TMath::Voigt((x[0]-par[1]),par[2],par[3])) + par[4] + par[5]*(x[0]-80));
}

int calConstsHengneCrossCheck()
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
    TH1D *meansDistBeforeP = new TH1D("meansDistBeforeP","Odd Evts, #mu_{E_{exp}/E_{obs}}, NT+, Before;#mu_{E_{exp}/E_{obs}};Events",40,0.8,1.2);
    TH1D *meansDistBeforeN = new TH1D("meansDistBeforeN","Odd Evts, #mu_{E_{exp}/E_{obs}}, NT-, Before;#mu_{E_{exp}/E_{obs}};Events",40,0.8,1.2);
    TH1D *meansDistAfterP = new TH1D("meansDistAfterP","Odd Evts, #mu_{E_{exp}/E_{obs}}, NT+, After;#mu_{E_{exp}/E_{obs}};Events",40,0.8,1.2);
    TH1D *meansDistAfterN = new TH1D("meansDistAfterN","Odd Evts, #mu_{E_{exp}/E_{obs}}, NT-, After;#mu_{E_{exp}/E_{obs}};Events",40,0.8,1.2);
    TH1D *zMassDistAfterHengne = new TH1D("zMassDistAfterHengne","M_{ee}, After Hengne's;M_{ee}",40,70,110);
    TH1D *zMassDistAfterAlexey = new TH1D("zMassDistAfterAlexey","M_{ee}, After Alexey's;M_{ee}",40,70,110);
    TH1D *zMassDistBeforeAll = new TH1D("zMassDistBeforeAll","M_{ee}, Before and After Calibration;M_{ee}",40,70,110);
    //and evens:
       
    //odds 2D hists:
    TH2D *meansBeforeN = new TH2D("meansBeforeN","Odd Evts. Mean Before, NT-;i_{x};i_{y}",40,30,70,40,30,70);
	TH2D *meansBeforeP = new TH2D("meansBeforeP","Odd Evts. Mean Before, NT+;i_{x};i_{y}",40,30,70,40,30,70);  
    TH2D *meansAfterN = new TH2D("meansAfterN","Odd Evts. Mean After, NT-;i_{x};i_{y}",40,30,70,40,30,70);
	TH2D *meansAfterP = new TH2D("meansAfterP","Odd Evts. Mean After, NT+;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D *zMassBeforeN = new TH2D("zMassBeforeN","Odd Evts. M_{ee} Before, NT-;i_{x};i_{y}",40,30,70,40,30,70);
	TH2D *zMassBeforeP = new TH2D("zMassBeforeP","Odd Evts. M_{ee} Before, NT+;i_{x};i_{y}",40,30,70,40,30,70);  
    TH2D *zMassAfterN = new TH2D("zMassAfterN","Odd Evts. M_{ee} After, NT-;i_{x};i_{y}",40,30,70,40,30,70);
	TH2D *zMassAfterP = new TH2D("zMassAfterP","Odd Evts. M_{ee} After, NT+;i_{x};i_{y}",40,30,70,40,30,70);   
    
    //cal. consts histograms:
    TH2D* calConstsN = new TH2D("calConstsN","Hegne's Cal. Consts, NT-;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D* calConstsP = new TH2D("calConstsP","Hegne's Cal. Consts, NT+;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D* myConstsN = new TH2D("myConstsN","Alexey's Cal. Consts, NT-;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D* myConstsP = new TH2D("myConstsP","Alexey's Cal. Consts, NT+;i_{x};i_{y}",40,30,70,40,30,70);
    //and their uncertainties (for pull)
    TH2D* sigmaN = new TH2D("sigmaN","Cal. Const. Uncertainty for Odd Evts. (found from ), NT-;i_{x};i_{y}",40,30,70,40,30,70);
    TH2D* sigmaP = new TH2D("sigmaP","Cal. Const. Uncertainty for Odd Evts. (found from ), NT+;i_{x};i_{y}",40,30,70,40,30,70);
    
    char line[1000];
    int cix, ciy, iz, dummy;
    float calConst, sigma;
    ifstream pFile;
    //read in  positive consts file :
	pFile.open("calibTable_out_7thP_toPrint.dat");
	if(!pFile.is_open())
	{
		std::cout<<"Failed to open P-constants file. Existing"<<std::endl;
		return 1;
	}
    do
    {
        pFile.getline(line, 1000)>>dummy>>cix>>ciy>>calConst>>sigma;
        calConstsP->SetBinContent(cix-30,ciy-30,calConst);
        sigmaP->SetBinContent(cix-30,ciy-30,sigma);
        
    } while(!pFile.eof());
    
    //now read in Negative :
    ifstream nFile;
	nFile.open("calibTable_out_7thN_toPrint.dat");
	if(!nFile.is_open())
	{
		std::cout<<"Failed to open N-constants file. Existing"<<std::endl;
		return 1;
	}
    do
    {
        nFile.getline(line, 1000)>>dummy>>cix>>ciy>>calConst>>sigma;
        calConstsN->SetBinContent(cix-30,ciy-30,calConst);
        sigmaN->SetBinContent(cix-30,ciy-30,sigma);
        
    } while(!nFile.eof());
    
    //read in My cal consts
    ifstream myFile;
	myFile.open("/home/grad/finkel/work/NoTrack/CMSSW_5_3_8/src/NoTrack/MakeZEffTree/test/Calibration/CalConsts_Data.txt");
	if(!myFile.is_open())
	{
		std::cout<<"Failed to open My constants file. Existing"<<std::endl;
		return 1;
	}
    do
    {
        myFile.getline(line, 1000)>>dummy>>cix>>ciy>>iz>>calConst>>sigma;
        if(iz==-1)myConstsN->SetBinContent(cix-30,ciy-30,calConst);
        if(iz==1)myConstsP->SetBinContent(cix-30,ciy-30,calConst);
        
    } while(!myFile.eof());
    
    std::cout<<"Consts read in."<<std::endl;
    //grab a data ntuple
    TFile* f1 = new TFile("/local/cms/user/finkel/NoTrack/Ntuple/DE_ReReco_2012Full_WithClusters_001.root");	
	if(f1 == NULL)
	{
		std::cout<<"Failed to open Root file. Exiting."<<std::endl;
		return 1;
	}
	
    //make data ntuple:
    ZEffTree* ze1 = new ZEffTree(*f1,false);
    
    //kinematic vars:
    int ix,iy;
    int hix, hiy;
    const double Mz = 91.1876; //nominal Z mass
    double expectedPt;
    double observedPt;
    double correction1, correction2;

    TLorentzVector elec1, elec2, theZ;
    float pt, eta, phi, E;
    
    //loop over  dataset-----------------------------------------------------------------------------------------------------
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
            correction1=0;//Hengne's
            correction2=0;//mine
            for(unsigned int k=0; k<(ze1->ixs->size()); k++ )
            {	
                hix = ze1->ixs->at(k);
                hiy = ze1->iys->at(k);
                if( (((hix-50.5)*(hix-50.5)+(hiy-50.5)*(hiy-50.5))  < 132.25) || (((hix-50.5)*(hix-50.5)+(hiy-50.5)*(hiy-50.5)) > 342.25) ) 
                {
                    correction1+=(ze1->hitEnergyFractions->at(k));
                    correction2+=(ze1->hitEnergyFractions->at(k));
                }
                if(ze1->reco.eta[1]>0)
                {
                    correction1+= ( (calConstsP->GetBinContent(hix-30,hiy-30) )*(ze1->hitEnergyFractions->at(k)) );//positive endcap 
                    correction2+= ( (myConstsP->GetBinContent(hix-30,hiy-30) )*(ze1->hitEnergyFractions->at(k)) );
                }
                else 
                {
                    correction1+=( (calConstsN->GetBinContent(hix-30,hiy-30))*(ze1->hitEnergyFractions->at(k)) );//negative endcap
                    correction2+=( (myConstsN->GetBinContent(hix-30,hiy-30))*(ze1->hitEnergyFractions->at(k)) );
                }
            }					
            //recalculate observed Pt using correction:
            observedPt = ze1->reco.pt[1]*correction1;

            //recalculate Z mass using found correction
            pt = ze1->reco.pt[0];    	//tracked electron
            eta = ze1->reco.eta[0];
            phi = ze1->reco.phi[0];
            E = pt*cosh(eta);
            elec1.SetPtEtaPhiE(pt,eta,phi,E);				
            pt = ze1->reco.pt[1]; 		//now untracked
            pt *= correction1;	
            eta = ze1->reco.eta[1];
            phi = ze1->reco.phi[1];    
            E = pt*cosh(eta);	
            elec2.SetPtEtaPhiE(pt,eta,phi,E);
            theZ = elec1+elec2;		//new Z vector
            zMassDistBeforeAll->Fill(ze1->reco.mz);
            zMassDistAfterHengne->Fill(theZ.M());
            //my calibration
            zMassDistAfterAlexey->Fill(theZ.M());
            pt = ze1->reco.pt[1]; 		//now untracked
            pt *= correction2;
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
    //end loop over events in ntuple-----------------------------------------------------------------------------------------
    //do fits;
    TF1 *ratioFit = new TF1("fit","gaus",0.4,1.6);
	TF1 *massFit = new TF1("fit","gaus",70,110);
    //TF1 *pullFit = new TF1("fit","gaus",-5,5);
    TF1 *voigtFit1 = new TF1("fit1",voigt,80,105,6);
    TF1 *voigtFit2 = new TF1("fit2",voigt,80,105,6);
    TF1 *voigtFit3 = new TF1("fit2",voigt,80,105,6);
    
    for(int i=30;i<71;i++)
    {
        for(int j=30;j<71;j++)
        {
            if( (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5))  < 132.25) || (((i-50.5)*(i-50.5)+(j-50.5)*(j-50.5)) > 342.25) ) continue;
            if(ratioMapBeforeN[std::make_pair(i,j)]->GetEntries() > 30)
            {
                ratioMapBeforeN[std::make_pair(i,j)]->Fit(ratioFit,"QLMN","",0.4,1.6);
                meansBeforeN->SetBinContent(i-30,j-30,ratioFit->GetParameter(1));
                meansDistBeforeN->Fill(ratioFit->GetParameter(1));
                ratioMapAfterN[std::make_pair(i,j)]->Fit(ratioFit,"QLMN","",0.4,1.6);
                meansAfterN->SetBinContent(i-30,j-30,ratioFit->GetParameter(1));
                meansDistAfterN->Fill(ratioFit->GetParameter(1));
                zMassMapBeforeN[std::make_pair(i,j)]->Fit(massFit,"QLMN","",80,110);
                zMassBeforeN->SetBinContent(i-30,j-30,massFit->GetParameter(1));
                zMassMapAfterN[std::make_pair(i,j)]->Fit(massFit,"QLMN","",75,105);
                zMassAfterN->SetBinContent(i-30,j-30,massFit->GetParameter(1));
            }
            if(ratioMapBeforeP[std::make_pair(i,j)]->GetEntries() > 30)
            {
                ratioMapBeforeP[std::make_pair(i,j)]->Fit(ratioFit,"QLMN","",0.4,1.6);
                meansBeforeP->SetBinContent(i-30,j-30,ratioFit->GetParameter(1));
                meansDistBeforeP->Fill(ratioFit->GetParameter(1));
                ratioMapAfterP[std::make_pair(i,j)]->Fit(ratioFit,"QLMN","",0.4,1.6);
                meansAfterP->SetBinContent(i-30,j-30,ratioFit->GetParameter(1));
                meansDistAfterP->Fill(ratioFit->GetParameter(1));
                zMassMapBeforeP[std::make_pair(i,j)]->Fit(massFit,"QLMN","",80,110);
                zMassBeforeP->SetBinContent(i-30,j-30,massFit->GetParameter(1));
                //zMassDistBeforeA->Fill(massFit->GetParameter(1));
                zMassMapAfterP[std::make_pair(i,j)]->Fit(massFit,"QLMN","",75,105);
                zMassAfterP->SetBinContent(i-30,j-30,massFit->GetParameter(1));
                //zMassDistAfterP->Fill(massFit->GetParameter(1));
            }
        }
    }
    // and make 2-D plots:
    meansBeforeN->GetZaxis()->SetRangeUser(0.75,1.25);
	meansBeforeN->GetZaxis()->SetLabelSize(0.03);
	meansBeforeN->Draw("colz");
	c1->Print("MeansBeforeN.png");
	c1->Clear();
    meansBeforeP->GetZaxis()->SetRangeUser(0.75,1.25);
	meansBeforeP->GetZaxis()->SetLabelSize(0.03);
	meansBeforeP->Draw("colz");
	c1->Print("MeansBeforeP.png");
	c1->Clear();
    calConstsN->GetZaxis()->SetRangeUser(0.75,1.25);
	calConstsN->GetZaxis()->SetLabelSize(0.03);
	calConstsN->Draw("colz");
	c1->Print("CalConstsN.png");
    calConstsP->GetZaxis()->SetRangeUser(0.75,1.25);
	calConstsP->GetZaxis()->SetLabelSize(0.03);
	calConstsP->Draw("colz");
	c1->Print("CalConstsNP.png");
	c1->Clear();
    meansAfterN->GetZaxis()->SetRangeUser(0.75,1.25);
	meansAfterN->GetZaxis()->SetLabelSize(0.03);
	meansAfterN->Draw("colz");
	c1->Print("MeansAfterN.png");
	c1->Clear();
    meansAfterP->GetZaxis()->SetRangeUser(0.75,1.25);
	meansAfterP->GetZaxis()->SetLabelSize(0.03);
	meansAfterP->Draw("colz");
	c1->Print("MeansAfterP.png");
	c1->Clear();
    zMassBeforeN->GetZaxis()->SetRangeUser(82,102);
	zMassBeforeN->GetZaxis()->SetLabelSize(0.03);
	zMassBeforeN->Draw("colz");
	c1->Print("ZMassBeforeN.png");
	c1->Clear();
    zMassBeforeP->GetZaxis()->SetRangeUser(82,102);
	zMassBeforeP->GetZaxis()->SetLabelSize(0.03);
	zMassBeforeP->Draw("colz");
	c1->Print("ZMassBeforeP.png");
	c1->Clear();
    zMassAfterN->GetZaxis()->SetRangeUser(82,102);
	zMassAfterN->GetZaxis()->SetLabelSize(0.03);
	zMassAfterN->Draw("colz");
	c1->Print("ZMassAfterN.png");
	c1->Clear();
    zMassAfterP->GetZaxis()->SetRangeUser(82,102);
	zMassAfterP->GetZaxis()->SetLabelSize(0.03);
	zMassAfterP->Draw("colz");
	c1->Print("ZMassAfterP.png");
	c1->Clear();
    //draw 1-D hists
    meansDistBeforeN->SetMarkerStyle(20);
	meansDistBeforeN->Draw("P");
	c1->Print("MeansDistBeforeN.png");
	c1->Clear();
    meansDistBeforeP->SetMarkerStyle(20);
	meansDistBeforeP->Draw("P");
	c1->Print("MeansDistBeforeP.png");
	c1->Clear();
    meansDistAfterN->SetMarkerStyle(20);
	meansDistAfterN->Draw("P");
	c1->Print("MeansDistAfterN.png");
	c1->Clear();
    meansDistAfterP->SetMarkerStyle(20);
	meansDistAfterP->Draw("P");
	c1->Print("MeansDistAfterP.png");
	c1->Clear();
    //zMassDistAfterP->SetMarkerStyle(20);
	//zMassDistAfterP->Draw("P");
	//c1->Print("ZMassDistAfterP.png");
	//c1->Clear();
        
    
    voigtFit1->SetParameter(1,90);
    voigtFit1->SetParLimits(1,90,100);
    voigtFit1->SetParameter(2,5);
    voigtFit1->SetParLimits(2,1,20);
    voigtFit1->FixParameter(3,4.9904);
    voigtFit1->SetParameter(4,1000);
    voigtFit1->SetParameter(5,-10);
    voigtFit1->SetLineWidth(4);
    zMassDistBeforeAll->SetMarkerStyle(20);
    zMassDistBeforeAll->Fit(voigtFit1,"LMN","",80,105);
    //sprintf(title,"#mu = %f\n#sigma = %f\nlg = %f",(float)voigtFit1->GetParameter(1),(float)voigtFit1->GetParameter(2),(float)voigtFit1->GetParameter(3)/2);
    zMassDistBeforeAll->Draw("P");
    sprintf(title,"#mu = %f",(float)voigtFit1->GetParameter(1) );
    l1->AddEntry(zMassDistBeforeAll,title,"l");
    sprintf(title,"#sigma = %f",(float)voigtFit1->GetParameter(2) );
    l1->AddEntry(zMassDistBeforeAll,title,"");
    sprintf(title,"lg = %f",(float)voigtFit1->GetParameter(3)/2 );
    l1->AddEntry(zMassDistBeforeAll,title,"");
    l1->Draw();
	c1->Print("ZMassDistBeforeAll.png");
	c1->Clear();
    l1->Clear();
    
    
    voigtFit1->SetParameter(1,90);
    voigtFit1->SetParLimits(1,90,100);
    voigtFit1->SetParameter(2,5);
    voigtFit1->SetParLimits(2,1,20);
    voigtFit1->FixParameter(3,4.9904);
    voigtFit1->SetParameter(4,1000);
    voigtFit1->SetParameter(5,-10);
    voigtFit1->SetLineWidth(4);
    zMassDistAfterHengne->SetMarkerStyle(20);
    zMassDistAfterHengne->Fit(voigtFit1,"LMN","",80,105);
    //sprintf(title,"#mu = %f\n#sigma = %f\nlg = %f",(float)voigtFit1->GetParameter(1),(float)voigtFit1->GetParameter(2),(float)voigtFit1->GetParameter(3)/2);
    zMassDistAfterHengne->Draw("P");
    sprintf(title,"#mu = %f",(float)voigtFit1->GetParameter(1) );
    l1->AddEntry(zMassDistBeforeAll,title,"l");
    sprintf(title,"#sigma = %f",(float)voigtFit1->GetParameter(2) );
    l1->AddEntry(zMassDistBeforeAll,title,"");
    sprintf(title,"lg = %f",(float)voigtFit1->GetParameter(3)/2 );
    l1->AddEntry(zMassDistBeforeAll,title,"");
    l1->Draw();
	c1->Print("ZMassDistAfterAll.png");
	c1->Clear();
    l1->Clear();
    
    //nice combination plot
    c1->SetMargin(0.12,0.02,0.09,0.025);
    int maxY = (int)zMassDistAfterHengne->GetBinContent(zMassDistAfterHengne->GetMaximumBin())*1.1;
    TH1D *dummy2 = new TH1D("dummy","",1,75,115);
    dummy2->SetTitle("Z-Peak Before and After Calibration;M_{Z} (GeV)");
    dummy2->GetYaxis()->SetRangeUser(0,maxY);
    dummy2->Draw(""); 
    
    zMassDistAfterHengne->SetTitle(";M_{ee} (GeV);Events/GeV");
    zMassDistAfterHengne->GetYaxis()->SetTitleOffset(1.4);
    zMassDistAfterHengne->GetYaxis()->SetTitleSize(0.045);
    zMassDistAfterHengne->SetMarkerStyle(21);
    zMassDistAfterHengne->SetLineColor(kBlue);
    zMassDistAfterHengne->SetMarkerColor(kBlue);
    zMassDistAfterHengne->SetMarkerSize(1.4);
    voigtFit1->SetParameter(1,90);
    voigtFit1->SetParLimits(1,90,100);
    voigtFit1->SetParameter(2,5);
    voigtFit1->SetParLimits(2,1,20);
    voigtFit1->FixParameter(3,4.9904);
    voigtFit1->SetParameter(4,1000);
    voigtFit1->SetParameter(5,-10);
    voigtFit1->SetLineWidth(6);
    voigtFit1->SetLineColor(kBlue);
    zMassDistAfterHengne->Fit(voigtFit1,"LMN","",80,105);
    sprintf(title,"After Hengne's Calibration" );
    l1->AddEntry(voigtFit1,title,"l");
    sprintf(title,"#mu = %.3g, #sigma = %.3g",(float)voigtFit1->GetParameter(1),(float)voigtFit1->GetParameter(2) );
    l1->AddEntry(voigtFit1,title,"");
    zMassDistAfterHengne->Draw("SAME P");
    voigtFit1->Draw("SAME C");
    
    zMassDistAfterAlexey->SetTitle(";M_{ee} (GeV);Events/GeV");
    zMassDistAfterAlexey->GetYaxis()->SetTitleOffset(1.4);
    zMassDistAfterAlexey->GetYaxis()->SetTitleSize(0.045);
    zMassDistAfterAlexey->SetMarkerStyle(22);
    zMassDistAfterAlexey->SetLineColor(kGreen+2);
    zMassDistAfterAlexey->SetMarkerColor(kGreen+2);
    zMassDistAfterAlexey->SetMarkerSize(1.3);
    voigtFit2->SetParameter(1,90);
    voigtFit2->SetParLimits(1,90,100);
    voigtFit2->SetParameter(2,5);
    voigtFit2->SetParLimits(2,1,20);
    voigtFit2->FixParameter(3,4.9904);
    voigtFit2->SetParameter(4,1000);
    voigtFit2->SetParameter(5,-10);
    voigtFit2->SetLineWidth(4);
    voigtFit2->SetLineColor(kGreen+2);
    zMassDistAfterAlexey->Fit(voigtFit2,"LMN","",80,105);
    sprintf(title,"After Alexey's Calibration" );
    l1->AddEntry(voigtFit2,title,"l");
    sprintf(title,"#mu = %.3g, #sigma = %.3g",(float)voigtFit2->GetParameter(1),(float)voigtFit2->GetParameter(2) );
    l1->AddEntry(voigtFit2,title,"");
    zMassDistAfterAlexey->Draw("SAME P");
    voigtFit2->Draw("SAME C");
    
    zMassDistBeforeAll->SetTitle(";;");
    zMassDistBeforeAll->SetMarkerStyle(20);
    zMassDistBeforeAll->SetLineColor(kRed+2);
    zMassDistBeforeAll->SetMarkerColor(kRed+2);
    zMassDistBeforeAll->SetMarkerSize(1.3);
    voigtFit3->SetParameter(1,90);
    voigtFit3->SetParLimits(1,90,100);
    voigtFit3->SetParameter(2,5);
    voigtFit3->SetParLimits(2,1,20);
    voigtFit3->FixParameter(3,4.9904);
    voigtFit3->SetParameter(4,1000);
    voigtFit3->SetParameter(5,-10);
    voigtFit3->SetLineWidth(4);
    voigtFit3->SetLineColor(kRed+2);
    zMassDistBeforeAll->Fit(voigtFit3,"LMN","",80,105);
    sprintf(title,"Before Calibration" );
    l1->AddEntry(voigtFit3,title,"l");
    sprintf(title,"#mu = %.3g, #sigma = %.3g",(float)voigtFit3->GetParameter(1),(float)voigtFit3->GetParameter(2) );
    l1->AddEntry(voigtFit3,title,"");
    zMassDistBeforeAll->Draw("SAME P");
    voigtFit3->Draw("SAME C");
    l1->Draw();
    //c1->Clear();//just in case    
    //zMassDistAfterHengne->DrawCopy("P");
    //voigtFit1->Draw("C");
    //zMassDistAfterHengne->DrawCopy(" SAME P");
    //zMassDistAfterAlexey->DrawCopy("SAME P");
    //voigtFit2->Draw("C");
    //zMassDistAfterAlexey->DrawCopy(" SAME P");
    //zMassDistBeforeAll->DrawCopy("SAME E");
    //voigtFit2->Draw("SAME");
    //zMassDistBeforeAll->DrawCopy("SAME E");
    //l1->Draw();
    
    //TLatex theTitle;
	//theTitle.SetTextSize(0.03);
	//theTitle.SetTextFont(42);
	//theTitle.SetNDC(true);    
    //theTitle.DrawLatex(0.16,0.9,"CMS preliminary");
    //theTitle.DrawLatex(0.16,0.85,"|#eta_{1}| < 2.5,  2.5 < |#eta_{2}| < 3.0");
    //theTitle.DrawLatex(0.16,0.8,"2012 -- #sqrt{s} = 8 TeV");
    //theTitle.DrawLatex(0.16,0.75,"#intL = 19 fb^{-1}, split sample calibration");
    
    //c1->Update();
    c1->Print("zMassCompare.png");
	c1->Clear();
    l1->Clear();
    
    c1->Close();

    
    return 0;
}