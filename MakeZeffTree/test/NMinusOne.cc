//this script will (hopefully) do N-1 hists for all the shower shape vars.


//the following may not be true anymore; but run in Root at your own risk
//**********************************************************************************************
//***WARNING***: this script does not compile from inside root! 
//               compile with g++ `root-config --cflags --glibs` NMinusOne.cc ../src/ZEffTree.cc
//               (run with ./a.out)
//**********************************************************************************************



#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <math.h>
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TColor.h"

#include "../src/ZEffTree.h"

bool passCuts(std::vector<double>& vars, std::vector<double>& lowCut, std::vector<double>& highCut,std::vector<bool>& cuts);

int NMinusOne()
{
	gROOT->SetStyle("Plain");	
	gStyle->SetErrorX(0);
	gStyle->SetOptStat(0);
	
	TCanvas* c1 = new TCanvas("C1", "c1", 1200, 900);
	c1->cd();
	
	TLegend* l1 = new TLegend(0.65,0.7,0.9,0.9);
	l1->SetFillColor(10);
	
	//grab MC ntuple
	TFile* f1 = new TFile("/local/cms/user/finkel/NoTrack/Ntuple/DE_MC2012_AlexTrgr3/DE_MC2012_AlexTrgr_newerCuts.root");	
	if(f1 == NULL)
	{
		std::cout<<"Failed to open MC file. Exiting."<<std::endl;
		return 1;
	}
	
	//grab a data ntuple
	TFile* f2 = new TFile("/local/cms/user/finkel/NoTrack/Ntuple/DE_ReReco_NewerCuts/DE_ReReco2012FULL_Jan13_AlexTrgr_newerCuts_001.root");	
	if(f2 == NULL)
	{
		std::cout<<"Failed to open Data file. Exiting."<<std::endl;
		return 1;
	}
	
	//make variable segmentation
	std::vector<double> lowBound(6,0), highBound(6,0), lowCut(6,0), highCut(6,0); //order is R9, sIeIe, H/EM, hIso, eIso, Pt
	
	lowBound[0] = 0.65;highBound[0]=1.15; lowCut[0] = 0.89; highCut[0] = 1.02;
	lowBound[1] = 0.015;highBound[1]=0.04; lowCut[1] = 0; highCut[1] = 0.029;
	lowBound[2] = 0;highBound[2]=0.1; lowCut[2] = 0; highCut[2] = 0.05;
	lowBound[3] = 0;highBound[3]=0.25; lowCut[3] = 0; highCut[3] = 0.11;
	lowBound[4] = 0;highBound[4]=0.05; lowCut[4] = 0; highCut[4] = 0.035;	
	lowBound[5] = 0;highBound[5]=100; lowCut[5] = 20; highCut[5] = 160;
	/*lowBound[0] = -0.1;highBound[0]=1.30; lowCut[0] = -1; highCut[0] = 2;
	lowBound[1] = -0.1;highBound[1]=0.5; lowCut[1] = -1; highCut[1] = 1;
	lowBound[2] = -0.1;highBound[2]=0.5; lowCut[2] = -1; highCut[2] = 1;
	lowBound[3] = -0.1;highBound[3]=0.5; lowCut[3] = -1; highCut[3] = 1;
	lowBound[4] = -0.1;highBound[4]=0.5; lowCut[4] = -1; highCut[4] = 1;	
	lowBound[5] = -0.1;highBound[5]=200; lowCut[5] = -1; highCut[5] = 200;*/
	
	//make hist collections:
	std::vector<TH1D*> VarNoCuts, VarNm1, VarFullCuts, VarOnlyCut;
	char name0[100], title0[200], name1[100], title1[200], name2[100], title2[200], name3[100], title3[200];
	for(int var =0; var<6; var++)
	{
		switch (var)
		{
			case 0:
				sprintf(name0,"rNine0");
				sprintf(title0,"R9;R9");
				sprintf(name1,"rNine1");
				sprintf(title1,"R9, N-1 Cuts;R9");
				sprintf(name2,"rNine2");
				sprintf(title2,"R9, All Cuts;R9");
				sprintf(name3,"rNine3");
				sprintf(title3,"R9, R9 Cut Only;R9");
				break;
			case 1:
				sprintf(name0,"sIeIe0");
				sprintf(title0,"#sigma_{i#eta,i#eta};#sigma_{i#eta,i#eta}");
				sprintf(name1,"sIeIe1");
				sprintf(title1,"#sigma_{i#eta,i#eta}, N-1 Cuts;#sigma_{i#eta,i#eta}");
				sprintf(name2,"sIeIe2");
				sprintf(title2,"#sigma_{i#eta,i#eta}, All Cuts;#sigma_{i#eta,i#eta}");
				sprintf(name3,"sIeIe3");
				sprintf(title3,"#sigma_{i#eta,i#eta}, #sigma_{i#eta,i#eta} Cut Only;#sigma_{i#eta,i#eta}");
				break;
			case 2:
				sprintf(name0,"HoEM0");
				sprintf(title0,"H/EM;H/EM");
				sprintf(name1,"HoEM1");
				sprintf(title1,"H/EM, N-1 Cuts;H/EM");
				sprintf(name2,"HoEM2");
				sprintf(title2,"H/EM, All Cuts;H/EM");
				sprintf(name3,"HoEM3");
				sprintf(title3,"H/EM, H/EM Cut Only;H/EM");
				break;
			case 3:
				sprintf(name0,"hIso0");
				sprintf(title0,"HcalIso;HcalIso");
				sprintf(name1,"hIso1");
				sprintf(title1,"HcalIso, N-1 Cuts;HcalIso");
				sprintf(name2,"hIso2");
				sprintf(title2,"HcalIso, All Cuts;HcalIso");
				sprintf(name3,"hIso3");
				sprintf(title3,"HcalIso, HcalIso Cut Only;HcalIso");
				break;
			case 4:
				sprintf(name0,"eIso0");
				sprintf(title0,"EcalIso;EcalIso");
				sprintf(name1,"eIso1");
				sprintf(title1,"EcalIso, N-1 Cuts;EcalIso");
				sprintf(name2,"eIso2");
				sprintf(title2,"EcalIso, All Cuts;EcalIso");
				sprintf(name3,"eIso3");
				sprintf(title3,"EcalIso, EcalIso Cut Only;EcalIso");
				break;
			case 5:
				sprintf(name0,"Pt0");
				sprintf(title0,"Pt;Pt");
				sprintf(name1,"Pt1");
				sprintf(title1,"Pt, N-1 Cuts;Pt");
				sprintf(name2,"Pt2");
				sprintf(title2,"Pt, All Cuts;Pt");
				sprintf(name3,"Pt3");
				sprintf(title3,"Pt, Pt Cut Only;Pt");
				break;
		}
		
		VarNoCuts.push_back( new TH1D(name0,title0,50,lowBound[var],highBound[var]) );
		VarNm1.push_back( new TH1D(name1,title1,50,lowBound[var],highBound[var]) );
		VarFullCuts.push_back( new TH1D(name2,title2,50,lowBound[var],highBound[var]) );
		VarOnlyCut.push_back( new TH1D(name3,title3,50,lowBound[var],highBound[var]) );
	}
	
	//mass comparison hists:
	TH1D* mBefore = new TH1D("mBefore",";M_{ee};Events/2GeV",50,50,150);
	TH1D* mAfter = new TH1D("mAfter",";M_{ee};Events/2GeV",50,50,150);
	TH1D* mMC = new TH1D("mMC",";M_{ee};Events/2GeV",50,50,150);
	
	ZEffTree *ze1;
	ze1 = new ZEffTree(*f1,false);
	for(int event=0; event<ze1->Entries(); event++ ) //fill the MC hists
	{
		if( (fabs(ze1->reco.eta[1])>2.5) && (fabs(ze1->reco.eta[1])<3.0) && (ze1->reco.isSelected(1,"NTLooseElectronId-EtaDet"))
			)
		{
			mMC->Fill(ze1->reco.mz);
		}
		ze1->GetNextEvent();
	}
	
	//make data ntuple:
	ZEffTree* ze2 = new ZEffTree(*f2,false);
	
	std::vector<double> vars(6,0);
	//Fill the data histograms:
	for(int event=0; event<ze2->Entries(); event++ ) //fill the data hists
	{
		if( (fabs(ze2->reco.eta[1])>2.5) && (fabs(ze2->reco.eta[1])<3.0) ) 
		{	//fill uncut variable hists
			vars[0] = ze2->reco.RNine[1];
			vars[1] = ze2->reco.Sieie[1];
			vars[2] = ze2->reco.HoEM[1];
			vars[3] = ze2->reco.Hiso[1];
			vars[4] = ze2->reco.Eiso[1];
			vars[5] = ze2->reco.pt[1];
			//fill mass hists
			mBefore->Fill(ze2->reco.mz);
			for(int var=0; var<(int)vars.size(); var++)
			{
				std::vector<bool> cuts(6,false);
				cuts[var] = true;
				if ( passCuts(vars,lowCut,highCut,cuts) )
				{
					VarOnlyCut[var]->Fill(vars[var]);
				}
				cuts = std::vector<bool>(6,true);
				cuts[var] = false;
				VarNoCuts[var]->Fill(vars[var]);
				if ( passCuts(vars,lowCut,highCut,cuts) )
				{
					VarNm1[var]->Fill(vars[var]);
				}
				cuts = std::vector<bool>(6,true);
				if ( passCuts(vars,lowCut,highCut,cuts) )
				{
					VarFullCuts[var]->Fill(vars[var]);
				}
			}
			//sanity check: explicit cuts vs. "isSelected"
			std::vector<bool> cuts(6,true);
			if( ze2->reco.isSelected(1,"NTLooseElectronId-EtaDet") ) mAfter->Fill(ze2->reco.mz);
		}
		ze2->GetNextEvent();
	}
	
	for(int var=0; var<(int)vars.size(); var++)
	{
		char filename[100];
		if(var==0) sprintf(filename,"Efficiencies/NMinusOne/Rnine.png");
		if(var==1) sprintf(filename,"Efficiencies/NMinusOne/Sieie.png");
		if(var==2) sprintf(filename,"Efficiencies/NMinusOne/HoEM.png");
		if(var==3) sprintf(filename,"Efficiencies/NMinusOne/Hiso.png");
		if(var==4) sprintf(filename,"Efficiencies/NMinusOne/Eiso.png");
		if(var==5) sprintf(filename,"Efficiencies/NMinusOne/Pt.png");
		VarNoCuts[var]->SetLineColor(kBlue);
		VarNoCuts[var]->SetLineWidth(2);
		VarNoCuts[var]->Draw("hist");
		VarNm1[var]->SetFillColor(kRed);
		VarNm1[var]->Draw("same hist");
		VarFullCuts[var]->SetFillColor(kGreen+2);
		VarFullCuts[var]->Draw("same hist");		
		l1->AddEntry(VarNoCuts[var],"No Cuts");
		l1->AddEntry(VarNm1[var],"N-1 Cuts");
		l1->AddEntry(VarFullCuts[var],"All Cuts");
		l1->Draw();		
		c1->Print(filename);
		c1->Clear();
		l1->Clear();
		
		if(var==0) sprintf(filename,"Efficiencies/NMinusOne/RnineOnly.png");
		if(var==1) sprintf(filename,"Efficiencies/NMinusOne/SieieOnly.png");
		if(var==2) sprintf(filename,"Efficiencies/NMinusOne/HoEMOnly.png");
		if(var==3) sprintf(filename,"Efficiencies/NMinusOne/HisoOnly.png");
		if(var==4) sprintf(filename,"Efficiencies/NMinusOne/EisoOnly.png");
		if(var==5) sprintf(filename,"Efficiencies/NMinusOne/PtOnly.png");
		VarNoCuts[var]->SetFillColor(kRed);
		VarNoCuts[var]->Draw("hist");
		VarOnlyCut[var]->SetFillColor(kBlue);
		VarOnlyCut[var]->Draw("same hist");		
		l1->AddEntry(VarNoCuts[var],"No Cuts");
		l1->AddEntry(VarOnlyCut[var],"N-1 Cuts");
		l1->Draw();		
		c1->Print(filename);
		c1->Clear();
		l1->Clear();
	}
	
	mBefore->SetMarkerSize(1.5);
	mBefore->SetMarkerColor(4);
	mBefore->SetMarkerStyle(20);
	mBefore->Draw("P");
	mAfter->SetMarkerSize(1.5);
	mAfter->SetMarkerColor(kRed+1);
	mAfter->SetMarkerStyle(20);
	mAfter->Draw("same P");
	mMC->SetLineWidth(2);
	mMC->Scale( mBefore->GetBinContent(mBefore->GetMaximumBin()) / mMC->GetBinContent(mMC->GetMaximumBin()) );
	mMC->Draw("same hist");
	
	std::cout<<mBefore->Integral(11,31)<<" events pass trigger\n"<<mAfter->Integral(11,31)<<" events pass NT selection."<<std::endl;
	
	l1->AddEntry(mBefore,"Before","P");
	l1->AddEntry(mAfter,"After","P");
	l1->AddEntry(mMC,"MC","L");
	l1->Draw();
	c1->Update();
	c1->Print("Efficiencies/NMinusOne/MtriggerVsNTcut.png");
	
	c1->Close();
	
	return 0;	
}




bool passCuts(std::vector<double>& vars, std::vector<double>& lowCut, std::vector<double>& highCut,std::vector<bool>& cuts)
{
	for(unsigned int i=0;i<cuts.size();i++)
	{
		if(!cuts[i]) continue;
		if( (vars[i]<lowCut[i]) || (vars[i]>highCut[i]) )
		{
			return false;
			//std::cout<<"Failed Event Found:\nR9="<<vars[0]<<" sieie="<<vars[1]<<" H/EM="<<vars[2]<<" hiso="<<vars[3]<<" eiso="<<vars[4]<<" pt="<<vars[5]<<std::endl;			
		}
	}	
	return true;
}




int main()
{
	return NMinusOne();
}

























	
	
	
	
