//This will be a test of whether I can use the Ntuple that Alexe's code made
//to make efficiency plots

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TF1.h"
//#include "TDirectory.h"
//#include "TObjArray.h"

#include "../src/ZEffTree.h"
//#include "tdrstyle.C"

void xSliceFitter(TH2*, TFile*);

int NtupleMassPlotter()
{
	gROOT->SetStyle("Plain");	
	//setTDRStyle();
	gStyle->SetErrorX(0);
	gStyle->SetOptStat(0);
	
	TCanvas* c1 = new TCanvas("C1", "c1", 800, 600);
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
	
	TFile* file = new TFile( "Zhists.root", "recreate" );
	if(file == NULL)
	{
		std::cout<<"Failed to create hists file. Exiting."<<std::endl;
		exit(1);
	}
	
	//define hists
	TH1* h1 = new TH1F("M_allNT", "All NT-SC Masses", 75, 0, 150 );
	TH1* h2 = new TH1F("M_goodNT", "Selected NT-SC Masses", 75, 0, 150 );
	
	TH1* h3 = new TH1F("MC_EBEB", "MC", 75, 0, 150 );
	TH1* h4 = new TH1F("MC_EBEE", "MC", 75, 0, 150 );
	TH1* h5 = new TH1F("MC_EBNT", "MC", 75, 0, 150 );
	TH1* h6 = new TH1F("MC_EENT", "MC", 75, 0, 150 );
	
	TH1* h7 = new TH1F("DT_EBEB", "Data", 75, 0, 150 );
	TH1* h8 = new TH1F("DT_EBEE", "Data", 75, 0, 150 );
	TH1* h9 = new TH1F("DT_EBNT", "Data", 75, 0, 150 );
	TH1* h10= new TH1F("DT_EENT", "Data", 75, 0, 150 );
	TH1* h11 = new TH1F("DT_GsfNT", "Data", 75, 0, 150 );
	TH1* h12 = new TH1F("DT_Wp80NTLoose", "Data", 75, 0, 150 );
	
	//define 2-D hists
	TH2* mcmvn1 = new TH2F("MC_EBEBMvN","MC, M vs Nvert EB-EB",30,60,120,30,0,30);
	TH2* mcmvn2 = new TH2F("MC_EBEEMvN","MC, M vs Nvert EB-EE",30,60,120,30,0,30);
	TH2* mcmvn3 = new TH2F("MC_EBNTMvN","MC, M vs Nvert EB-NT;nVert;Peak Location",30,60,120,30,0,30);
	TH2* mcmvn4 = new TH2F("MC_EENTMvN","MC, M vs Nvert EE-NT;nVert;Peak Location",30,60,120,30,0,30);
	
	TH2* dtmvn1 = new TH2F("DT_EBEBMvN","Data, M vs Nvert EB-EB",30,60,120,30,0,30);
	TH2* dtmvn2 = new TH2F("DT_EBEEMvN","Data, M vs Nvert EB-EE",30,60,120,30,0,30);
	TH2* dtmvn3 = new TH2F("DT_EBNTMvN","Data, M vs Nvert EB-NT;M_{ee};nVert",30,60,120,30,0,30);
	TH2* dtmvn4 = new TH2F("DT_EENTMvN","Data, M vs Nvert EE-NT;M_{ee};nVert",30,60,120,30,0,30);
	
	const std::string NT = "NTLooseElectronId-EtaDet";
	const std::string WP80 = "WP80";
	
	ZEffTree *ze1, *ze2;
    ze1 = new ZEffTree(*f1,false);
    ze2 = new ZEffTree(*f2,false);
    
    //run over MC ntuple
    bool run = true;
    while (run)
    {
    	ze1->Entries();
    	
        h1->Fill( ze1->reco.mz );
        
        if( ze1->reco.isSelected(1,NT) )
        {
        	h2->Fill( ze1->reco.mz );
        }
        
        if( fabs(ze1->reco.eta[0]) < 1.4442 )
        {
        	if( (fabs(ze1->reco.eta[1])<1.4442) && (ze1->reco.isSelected(1,WP80)) )
        	{
        		h3->Fill( ze1->reco.mz );
        		mcmvn1->Fill(ze1->reco.mz,ze1->reco.nverts);
        	}
        	
        	if( (fabs(ze1->reco.eta[1])>1.566) && (fabs(ze1->reco.eta[1])<2.5) && (ze1->reco.isSelected(1,WP80)) )
        	{
        		h4->Fill( ze1->reco.mz );
        		mcmvn2->Fill(ze1->reco.mz,ze1->reco.nverts);
        		//std::cout<<"Eta1 = "<<(float)abs(ze->reco.eta[1])<<std::endl;
        	}
        	
        	if( (fabs(ze1->reco.eta[1])>2.5) && (fabs(ze1->reco.eta[1])<3.0) && (ze1->reco.isSelected(1,NT)) )
        	{
        		h5->Fill( ze1->reco.mz );
        		mcmvn3->Fill(ze1->reco.mz,ze1->reco.nverts);
        	}
        }
        else if( (fabs(ze1->reco.eta[0])>1.566) && (fabs(ze1->reco.eta[0])<2.5) )
        {	
        	if( (fabs(ze1->reco.eta[1])>2.5) && (fabs(ze1->reco.eta[1])<3.0) && (ze1->reco.isSelected(1,NT)) )
        	{
        		h6->Fill( ze1->reco.mz );
        		mcmvn4->Fill(ze1->reco.mz,ze1->reco.nverts);
        	}
        }
        
    	run = ze1->GetNextEvent();
    }
    
    //now over Data
    run = true;
    while (run)
    {
        ze2->Entries();
        
        h11->Fill( ze2->reco.mz );
        
        if( ze2->reco.isSelected(1,NT) )
        {
        	h12->Fill( ze2->reco.mz );
        }
        
        if( fabs(ze2->reco.eta[0]) < 1.4442 )
        {
        	if( (fabs(ze2->reco.eta[1])<1.4442) && (ze2->reco.isSelected(1,WP80)) )
        	{
        		h7->Fill( ze2->reco.mz );
        		dtmvn1->Fill(ze2->reco.mz,ze2->reco.nverts);
        	}
        	
        	if( (fabs(ze2->reco.eta[1])>1.566) && (fabs(ze2->reco.eta[1])<2.5) && (ze2->reco.isSelected(1,WP80)) )
        	{
        		h8->Fill( ze2->reco.mz );
        		dtmvn2->Fill(ze2->reco.mz,ze2->reco.nverts);
        		//std::cout<<"Eta1 = "<<(float)abs(ze->reco.eta[1])<<std::endl;
        	}
        	
        	if( (fabs(ze2->reco.eta[1])>2.5) && (fabs(ze2->reco.eta[1])<3.0) && (ze2->reco.isSelected(1,NT)) )
        	{
        		h9->Fill( ze2->reco.mz );
        		dtmvn3->Fill(ze2->reco.mz,ze2->reco.nverts);
        	}
        }
        else if( (fabs(ze2->reco.eta[0])>1.566) && (fabs(ze2->reco.eta[0])<2.5) )
        {	
        	if( (fabs(ze2->reco.eta[1])>2.5) && (fabs(ze2->reco.eta[1])<3.0) && (ze2->reco.isSelected(1,NT)) )
        	{
        		h10->Fill( ze2->reco.mz );
        		dtmvn4->Fill(ze2->reco.mz,ze2->reco.nverts);
        	}
        }
        
    	run = ze2->GetNextEvent();
    }
    
    //write hists to file
    file->Write();
    
    h1->SetLineWidth(3);
    h1->SetLineColor(2);
    h2->SetLineWidth(3);
    h2->SetLineColor(2);
    h3->SetLineWidth(3);
    h3->SetLineColor(2);
    h4->SetLineWidth(3);
    h4->SetLineColor(2);
    h5->SetLineWidth(3);
    h5->SetLineColor(2);
    h6->SetLineWidth(3);
    h6->SetLineColor(2);
    
    //h7 = (TH1*)h7->Clone();
    h7->SetMarkerStyle(20);
    h7->SetMarkerColor(4);
    h7->SetMarkerSize(1);
    h8->SetMarkerStyle(20);
    h8->SetMarkerColor(4);
    h8->SetMarkerSize(1);
    h9->SetMarkerStyle(20);
    h9->SetMarkerColor(4);
    h9->SetMarkerSize(1);
    h10->SetMarkerStyle(20);
    h10->SetMarkerColor(4);
    h10->SetMarkerSize(1);
    h11->SetMarkerStyle(20);
    h11->SetMarkerColor(4);
    h11->SetMarkerSize(1);
    h12->SetMarkerStyle(20);
    h12->SetMarkerColor(4);
    h12->SetMarkerSize(1);
    
    /*c2->cd();
    h10->Draw("");
    c2->Print("sdla.png");*/
    
    //I really feel like this part should be done in a loop... But it's not.
    
    l1->AddEntry(h3, "MC");
    l1->AddEntry(h7, "Data");    
    h3->Scale(h7->Integral(35,55)/(h3->Integral(35,55)));
    h3->GetXaxis()->SetRangeUser(50, 150);
    h3->GetYaxis()->SetRangeUser(0, 1.3*h3->GetBinContent(h3->GetMaximumBin()));
    h3->SetTitle("Z mass EB+EB(WP80); Mass (GeV)");
    h3->Draw("hist");
    h7->Draw("same E1");
    l1->Draw();
    c1->Update();
    c1->Draw();
    c1->Print("Efficiencies/ZMassPlots/ZmassMCvsData_EBEBWP80.png");
    l1->Clear();
    c1->Clear();
    
    l1->AddEntry(h4, "MC");
    l1->AddEntry(h8, "Data");    
    h4->Scale((h8->Integral(30,60))/(h4->Integral(30,60)));
    h4->GetXaxis()->SetRangeUser(50, 150);
    h4->GetYaxis()->SetRangeUser(0, 1.3*h4->GetBinContent(h4->GetMaximumBin()));
    h4->SetTitle("Z mass EB+EE(WP80); Mass (GeV)");
    h4->Draw("hist");
    h8->Draw("same E1");
    l1->Draw();
    c1->Update();
    c1->Print("Efficiencies/ZMassPlots/ZmassMCvsData_EBEEWP80.png");
    l1->Clear();
    c1->Clear();
    
    l1->AddEntry(h5, "MC");
    l1->AddEntry(h9, "Data");    
    h5->Scale((h9->Integral(30,60))/(h5->Integral(30,60)));
    h5->GetXaxis()->SetRangeUser(50, 150);
    h5->GetYaxis()->SetRangeUser(0, 1.3*h5->GetBinContent(h5->GetMaximumBin()));
    h5->SetTitle("Z mass EB+NTLoose; Mass (GeV)");
    h5->Draw("hist");
    h9->Draw("same E");
    l1->Draw();
    c1->Update();
    c1->Print("Efficiencies/ZMassPlots/ZmassMCvsData_EBNTLoose.png");
    l1->Clear();
    c1->Clear();
    
    l1->AddEntry(h6, "MC");
    l1->AddEntry(h10, "Data");    
    h6->Scale((h10->Integral(30,60))/(h6->Integral(30,60)));
    h6->GetXaxis()->SetRangeUser(50, 150);
    h6->GetYaxis()->SetRangeUser(0, 1.3*h6->GetBinContent(h6->GetMaximumBin()));
    h6->SetTitle("Z mass EE+NTLoose; Mass (GeV)");
    h6->Draw("hist");
    h10->Draw("same E");
    l1->Draw();
    c1->Update();
    c1->Print("Efficiencies/ZMassPlots/ZmassMCvsData_EENTLoose.png");
    l1->Clear();
    c1->Clear();
    
    l1->AddEntry(h1, "MC");
    l1->AddEntry(h11, "Data");    
    h1->Scale((h11->Integral(30,60))/(h1->Integral(30,60)));
    h1->GetXaxis()->SetRangeUser(50, 150);
    h1->GetYaxis()->SetRangeUser(0, 1.3*h1->GetBinContent(h1->GetMaximumBin()));
    h1->SetTitle("Z mass GSF+NT(acc.); Mass (GeV)");
    h1->Draw("hist");
    h11->Draw("same E");
    l1->Draw();
    c1->Update();
    c1->Print("Efficiencies/ZMassPlots/ZmassMCvsData_GsfNT.png");
    l1->Clear();
    c1->Clear();
    
    l1->AddEntry(h2, "MC");
    l1->AddEntry(h12, "Data");    
    h2->Scale((h12->Integral(30,60))/(h2->Integral(30,60)));
    h2->GetXaxis()->SetRangeUser(50, 150);
    h2->GetYaxis()->SetRangeUser(0, 1.3*h2->GetBinContent(h2->GetMaximumBin()));
    h2->SetTitle("Z mass WP80+NTLoose; Mass (GeV)");
    h2->Draw("hist");
    h12->Draw("same E");
    l1->Draw();
    c1->Update();
    c1->Print("Efficiencies/ZMassPlots/ZmassMCvsData_WP80NTLoose.png");
    l1->Clear();
    c1->Clear();
    
    mcmvn1->Draw("COLZ");
    c1->Print("Efficiencies/ZMassPlots/MZcsNvert_EBEB_MC.png");
    c1->Clear();
    mcmvn2->Draw("COLZ");
    c1->Print("Efficiencies/ZMassPlots/MZcsNvert_EBEE_MC.png");
    c1->Clear();
    mcmvn3->Draw("COLZ");
    c1->Print("Efficiencies/ZMassPlots/MZcsNvert_EBNT_MC.png");
    c1->Clear();
    mcmvn4->Draw("COLZ");
    c1->Print("Efficiencies/ZMassPlots/MZcsNvert_EENT_MC.png");
    c1->Clear();
    dtmvn1->Draw("COLZ");
    c1->Print("Efficiencies/ZMassPlots/MZcsNvert_EBEB_Data.png");
    c1->Clear();
    dtmvn2->Draw("COLZ");
    c1->Print("Efficiencies/ZMassPlots/MZcsNvert_EBEE_Data.png");
    c1->Clear();
    dtmvn3->Draw("COLZ");
    c1->Print("Efficiencies/ZMassPlots/MZcsNvert_EBNT_Data.png");
    c1->Clear();
    dtmvn4->Draw("COLZ");
    c1->Print("Efficiencies/ZMassPlots/MZcsNvert_EENT_Data.png");
    c1->Clear();
    
    //let's try to make some profiles
    TProfile* profile;
    profile = mcmvn1->ProfileY();
    profile->SetMarkerStyle(20);
    profile->GetYaxis()->SetRangeUser(80,100);
    profile->Draw("E");
    c1->Print("Efficiencies/ZMassPlots/EBEB_profile_MC.png");
    c1->Clear();
    
    profile = mcmvn2->ProfileY();
    profile->SetMarkerStyle(20);
    profile->GetYaxis()->SetRangeUser(80,100);
    profile->Draw("E");
    c1->Print("Efficiencies/ZMassPlots/EBEE_profile_MC.png");
    c1->Clear();
    
    profile = mcmvn3->ProfileY();
    profile->SetMarkerStyle(20);
    profile->GetYaxis()->SetRangeUser(80,100);
    profile->Draw("E");
    c1->Print("Efficiencies/ZMassPlots/EBNT_profile_MC.png");
    c1->Clear();
    
    profile = mcmvn4->ProfileY();
    profile->SetMarkerStyle(20);
    profile->GetYaxis()->SetRangeUser(80,100);
    profile->Draw("E");
    c1->Print("Efficiencies/ZMassPlots/EENT_profile_MC.png");
    c1->Clear();
    
    profile = dtmvn1->ProfileY();
    profile->SetMarkerStyle(20);
    profile->GetYaxis()->SetRangeUser(80,100);
    profile->Draw("E");
    c1->Print("Efficiencies/ZMassPlots/EBEB_profile_Data.png");
    c1->Clear();
    
    profile = dtmvn2->ProfileY();
    profile->SetMarkerStyle(20);
    profile->GetYaxis()->SetRangeUser(80,100);
    profile->Draw("E");
    c1->Print("Efficiencies/ZMassPlots/EBEE_profile_Data.png");
    c1->Clear();
    
    profile = dtmvn3->ProfileY();
    profile->SetMarkerStyle(20);
    profile->GetYaxis()->SetRangeUser(80,100);
    profile->GetYaxis()->SetTitle("Peak Location");
    profile->Draw("E");
    c1->Print("Efficiencies/ZMassPlots/EBNT_profile_Data.png");
    c1->Clear();
    
    profile = dtmvn4->ProfileY();
    profile->SetMarkerStyle(20);
    profile->GetYaxis()->SetRangeUser(80,100);
    profile->GetYaxis()->SetTitle("Peak Location");
    profile->Draw("E");
    c1->Print("Efficiencies/ZMassPlots/EENT_profile_Data.png");
    c1->Clear();
    
    //even better, let's do some slice fitting
    //turns out, impossible to get built in slice-fitter to work, so going to write my own.
    
    //c1->Print("EENT_MeanVsNvert_Data.png");
    //c1->Clear();
    
    c1->Close();
    
    file->Write();
	
	return 0;
}



void xSliceFitter(TH2* h2, TFile* file)
{
	TH1** hists = new TH1*[h2->GetNbinsY()];
	TF1** fits  = new TF1*[h2->GetNbinsY()];
	
	TH1* Mean = new TH1F("GFMean","Gauss Fit Mean vs. Nvert",h2->GetNbinsY(),0,h2->GetNbinsY());
	TH1* Sigma = new TH1F("GFSigma","Gauss Fit Sigma vs. Nvert",h2->GetNbinsY(),0,h2->GetNbinsY());
	
	for(int i=0; i < h2->GetNbinsY(); i++ )
	{
		char temp[128];
		sprintf(temp,"slice_%i", i);
		hists[i] = (TH1*)h2->ProjectionX(temp, i, i);
		sprintf(temp,"fit_%i", i);
		fits[i] = new TF1(temp, "gaus", 60, 120);
		hists[i]->Fit(fits[i], "", "", 60, 120);
		Mean->SetBinContent(i,fits[i]->GetParameter(1));
		Sigma->SetBinContent(i,fits[i]->GetParameter(2));
		//if(i==10) break;
	}
	Mean->Draw();
}



















