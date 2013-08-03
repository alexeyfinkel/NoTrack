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
#include "TProfile.h"
#include "TF1.h"

#include "../src/ZEffTree.h"

void ySliceFitter(TH2D*, TH1**);

int pileupDep()
{
	gROOT->SetStyle("Plain");	
	gStyle->SetErrorX(0);
	gStyle->SetOptStat(0);
	
	TCanvas* c1 = new TCanvas("C1", "c1", 1200, 900);
	c1->cd();
	//c1->SetLeftMargin(0.08);
	c1->SetRightMargin(0.12);
	
	TLegend* l1 = new TLegend(0.65,0.7,0.9,0.9);
	l1->SetFillColor(10);
	
	//grab MC ntuple (just in case I need it later)
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
	
	//make the var-nV hists:
	TH2D* rNineNv = new TH2D("rNineNv","R9-nVertex;nVertex;R9",30,0,30,40,0.85,1.05);
	TH2D* sIeIeNv = new TH2D("sIeIeNv","#sigma_{i#etai#eta}-nVertex;nVertex;#sigma_{i#etai#eta}",30,0,30,40,0.015,0.035);
	TH2D* HoEmNv = new TH2D("HoEmNv","H/EM-nVertex;nVertex;H/EM",30,0,30,20,0,0.1);
	TH2D* hIsoNv = new TH2D("hIsoNv","HcalIso-nVertex;nVertex;HcalIso",30,0,30,30,0,0.15);
	TH2D* eIsoNv = new TH2D("eIsoNv","EcalIso-nVertex;nVertex;EcalIso",30,0,30,40,0,0.05);
	
	//make data ntuple:
	ZEffTree* ze2 = new ZEffTree(*f2,false);
	//Fill the data histograms:
	for(int event=0; event<ze2->Entries(); event++ ) //fill the data hists
	{
		if( (fabs(ze2->reco.eta[1])>2.5) && (fabs(ze2->reco.eta[1])<3.0) ) 
		{	//fill uncut variable hists
			rNineNv->Fill(ze2->reco.nverts,ze2->reco.RNine[1]);
			sIeIeNv->Fill(ze2->reco.nverts,ze2->reco.Sieie[1]);
			HoEmNv->Fill(ze2->reco.nverts,ze2->reco.HoEM[1]);
			hIsoNv->Fill(ze2->reco.nverts,ze2->reco.Hiso[1]);
			eIsoNv->Fill(ze2->reco.nverts,ze2->reco.Eiso[1]);
		}
		ze2->GetNextEvent();
	}
	
	rNineNv->GetYaxis()->SetLabelSize(0.03);
	rNineNv->GetZaxis()->SetLabelSize(0.03);
	rNineNv->Draw("colz");
	c1->Print("Efficiencies/PUdependence/rNineNv.png");
	c1->Clear();
	sIeIeNv->GetYaxis()->SetTitleSize(0.05);
	sIeIeNv->GetYaxis()->SetLabelSize(0.03);
	sIeIeNv->GetZaxis()->SetLabelSize(0.03);
	sIeIeNv->Draw("colz");
	c1->Print("Efficiencies/PUdependence/sIeIeNv.png");
	c1->Clear();
	HoEmNv->GetYaxis()->SetLabelSize(0.03);
	HoEmNv->GetZaxis()->SetLabelSize(0.03);
	HoEmNv->Draw("colz");
	c1->Print("Efficiencies/PUdependence/HoEmNv.png");
	c1->Clear();
	hIsoNv->GetYaxis()->SetLabelSize(0.03);
	hIsoNv->GetZaxis()->SetLabelSize(0.03);
	hIsoNv->Draw("colz");
	c1->Print("Efficiencies/PUdependence/hIsoNv.png");
	c1->Clear();
	eIsoNv->GetYaxis()->SetLabelSize(0.03);
	eIsoNv->GetZaxis()->SetLabelSize(0.03);
	eIsoNv->Draw("colz");
	c1->Print("Efficiencies/PUdependence/eIsoNv.png");
	c1->Clear();
	
	//////////profiles/////////
	
	TH1D *MeanR9 = new TH1D();
	TProfile* profile;
    profile = rNineNv->ProfileX();
    profile->SetTitle(";nVertex;Peak R9");
    profile->SetMarkerStyle(20);
    profile->GetYaxis()->SetRangeUser(0.85,1.05);
    profile->GetYaxis()->SetLabelSize(0.03);
    profile->Draw("E");
    c1->Print("Efficiencies/PUdependence/PeakR9vsPU.png");
    c1->Clear();
    ySliceFitter(rNineNv,(TH1**)&MeanR9);
    MeanR9->SetTitle(";nVertex;Gauss Mean R9");
    MeanR9->SetMarkerStyle(20);
    MeanR9->GetYaxis()->SetRangeUser(0.85,1.1);
    MeanR9->GetYaxis()->SetLabelSize(0.03);
    MeanR9->Draw("E");
    c1->Print("Efficiencies/PUdependence/MeanR9vsPU.png");
    c1->Clear();
    
    TH1D *MeanSieie = new TH1D();
    profile = sIeIeNv->ProfileX();
    profile->SetTitle(";nVertex;Peak #sigma_{i#etai#eta}");
    profile->SetMarkerStyle(20);
    profile->GetYaxis()->SetRangeUser(0.015,0.029);
    profile->GetYaxis()->SetLabelSize(0.03);
    profile->Draw("E");
    c1->Print("Efficiencies/PUdependence/PeakSieieVsPU.png");
    c1->Clear();
    ySliceFitter(sIeIeNv,(TH1**)&MeanSieie);
    MeanSieie->SetTitle(";nVertex;Gauss Mean #sigma_{i#etai#eta}");
    MeanSieie->SetMarkerStyle(20);
    MeanSieie->GetYaxis()->SetRangeUser(0.02,0.035);
    MeanSieie->GetYaxis()->SetLabelSize(0.03);
    MeanSieie->Draw("E");
    c1->Print("Efficiencies/PUdependence/MeanSisisVsPU.png");
    c1->Clear();
    
    TH1D *MeanHoEM = new TH1D();
    profile = HoEmNv->ProfileX();
    profile->SetTitle(";nVertex;Peak H/EM");
    profile->SetMarkerStyle(20);
    profile->GetYaxis()->SetRangeUser(0.0,0.05);
    profile->GetYaxis()->SetLabelSize(0.03);
    profile->Draw("E");
    c1->Print("Efficiencies/PUdependence/PeakHoEmVsPU.png");
    c1->Clear();
    ySliceFitter(sIeIeNv,(TH1**)&MeanHoEM);
    MeanHoEM->SetTitle(";nVertex;Gauss Mean H/EM");
    MeanHoEM->SetMarkerStyle(20);
    MeanHoEM->GetYaxis()->SetRangeUser(0.0,0.05);
    MeanHoEM->GetYaxis()->SetLabelSize(0.03);
    MeanHoEM->Draw("E");
    c1->Print("Efficiencies/PUdependence/MeanHoEMvsPU.png");
    c1->Clear();
    
    TH1D *MeanHiso = new TH1D();
    profile = hIsoNv->ProfileX();
    profile->SetTitle(";nVertex;Peak Hcal Iso");
    profile->SetMarkerStyle(20);
    profile->GetYaxis()->SetRangeUser(0,0.11);
    profile->GetYaxis()->SetLabelSize(0.03);
    profile->Draw("E");
    c1->Print("Efficiencies/PUdependence/PeakHisoVsPU.png");
    c1->Clear();
    ySliceFitter(sIeIeNv,(TH1**)&MeanHiso);
    MeanHiso->SetTitle(";nVertex;Gauss Mean Hcal Iso");
    MeanHiso->SetMarkerStyle(20);
    MeanHiso->GetYaxis()->SetRangeUser(0,0.11);
    MeanHiso->GetYaxis()->SetLabelSize(0.03);
    MeanHiso->Draw("E");
    c1->Print("Efficiencies/PUdependence/MeanHisoVsPU.png");
    c1->Clear();
    
    TH1D *MeanEiso = new TH1D();
    profile = eIsoNv->ProfileX();
    profile->SetTitle(";nVertex;Peak Hcal Iso");
    profile->SetMarkerStyle(20);
    profile->GetYaxis()->SetRangeUser(0,0.035);
    profile->GetYaxis()->SetLabelSize(0.03);
    profile->Draw("E");
    c1->Print("Efficiencies/PUdependence/PeakEisoVsPU.png");
    c1->Clear();
    ySliceFitter(sIeIeNv,(TH1**)&MeanEiso);
    MeanEiso->SetTitle(";nVertex;Gauss Mean Ecal Iso");
    MeanEiso->SetMarkerStyle(20);
    MeanEiso->GetYaxis()->SetRangeUser(0,0.035);
    MeanEiso->GetYaxis()->SetLabelSize(0.03);
    MeanEiso->Draw("E");
    c1->Print("Efficiencies/PUdependence/MeanEisoVsPU.png");
    c1->Clear();
	
	c1->Close();
	return 0;
}

void ySliceFitter(TH2D* h2, TH1** Mean)
{
	TF1* fit  = new TF1("fit", "gaus",0,2 );	
	*Mean = new TH1D("Mean","Name",h2->GetNbinsX(),0,h2->GetNbinsX());	
	for(int i=1; i < h2->GetNbinsX(); i++ )
	{
		TH1D* hist;
		hist = h2->ProjectionY("HOTAHAKER", i+1, i+1);
		hist->Fit(fit, "QN", "", hist->GetBinCenter(1),hist->GetBinCenter(hist->GetNbinsX()));
		(*Mean)->SetBinContent(i,fit->GetParameter(1));
		(*Mean)->SetBinError(i,fit->GetParameter(2));
		hist->Delete();
	}
}

int main()
{
	pileupDep();
}














