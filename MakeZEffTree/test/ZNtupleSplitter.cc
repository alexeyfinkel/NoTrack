//this script is intended to take a Z Ntuple and split it into two for validation (even/odd events)


//***WARNING***: this script does not compile from inside root! 
//               compile with g++ `root-config --cflags --glibs` ZNtupleSplitter.cc ../src/ZEffTree.cc
//               (run with ./a.out)
//**********************************************************************************************

#include <stdio.h>
#include <iostream>
#include <string>
#include <TFile.h>
#include <TH1.h>

#include "TTree.h"
//#include "TVirtualIndex.h"

//#include "../src/ZEffTree.h"

//void copyZ(const ZEffTree* oldZ, ZEffTree* newZ);

int main()
{
    //make product TFiles:
     TFile* f0 = new TFile("/local/cms/user/finkel/NoTrack/Ntuple/DE_ReReco_2012Full_WithClusters_001.root");	
	if(!f0->IsOpen())
	{
		std::cout<<"Failed to open Original file. Exiting."<<std::endl;
		return 1;
	}
    
    TObject* obj = f0->Get("ZEffs");

    TFile* f1 = new TFile("/local/cms/user/finkel/NoTrack/Ntuple/DE_ReReco_2012Full_WithClusters_001_ODDS.root", "recreate");	
	if(!f1->IsOpen())
	{
		std::cout<<"Failed to create Odds file. Exiting."<<std::endl;
		return 1;
	}
    TTree* ot1 = ((TTree*)obj)->CloneTree(0);
    ot1->SetAutoSave(((TTree*)obj)->GetEntries() + 1);
    ((TTree*)obj)->GetListOfClones()->Remove(ot1);
    ((TTree*)obj)->ResetBranchAddresses();
    ot1->ResetBranchAddresses();
    
    TTree* tm = (TTree*)ot1;
    TTree* ts = (TTree*)obj;
    tm->CopyAddresses(ts);
    for (int i = 0; i < ts->GetEntries(); i += 2)
    {
        ts->GetEntry(i);
        tm->Fill();
    }
    ts->ResetBranchAddresses();
    
    f1->cd();
    ot1->Write();
    f1->Close();

    TFile* f2 = new TFile("/local/cms/user/finkel/NoTrack/Ntuple/DE_ReReco_2012Full_WithClusters_001_EVENS.root", "recreate");	
	if(!f2->IsOpen())
	{
		std::cout<<"Failed to create Evens file. Exiting."<<std::endl;
		return 1;
    }
    TTree* ot2 = ((TTree*)obj)->CloneTree(0);
    ot2->SetAutoSave(((TTree*)obj)->GetEntries() + 1);
    ((TTree*)obj)->GetListOfClones()->Remove(ot2);
    ((TTree*)obj)->ResetBranchAddresses();
    ot2->ResetBranchAddresses();
    
    tm = (TTree*)ot2;
    ts = (TTree*)obj;
    tm->CopyAddresses(ts);
    for (int i = 1; i < ts->GetEntries(); i += 2)
    {
        ts->GetEntry(i);
        tm->Fill();
    }
    ts->ResetBranchAddresses();
    
    f2->cd();
    ot2->Write();
    f2->Close();

    
    //ZEffTree* ze0 = new ZEffTree(*f0,false);
    //ZEffTree* ze1 = new ZEffTree(*f1,true);
    //ZEffTree* ze2 = new ZEffTree(*f2,true);
    //
    //for(int event=0; event<ze0->Entries(); event++ ) //fill the data hists
    //{
    //    if(event%2)
    //    {
    //        copyZ(ze0, ze1);
    //        ze1->Fill();
    //    }
    //    else 
    //    {
    //        copyZ(ze0, ze2);
    //        ze2->Fill();
    //    }
    //    
    //    ze2->GetNextEvent();
    //}
    f0->Close();
    
    return 0;
}

//void copyZ(const ZEffTree* oldZ, ZEffTree* newZ)
//{
//    newZ->gen = oldZ->gen;
//    newZ->reco = oldZ->reco;
//    *(newZ->ixs) = *(oldZ->ixs);
//    *(newZ->iys) = *(oldZ->iys);
//    *(newZ->hitEnergyFractions) = *(oldZ->hitEnergyFractions);
//}