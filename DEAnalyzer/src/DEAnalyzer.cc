// -*- C++ -*-
//
// Package:    DEAnalyzer
// Class:      DEAnalyzer
// 
/**\class DEAnalyzer DEAnalyzer.cc NoTrack/DEAnalyzer/src/DEAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Alexey Finkel
//         Created:  Mon Sep 10 18:37:03 CDT 2012
// $Id$
//
//


// system include files
#include <memory>
#include <iostream>
#include <algorithm>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
// Needed for 39X
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/FileBlock.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Math/interface/Point3D.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile2D.h"
#include "TVector3.h"
#include "TRandom.h"
#include "TLorentzVector.h"

//
// class declaration
//

class DEAnalyzer : public edm::EDAnalyzer {
   public:
      explicit DEAnalyzer(const edm::ParameterSet&);
      ~DEAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      
  	  int evtCounter;
  	  bool firstEvent_;
  	  
  	  edm::InputTag elecTag_; 
      edm::InputTag photTag_;
      
      struct DEevent
      {
      	double M;
      	double ElecPt;
      	double PhotPt;
      	double ElecEta;
      	double PhotEta;
      	double E1x5;
      	double E2x5;
      	double E3x3;
      	double E5x5;
      	double R9;
      	double HoEM;
      	double Sieie;
      	double ElEta;
      	double PhotonIso;
      	int No50;
      	DEevent():M(0), ElecPt(0), PhotPt(0), ElecEta(0), PhotEta(0), E1x5(0), E2x5(0), E3x3(0), E5x5(0), R9(0), HoEM(0), Sieie(0), ElEta(0), PhotonIso(0), No50(0) {}
      } ;

      // ----------member data ---------------------------
      
      class CutLevels
	  {
		public:		
			void book(TFileDirectory *, const std::string& );
			void fill(DEevent& dev);
			TH1 *EveryMass, *Nover50, *ElecPt, *ElecEta, *PhotPt, *PhotEta, *PhotE1x5, *PhotE2x5, *PhotE3x3, *PhotE5x5, *PhotR9, *PhotSieie, *PhotHoEM, *PhotonIso;
	  } ;
      
      TFileDirectory *rundir;
      
      TH1* ElectronEta;
      
      CutLevels PtEta, HoEM, Sieie, R9, EMxN, GSFelectron; 
      
      reco::Particle::LorentzVector V1;
};



void DEAnalyzer::CutLevels::book(TFileDirectory *mydir, const std::string& post)
{
	std::string t, T; // histogram title string;
	TH1::SetDefaultSumw2();
	t = post + "_AllM";
	T = post + " Highest Mass";
	EveryMass = mydir->make<TH1D>(t.c_str(), T.c_str(),100,0,200);
	t = post + "_PhotPt";
	T = post + " Photon Pt";
	PhotPt = mydir->make<TH1D>(t.c_str(), T.c_str(),100,0,100);
	t = post + "_ElecPt";
	T = post + " Electron Pt";
	ElecPt = mydir->make<TH1D>(t.c_str(), T.c_str(),100,0,100);
	t = post + "_PhotEta";
	T = post + " Photon Eta";
	PhotEta = mydir->make<TH1D>(t.c_str(), T.c_str(),100,-5,5);
	t = post + "_ElecEta";
	T = post + " Electron Eta";
	ElecEta = mydir->make<TH1D>(t.c_str(), T.c_str(),100,-5,5);
	t = post + "_PhotE1x5";
	T = post + " Photon E1x5";
	PhotE1x5 = mydir->make<TH1D>(t.c_str(), T.c_str(),100,0,1);
	t = post + "_PhotE2x5";
	T = post + " Photon E2x5";
	PhotE2x5 = mydir->make<TH1D>(t.c_str(), T.c_str(),100,0,1);
	t = post + "_PhotE3x3";
	T = post + " Photon E3x3";
	PhotE3x3 = mydir->make<TH1D>(t.c_str(), T.c_str(),100,0,1);
	t = post + "_PhotE5x5";
	T = post + " Photon E5x5";
	PhotE5x5 = mydir->make<TH1D>(t.c_str(), T.c_str(),100,0,1);
	t = post + "_PhotR9";
	T = post + " Photon R9";
	PhotR9 = mydir->make<TH1D>(t.c_str(), T.c_str(),100,0,1);
	t = post + "_PhotHoEM";
	T = post + " NT H/EM";
	PhotHoEM = mydir->make<TH1D>(t.c_str(), T.c_str(),50,0,0.5);
	t = post + "_PhotSiEiE";
	T = post + " Photon #sigma_{i#eta-i#eta}";
	PhotSieie = mydir->make<TH1D>(t.c_str(), T.c_str(),100,0,0.1);
	t = post + "_PhotonIso";
	T = post + " Photon Isolation";
	PhotonIso = mydir->make<TH1D>(t.c_str(), T.c_str(),100,0,1);
}

void DEAnalyzer::CutLevels::fill(DEevent& dev)
{
	EveryMass->Fill( dev.M );
	ElecPt->Fill( dev.ElecPt );
	ElecEta->Fill( dev.ElecEta );
	PhotPt->Fill( dev.PhotPt );
	PhotEta->Fill( dev.PhotEta );
	PhotHoEM->Fill( dev.HoEM );
	PhotR9->Fill( dev.R9 );   			
	PhotE1x5->Fill( dev.E1x5 );
	PhotE2x5->Fill( dev.E2x5 );
	PhotE3x3->Fill( dev.E3x3 );
	PhotE5x5->Fill( dev.E5x5 );
	PhotSieie->Fill( dev.Sieie );
	PhotonIso->Fill( dev.PhotonIso );
}



//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DEAnalyzer::DEAnalyzer(const edm::ParameterSet& iConfig)

{
   elecTag_ = iConfig.getParameter< edm::InputTag > ("electronTag");
   photTag_ = iConfig.getParameter< edm::InputTag > ("photonTag");

   // ==================== Book the histos ====================
   //
    
   edm::Service<TFileService> fs;
   
   PtEta.book(new TFileDirectory(fs->mkdir("PtEta")),"PtEta");
   HoEM.book(new TFileDirectory(fs->mkdir("HoEM")),"HoEM");
   Sieie.book(new TFileDirectory(fs->mkdir("Sieie")),"Sieie");
   R9.book(new TFileDirectory(fs->mkdir("R9")),"R9");
   EMxN.book(new TFileDirectory(fs->mkdir("EMxN")),"EMxN");
   GSFelectron.book(new TFileDirectory(fs->mkdir("GSFelectron")),"GSFelectron");
   
   ElectronEta = fs->make<TH1D>("ElEta","Electron Eta",35,3.5,3.5);
   
}


DEAnalyzer::~DEAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DEAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    
    DEevent dev;

    evtCounter++;

    edm::Handle<reco::GsfElectronCollection> recoElectrons;
    iEvent.getByLabel(elecTag_, recoElectrons);

    edm::Handle<reco::PhotonCollection> recoGammas;
    iEvent.getByLabel(photTag_, recoGammas);
    
    if(firstEvent_)
    {
        firstEvent_ = false;
        //std::cout<<"================================================================="<<std::endl;
    }

	if (!recoElectrons.isValid() || !recoGammas.isValid() ) return ;
	
	//int Nover = 0;
    //double Mass = 0;
    
    for (unsigned int i=0; i<recoElectrons->size(); i++)
    {
    	if( recoElectrons->at(i).p4().Pt() < 5 ) continue;    
    	if( recoElectrons->at(i).dr03EcalRecHitSumEt()/recoElectrons->at(i).p4().Et() >0.3 ) continue;
    	if( recoElectrons->at(i).hadronicOverEm() >0.2 ) continue;
    	if( recoElectrons->at(i).sigmaIetaIeta() >0.1 ) continue;
    	//if( recoElectrons->at(i).r9() < 0.7 ) continue;
    	ElectronEta->Fill( recoElectrons->at(i).p4().Eta() );
    	for (unsigned int j=0; j<recoGammas->size(); j++)
    	{
    		if ( fabs(recoGammas->at(j).p4().Eta()) < 2.5 || fabs(recoGammas->at(j).p4().Eta()) > 3 ) continue;
    		if( recoGammas->at(j).p4().Pt() < 10 ) continue;
    		V1 = recoElectrons->at(i).p4() + recoGammas->at(j).p4() ;
    		dev.M = V1.M();
    		dev.ElecPt = recoElectrons->at(i).p4().Pt();
    		dev.ElecEta = recoElectrons->at(i).p4().Eta();
    		dev.PhotPt = recoGammas->at(j).p4().Pt();
    		dev.PhotEta = recoGammas->at(j).p4().Eta();
    		dev.E1x5 = recoGammas->at(j).e1x5()/recoGammas->at(j).superCluster()->energy();
    		dev.E2x5 = recoGammas->at(j).e2x5()/recoGammas->at(j).superCluster()->energy();
    		dev.E3x3 = recoGammas->at(j).e3x3()/recoGammas->at(j).superCluster()->energy(); 
    		dev.E5x5 = recoGammas->at(j).e5x5()/recoGammas->at(j).superCluster()->energy();
    		dev.R9 = recoGammas->at(j).r9();
    		dev.HoEM = recoGammas->at(j).hadronicOverEm();
    		dev.Sieie = recoGammas->at(j).sigmaIetaIeta();
    		dev.PhotonIso = recoGammas->at(j).ecalRecHitSumEtConeDR03()/recoGammas->at(j).p4().Et();
    		
    		if( dev.M < 10 ) continue;
    		PtEta.fill(dev);
    		if( dev.HoEM > 0.05 ) continue;
    		HoEM.fill(dev);
    		if( dev.Sieie > 0.028 ) continue; 
    		Sieie.fill(dev);
    		if( dev.R9 > 0.99 || dev.R9 < 0.93 ) continue;
    		R9.fill(dev);
    		if( dev.E5x5 > 0.99 || dev.E2x5 > 0.96 || dev.E5x5 <0.86 || dev.E2x5 < 0.80 || dev.E1x5 > 0.9 )continue;
    		EMxN.fill(dev);	
    		
    	}
    }

	//now to get an idea of what those plots look like for GSF electrons
	
	for (unsigned int i=0; i<recoElectrons->size(); i++)
    {
    	if( recoElectrons->at(i).p4().Pt() < 5 ) continue;
    	if( recoElectrons->at(i).dr03EcalRecHitSumEt()/recoElectrons->at(i).p4().Et() >0.3 ) continue;
    	if( recoElectrons->at(i).r9() < 0.8 ) continue;
    	if( recoElectrons->at(i).hadronicOverEm() >0.06 ) continue;
    	if( recoElectrons->at(i).sigmaIetaIeta() >0.06 ) continue;
    	for (unsigned int j=0; j < recoElectrons->size(); j++)
    	{
    		if( j<=i ) continue;
    		if( recoElectrons->at(j).p4().Pt() < 7 ) continue;
    		if ( fabs(recoElectrons->at(j).p4().Eta()) < 2 || fabs(recoElectrons->at(j).p4().Eta()) > 2.5 ) continue;
    		if ((recoElectrons->at(i).charge()*recoElectrons->at(j).charge()) > 0) continue ;
    		V1 = recoElectrons->at(i).p4() + recoElectrons->at(j).p4() ;
    		dev.M = V1.M();
    		if( dev.M < 80 || dev.M >100 ) continue;
    		dev.ElecPt = recoElectrons->at(i).p4().Pt();
    		dev.ElecEta = recoElectrons->at(i).p4().Eta();
    		dev.PhotPt = recoElectrons->at(j).p4().Pt();
    		dev.PhotEta = recoElectrons->at(j).p4().Eta();
    		dev.E1x5 = recoElectrons->at(j).e1x5()/recoElectrons->at(j).superCluster()->energy();
    		dev.E2x5 = recoElectrons->at(j).e2x5Max()/recoElectrons->at(j).superCluster()->energy();
    		//dev.E3x3 = recoElectrons->at(j).e3x3()/recoElectrons->at(j).superCluster()->energy(); //GSF doesn't have this one...
    		dev.E5x5 = recoElectrons->at(j).e5x5()/recoElectrons->at(j).superCluster()->energy();
    		dev.R9 = recoElectrons->at(j).r9();
    		dev.HoEM = recoElectrons->at(j).hadronicOverEm();
    		dev.Sieie = recoElectrons->at(j).sigmaIetaIeta();
    		
    		GSFelectron.fill(dev);
    	}
    
    }
}


// ------------ method called once each job just before starting event loop  ------------
void 
DEAnalyzer::beginJob()
{
	firstEvent_ = true;
    evtCounter = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DEAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
DEAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
DEAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
DEAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
DEAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DEAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DEAnalyzer);
