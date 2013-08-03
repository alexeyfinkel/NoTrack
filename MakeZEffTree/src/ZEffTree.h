#include "TFile.h"
#include "TBranch.h"
#include "TTree.h"
#include <map>
#include <string>
//#include "DataFormats/DetId/interface/DetId.h"
//#include "DataFormats/EcalDetId/interface/EEDetId.h"

#ifndef ZEffTree_h_included
#define ZEffTree_h_included


class ZEffTree {
    public:

        ZEffTree(TFile& f, bool writable);

        static int cutToBit(const std::string& cut);

        struct ZInfo {
            float eta[2];
            float phi[2];
            float pt[2];
            float mz, yz, qtz;
            float Eiso[2];
            float Hiso[2];
            float HoEM[2];
            float Sieie[2];
            float RNine[2];
            //float E1x5[2];
      		//float E2x5[2];
      		//float E3x3[2];
      		//float E5x5[2];
      		int bits[2];
            int nverts;
            bool isSelected(int ielec, const std::string& bitname) const;
            void setBit(int ielec, const std::string& bitname, bool val);
            int charge[2];
            int ix[2], iy[2];
            float phistar;
            
            //std::vector<std::pair<EEDetId,float> > *hitsAndFractions;
            
        } gen, reco;
        
        //std::vector<int> *ixs, *iys;
        //std::vector<float> *hitEnergyFractions;

	std::vector<int> *ixs; 
	std::vector<int> *iys;
	std::vector<float> *hitEnergyFractions;

        void Fill() { 
            m_tree->Fill();
        }

        void Write() { 
            m_file.Write();
        }

        int Entries() { 
            return m_tree->GetEntries(); 
        }

        bool GetNextEvent() {
            if (nevt==Entries()) {
                return false; 
            }
            m_tree->GetEntry(nevt); 
            nevt++; 
            return true; 
        }

        void Clear() { 
            gen.eta[0] = 12; gen.eta[1] = 12; gen.phi[0] = 10; gen.phi[1] = 10; gen.pt[0] = -10; gen.pt[1] = -10; gen.mz = -10; gen.yz = -10; gen.qtz = 0; gen.bits[0] = 0; gen.bits[1] = 0; gen.nverts = -10; gen.charge[0] = 5; gen.charge[1] = 5; gen.ix[0]=0; gen.ix[1]=0; gen.iy[0]=0; gen.iy[1]=0; gen.phistar = -999.; // gen.Eiso[0] = -10; gen.Eiso[1] = -10; gen.Hiso[0] = -10; gen.Hiso[1] = -10; gen.R9[0] = -10; gen.R9[1] = -10; gen.HoEM[0] = -10; gen.HoEM[1] = -10; gen.Sieie[0] = -10; gen.Sieie[1] = -10; gen.E1x5[0] = -10; gen.E1x5[1] = -10; gen.E2x5[0] = -10; gen.E2x5[1] = -10; gen.E3x3[0] = -10; gen.E3x3[1] = -10; gen.E5x5[0] = -10; gen.E5x5[1] = -10;  
            
            reco.eta[0] = 10; reco.eta[1] = 10; reco.phi[0] = 10; reco.phi[1] = 10; reco.pt[0] = 10; reco.pt[1] = -10; reco.mz = -10; reco.yz = -10; reco.qtz = 0; reco.bits[0] = 0; reco.bits[1] = 0; reco.nverts = -10; reco.charge[0] = 5; reco.charge[1] = 5; reco.ix[0]=-10; reco.ix[1]=-10; reco.iy[0]=-10; reco.iy[1]=-10; reco.phistar = -999.; reco.Eiso[0] = 10; reco.Eiso[1] = 10; reco.Hiso[0] = 10; reco.Hiso[1] = 10; reco.RNine[0] = -10; reco.RNine[1] = -10; reco.HoEM[0] = 10; reco.HoEM[1] = 10; reco.Sieie[0] = 10; reco.Sieie[1] = 10;
            ixs->clear(); iys->clear(); hitEnergyFractions->clear(); // reco.hitsAndFractions->clear(); reco.E1x5[0] = -10; reco.E1x5[1] = -10; reco.E2x5[0] = -10; reco.E2x5[1] = -10; reco.E3x3[0] = -10; reco.E3x3[1] = -10; reco.E5x5[0] = -10; reco.E5x5[1] = -10; 
        }

    private:
        TFile& m_file;
        TTree* m_tree;
        int nevt;
        TBranch* br_gen, *br_reco, *br_ixs, *br_iys, *br_fracs;

        void makeBranches(bool writable) {
            nevt=0;
	    ixs = new std::vector<int>();
	    iys = new std::vector<int>();
	    hitEnergyFractions = new std::vector<float>();

            if (writable) {
                m_file.cd();
                m_tree=new TTree("ZEffs","Minnesota ZEffs");
                br_gen=m_tree->Branch("gen",&gen,"eta0/f:eta1:phi0:phi1:pt0:pt1:mz:yz:qtz:bits0/I:bits1:nverts:charge0:charge1:ix0:ix1:iy0:iy1:phistar/f");
                br_reco=m_tree->Branch("reco",&reco,"eta0/f:eta1:phi0:phi1:pt0:pt1:mz:yz:qtz:Eiso0:Eiso1:Hiso0:Hiso1:RNine0:RNine1:HoEM0:HoEM1:Sieie0:Sieie1:bits0/I:bits1:nverts:charge0:charge1:ix0:ix1:iy0:iy1:phistar/f");
                br_fracs=m_tree->Branch("hitEnergyFractions", hitEnergyFractions);
                br_ixs=m_tree->Branch("ixs", ixs);
                br_iys=m_tree->Branch("iys", iys);
            } else {
                m_tree=(TTree*)m_file.Get("ZEffs");
                m_tree->SetBranchAddress("gen",&gen);
                m_tree->SetBranchAddress("reco",&reco);
                m_tree->SetBranchAddress("ixs", &ixs);//, &br_ixs
                m_tree->SetBranchAddress("iys", &iys);//, &br_iys
                m_tree->SetBranchAddress("hitEnergyFractions", &hitEnergyFractions);//, &br_fracs
            }
        }

        /* Map from names to bitnums  */
        static std::map<std::string, int> cutToBits_;

        void prepBitmap();

};

#endif
