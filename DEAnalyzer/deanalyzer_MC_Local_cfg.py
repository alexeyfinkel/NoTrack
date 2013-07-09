import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/local/cms/phedex/store/mc/Summer11/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/00F71C3D-88A8-E011-AC48-001D0967D0DA.root',
        'file:/local/cms/phedex/store/mc/Summer11/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/0257F93D-85A8-E011-9244-0024E8768D34.root',
        'file:/local/cms/phedex/store/mc/Summer11/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/02E6BF8C-09A8-E011-AE29-00A0D1EE8D6C.root',
        'file:/local/cms/phedex/store/mc/Summer11/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/047818FC-A8A7-E011-BCFF-00151796C088.root',
        'file:/local/cms/phedex/store/mc/Summer11/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/04EB52AF-5BA8-E011-A4C2-00266CFABAF0.root',
        'file:/local/cms/phedex/store/mc/Summer11/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/065B5891-11A8-E011-987C-0015178C48FC.root',
        'file:/local/cms/phedex/store/mc/Summer11/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/083084E5-32A8-E011-B2CC-00A0D1EEF6B8.root',
        'file:/local/cms/phedex/store/mc/Summer11/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/08861A90-82A8-E011-8FAC-00151796D76C.root',
        'file:/local/cms/phedex/store/mc/Summer11/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/0A3F4864-A7A8-E011-A2D2-00151796C450.root',
        'file:/local/cms/phedex/store/mc/Summer11/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/0C3C9845-9AA8-E011-BBC7-00266CF9AB9C.root'
       
    )
)

process.DoubleElectron = cms.EDAnalyzer('DEAnalyzer',
	electronTag      = cms.InputTag( "gsfElectrons" ),
	photonTag        = cms.InputTag( "photons" )
)


process.p = cms.Path(process.DoubleElectron)

#--- Output histgram file ---#
process.TFileService = cms.Service("TFileService",
       fileName = cms.string("ToyDE_Local_001.root"),
)
