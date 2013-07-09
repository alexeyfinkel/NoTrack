import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_000.root',
       
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
