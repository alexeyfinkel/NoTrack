import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring( )
)

process.DoubleElectron = cms.EDAnalyzer('DEAnalyzer',
	electronTag      = cms.InputTag( "gsfElectrons" ),
	photonTag        = cms.InputTag( "photons" )
)


process.p = cms.Path(process.DoubleElectron)

#--- Output histgram file ---#
process.TFileService = cms.Service("TFileService",
       fileName = cms.string("SingleElectronTrigger_Local_001.root"),
)
