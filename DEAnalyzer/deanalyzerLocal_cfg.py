import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_000.root',
        'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_001.root'
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_002.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_003.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_004.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_005.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_006.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_007.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_008.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_009.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_010.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_011.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_012.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_013.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_014.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_015.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_016.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_017.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_018.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_019.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_020.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_021.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_022.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_023.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_024.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_025.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_026.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_027.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_028.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_029.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_030.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_031.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_032.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_033.root',
        #'file:/local/cms/user/finkel/NoTrack/Skim/skimReReco/skimReReco_034.root'
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
