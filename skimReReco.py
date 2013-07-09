import FWCore.ParameterSet.Config as cms

process = cms.Process("Work2")
process.load("RecoEgamma.EgammaHFProducers.hfEMClusteringSequence_cff")

#process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cff")

#process.load("Geometry.CaloEventSetup.CaloGeometry_cff")

#process.load("Configuration.EventContent.EventContent_cff")
process.load('FWCore/MessageService/MessageLogger_cfi')

## Load additional processes
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
## global tags:
#process.GlobalTag.globaltag = cms.string('GR_R_42_V19::All')
process.GlobalTag.globaltag = cms.string('GR_R_52_V8::All')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Services_cff')


process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options   = cms.untracked.PSet( 
        wantSummary = cms.untracked.bool(True) 
        )

process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
        )

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring( 'file:/local/cms/phedex/store/data/Run2012A/DoubleElectron/AOD/22Jan2013-v1/20000/003EC246-5E67-E211-B103-00259059642E.root')
        )

process.out = cms.OutputModule( "PoolOutputModule",
        fileName = cms.untracked.string("/local/cms/user/finkel/NoTrack/Skim/SeSkimTest_01.root"),
        SelectEvents = cms.untracked.PSet(
            SelectEvents = cms.vstring('p1'),
            ),
        outputCommands = cms.untracked.vstring(
            'keep *',
            )
        )

## Skimmer
process.hltPickTriggered = cms.EDFilter('TriggerResultsFilter',
        hltResults              = cms.InputTag('TriggerResults','','HLT'),   # HLT results   - set to empty to ignore HLT
        l1tResults              = cms.InputTag(''),                 # L1 GT results - set to empty to ignore L1
        l1tIgnoreMask           = cms.bool(False),                  # use L1 mask
        l1techIgnorePrescales   = cms.bool(False),                  # read L1 technical bits from PSB#9, bypassing the prescales
        daqPartitions           = cms.uint32(0x01),                 # used by the definition of the L1 mask
        throw                   = cms.bool(False),                  # throw exception on unknown trigger names
        triggerConditions       = cms.vstring( 
            'HLT_Ele27_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele15_CaloIdT_CaloIsoVL_trackless_v*',
            'HLT_Ele27_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_HFT15_v*',
            'HLT_Ele23_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_HFT30_v*',
            'HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v*'
            )
        )

process.p1 = cms.Path( 
        process.hltPickTriggered
        )

process.final = cms.EndPath(
        process.out
        )
