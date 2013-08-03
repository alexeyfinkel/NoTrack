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
    fileNames=cms.untracked.vstring(
         '/store/data/Run2012A/DoubleElectron/AOD/PromptReco-v1/000/190/641/CE221EB2-6A82-E111-944A-001D09F23C73.root',
         '/store/data/Run2012A/DoubleElectron/AOD/PromptReco-v1/000/190/595/32C3878A-2982-E111-97DC-003048F1183E.root',
         '/store/data/Run2012A/DoubleElectron/AOD/PromptReco-v1/000/190/595/34229116-1F82-E111-8CC5-0030486780E6.root',
         '/store/data/Run2012A/DoubleElectron/AOD/PromptReco-v1/000/190/595/4268FD63-2582-E111-B6D9-002481E0D7D8.root',
         '/store/data/Run2012A/DoubleElectron/AOD/PromptReco-v1/000/190/595/56DC73FB-2582-E111-9FE0-001D09F29114.root',
         '/store/data/Run2012A/DoubleElectron/AOD/PromptReco-v1/000/190/595/8245B417-3982-E111-9636-003048F118C6.root',
         '/store/data/Run2012A/DoubleElectron/AOD/PromptReco-v1/000/190/595/AAB48925-2882-E111-A6C0-001D09F24664.root',
         '/store/data/Run2012A/DoubleElectron/AOD/PromptReco-v1/000/190/595/C22AA016-1F82-E111-B790-003048D2BB90.root',
         '/store/data/Run2012A/DoubleElectron/AOD/PromptReco-v1/000/190/628/5AF5EE68-2582-E111-896D-001D09F2932B.root',
         '/store/data/Run2012A/DoubleElectron/AOD/PromptReco-v1/000/190/482/1A507AD9-1981-E111-9018-003048F11C5C.root'
     )
        )

process.out = cms.OutputModule( "PoolOutputModule",
       fileName = cms.untracked.string("/local/cms/user/gude/2012_single_electron_skim//2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_001.root"),
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
        throw                   = cms.bool(True),                   # throw exception on unknown trigger names
        triggerConditions       = cms.vstring( 
            'HLT_Ele17_CaloIdL_CaloIsoVL_v*',
            'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*',
            'HLT_Ele8_CaloIdL_CaloIsoVL_v*',
            'HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*',
            'HLT_Ele8_CaloIdT_TrkIdVL_v*'
            )
        )

process.p1 = cms.Path( 
        process.hltPickTriggered
        )

process.final = cms.EndPath(
        process.out
        )
