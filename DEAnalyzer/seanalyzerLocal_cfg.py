import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.options = cms.untracked.PSet(
wantSummary = cms.untracked.bool(True)
)

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_000.root',
        'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_001.root',
        'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_002.root',
        'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_003.root',
        'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_004.root',
        'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_005.root',
        'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_006.root',
        'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_007.root',
        'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_008.root',
        'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_009.root',
        #'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_010.root',
        #'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_011.root',
        #'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_012.root',
        #'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_013.root',
        #'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_014.root',
        #'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_015.root',
        #'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_016.root',
        #'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_017.root',
        #'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_018.root',
        #'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_019.root',
        #'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_020.root',
        #'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_021.root',
        #'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_022.root',
        #'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_023.root',
        #'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_024.root',
        #'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_025.root',
        #'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_026.root',
        #'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_027.root',
        #'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_028.root',
        #'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_029.root',
        #'file:/local/cms/user/gude/2012_single_electron_skim/2012A_Single_Electron_Trigger_Skim/2012A_Single_Electron_Trigger_Skim_030.root'
    )
)

process.DoubleElectron = cms.EDAnalyzer('DEAnalyzer',
	electronTag      = cms.InputTag( "gsfElectrons" ),
	photonTag        = cms.InputTag( "photons" )
)


process.p = cms.Path(process.DoubleElectron)

#--- Output histgram file ---#
process.TFileService = cms.Service("TFileService",
       fileName = cms.string("SingleElectronTrigger_Local_002.root"),
)
