import FWCore.ParameterSet.Config as cms

#------ Setup ------#

#initialize the process
process = cms.Process("ExampleAnalyzer")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("Configuration.EventContent.EventContent_cff")

#load input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/cmst3/user/gpetrucc/miniAOD/v1/TT_Tune4C_13TeV-pythia8-tauola_PU_S14_PAT.root'
    )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

#TFileService for output (currently output of "ad-hoc" root files via CRAB3 is not supported)
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("example.root"),
    closeFileFast = cms.untracked.bool(True)
)

#load run conditions
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

#------ Declare the correct global tag ------#

#global tag for CSA14 25ns (asymptotic alignment and calibration) scenario
process.GlobalTag.globaltag = 'PLS170_V7AN1::All'
#global tag for CSA14 50ns (more pessimistic alignment and calibration) scenario
#process.GlobalTag.globaltag = 'PLS170_V6AN1::All'

#------ Analysis-specific imports ------#

process.load('RunTwoTools.MiniSelector.tools.ak5CHSJetSequence')

#add the raw PFMET producer (uncorrected MET not included by default in miniAOD)
from RecoMET.METProducers.PFMET_cfi import pfMet
process.pfMet = pfMet.clone(src = "packedPFCandidates")
process.pfMet.calculateSignificance = False # this can't be easily implemented on packed PF candidates at the moment

#------ Analyzer ------#

#list input collections
process.example = cms.EDAnalyzer('ExampleAnalyzer', 

    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    taus = cms.InputTag("slimmedTaus"),
    photons = cms.InputTag("slimmedPhotons"),
    jets = cms.InputTag("slimmedJets"),
    fatjets = cms.InputTag("slimmedJetsAK8"),
    mets = cms.InputTag("slimmedMETs"),
    pfCands = cms.InputTag("packedPFCandidates"),

    packed = cms.InputTag("packedGenParticles"),
    pruned = cms.InputTag("prunedGenParticles"),
    genjets = cms.InputTag("slimmedGenJets", "", "PAT"),

    bits = cms.InputTag("TriggerResults","","HLT"),
    prescales = cms.InputTag("patTrigger"),
    objects = cms.InputTag("selectedPatTrigger"),
    metBits = cms.InputTag("TriggerResults", "", "PAT"),

    rawPfMet = cms.InputTag("pfMet"),
    ak5PFJetsCHS = cms.InputTag("ak5PFJetsCHS")
)

#run
process.p = cms.Path(
        process.pfMet *
        process.ak5CHSJetSequence *
        process.example)
