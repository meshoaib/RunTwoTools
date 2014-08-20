import FWCore.ParameterSet.Config as cms

#------ Setup ------#

#initialize the process
process = cms.Process("RazorJetAnalyzer")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("Configuration.EventContent.EventContent_cff")

#load input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'/store/mc/Spring14miniaod/SMS-T1tttt_2J_mGl-1500_mLSP-100_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/0A59FC95-330F-E411-94EB-E0CB4E29C4CA.root' #PU20, T1tttt
        ' /store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/003E832C-8AFC-E311-B7AA-002590596490.root' #PU20, TTBar
        #'/store/results/top/StoreResults/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/USER/Spring14dr_PU_S14_POSTLS170_V6AN1_miniAOD706p1_814812ec83fce2f620905d2bb30e9100-v2/00000/0012F41F-FA17-E411-A1FF-0025905A48B2.root' #PU40, TTBar
    )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

#TFileService for output (currently output of "ad-hoc" root files via CRAB3 is not supported)
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("razorJetAnalysis.root"),
    closeFileFast = cms.untracked.bool(True)
)

#load run conditions
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

#------ Declare the correct global tag ------#

#global tag for CSA14 25ns (asymptotic alignment and calibration) scenario
#process.GlobalTag.globaltag = 'PLS170_V7AN1::All'
#global tag for CSA14 50ns (more pessimistic alignment and calibration) scenario
process.GlobalTag.globaltag = 'PLS170_V6AN1::All'

#------ Analysis-specific imports ------#
process.load('RunTwoTools.MiniSelector.tools.ak4PFPuppiJetSequence')
process.load('RunTwoTools.MiniSelector.tools.ak8PFPuppiJetSequence')
from RecoJets.Configuration.RecoGenJets_cff import ak8GenJets
process.ak8GenJets = ak8GenJets.clone(src = 'packedGenParticles')

#------ Analyzer ------#

#list input collections
process.analyzer = cms.EDAnalyzer('RazorJetAnalyzer', 

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

    puppiJets = cms.InputTag("ak4PFPuppiJets"),
    puppiJetsAk8 = cms.InputTag("ak8PFPuppiJets"),
    fatGenJets = cms.InputTag("ak8GenJets")

)

#when running PAT, we have no choice but to use unscheduled execution
#process.options = cms.untracked.PSet(
#        allowUnscheduled = cms.untracked.bool(True)
#)

#run
process.p = cms.Path(
        process.ak4PFPuppiJetSequence *
        process.ak8PFPuppiJetSequence *
        process.ak8GenJets *
        process.analyzer)
