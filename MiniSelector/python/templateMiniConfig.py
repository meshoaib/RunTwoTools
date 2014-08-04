#Config file for MiniSelector analyzer jobs

#   To use: copy this file, replace all instances of "MiniSelector"
#   with the name of your analyzer, enter the correct global tag info, 
#   and import any needed CMSSW configuration files.  
#   Add inputs in the "Analyzer" section, and add modules to the path as needed.

import FWCore.ParameterSet.Config as cms

#------ Setup ------#

#initialize the process
process = cms.Process("MiniSelector")
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
    fileName = cms.string("selectedInfo.root"),
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

#------ Analyzer ------#

#list input collections
process.template = cms.EDAnalyzer('MiniSelector', 

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

    #uncomment any of these lines
    #(make sure to add the corresponding EDM Handles to your analyzer)
    
    #genInfo = cms.InputTag("generator", "", "SIM"),
    #puInfo = cms.InputTag("addPileupInfo", "", "HLT"),
    #hcalNoiseInfo = cms.InputTag("hcalnoise", "", "RECO"),

    #secondaryVertices = cms.InputTag("slimmedSecondaryVertices", "", "PAT"),

    #rhoAll = cms.InputTag("fixedGridRhoAll", "", "RECO"),
    #rhoFastjetAll = cms.InputTag("fixedGridRhoFastjetAll", "", "RECO"),
    #rhoFastjetAllCalo = cms.InputTag("fixedGridRhoFastjetAllCalo", "", "RECO"),
    #rhoFastjetCentralCalo = cms.InputTag("fixedGridRhoFastjetCentralCalo", "", "RECO"),
    #rhoFastjetCentralChargedPileUp = cms.InputTag("fixedGridRhoFastjetCentralChargedPileUp", "", "RECO"),
    #rhoFastjetCentralNeutral = cms.InputTag("fixedGridRhoFastjetCentralNeutral", "", "RECO"),

    #beamSpot = cms.InputTag("offlineBeamSpot", "", "RECO"),

    #ebRecHits = cms.InputTag("reducedEgamma", "reducedEBRecHits", "PAT"),
    #eeRecHits = cms.InputTag("reducedEgamma", "reducedEERecHits", "PAT"),
    #esRecHits = cms.InputTag("reducedEgamma", "reducedESRecHits", "PAT"),
    #ebeeClusters = cms.InputTag("reducedEgamma", "reducedEBEEClusters", "PAT"),
    #esClusters = cms.InputTag("reducedEgamma", "reducedESClusters", "PAT"),
    #conversions = cms.InputTag("reducedEgamma", "reducedConversions", "PAT"),
    #singleLegConversions = cms.InputTag("reducedEgamma", "reducedSingleLegConversions", "PAT"),
    #gedGsfElectronCores = cms.InputTag("reducedEgamma", "reducedGedGsfElectronCores", "PAT"),
    #gedPhotonCores = cms.InputTag("reducedEgamma", "reducedGedPhotonCores", "PAT"),
    #superClusters = cms.InputTag("reducedEgamma", "reducedSuperClusters", "PAT"),

    #lostTracks = cms.InputTag("lostTracks", "", "PAT")
)

#run
process.p = cms.Path(
        process.template)
