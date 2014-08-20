import FWCore.ParameterSet.Config as cms

puppiPFCandidates = cms.EDProducer("PuppiProducer",
    src      = cms.InputTag("packedPFCandidates")
)
