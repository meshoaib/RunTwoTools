import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.ak5PFJets_cfi import ak5PFJets
from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets

# Select candidates that would pass CHS requirements
chs = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))

#makes chs ak5 jets   (instead of ak4 that are default in miniAOD 70X)
ak5PFJetsCHS = ak5PFJets.clone(src = 'chs')
#let's also make non-chs ak5 jets
ak5PFJets = ak5PFJets.clone(src = 'packedPFCandidates') 
ak5GenJets = ak5GenJets.clone(src = 'packedGenParticles')

ak5CHSJetSequence = cms.Sequence(chs * ak5PFJetsCHS * ak5PFJets * ak5GenJets)
