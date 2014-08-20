import FWCore.ParameterSet.Config as cms

from RunTwoTools.MiniSelector.tools.puppiPFCandidates import puppiPFCandidates
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets

#makes ak4 jets from the PUPPI particle collection
ak4PFPuppiJets = ak4PFJets.clone(src = 'puppiPFCandidates')

ak4PFPuppiJetSequence = cms.Sequence(puppiPFCandidates * ak4PFPuppiJets)
