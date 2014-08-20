import FWCore.ParameterSet.Config as cms

from RunTwoTools.MiniSelector.tools.puppiPFCandidates import puppiPFCandidates
from RecoJets.Configuration.RecoPFJets_cff import ak8PFJets

#makes ak8 jets from the PUPPI particle collection
ak8PFPuppiJets = ak8PFJets.clone(src = 'puppiPFCandidates')

ak8PFPuppiJetSequence = cms.Sequence(puppiPFCandidates * ak8PFPuppiJets)

