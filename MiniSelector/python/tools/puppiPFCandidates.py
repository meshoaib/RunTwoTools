import FWCore.ParameterSet.Config as cms

#the new PUPPI implementation is found on github, https://github.com/violatingcp/Dummy
from Dummy.Puppi.Puppi_cff import puppi
puppiPFCandidates = puppi.clone()
