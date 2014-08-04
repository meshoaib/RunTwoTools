// -*- C++ -*-
// Class:      MiniSelector
/*
Description: Base class for miniAOD analysis with CRAB
*/
// Original Author:  Dustin James Anderson
//         Created:  Thu, 17 Jul 2014 15:00:06 GMT

#ifndef MINISELECTOR_H
#define MINISELECTOR_H

// system include files
#include <memory>
#include <string>
#include <vector>

using namespace std;

// CMSSW framework includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//CMSSW package includes
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

//ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

//------ Class declaration ------//

class MiniSelector : public edm::EDAnalyzer {
    public:
        //analyzer constructor and destructor
        explicit MiniSelector(const edm::ParameterSet&);
        ~MiniSelector();

        //------ HELPER FUNCTIONS ------//
        
        //splits jets into two hemisperes for razor variable calculation
        //(minimizes sum of mass^2's of hemispheres)
        vector<TLorentzVector> getHemispheres(vector<TLorentzVector> jets);
        //compute M_R using two hemispheres
        double computeMR(TLorentzVector hem1, TLorentzVector hem2);
        //compute R^2 using two hemispheres and MET vector
        double computeR2(TLorentzVector hem1, TLorentzVector hem2, TLorentzVector pfMet);
        //returns true if particle 1 is an ancestor of particle 2, false otherwise
        //(takes two members of prunedGenParticles)
        bool isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle);

        void loadEvent(const edm::Event& iEvent); //call at the beginning of each event to get input handles from the python config

    protected:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        //----- Member data ------//

        //primary vertices
        edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
        //particle info
        edm::EDGetTokenT<pat::MuonCollection> muonToken_;
        edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
        edm::EDGetTokenT<pat::TauCollection> tauToken_;
        edm::EDGetTokenT<pat::PhotonCollection> photonToken_;
        edm::EDGetTokenT<pat::JetCollection> jetToken_;
        edm::EDGetTokenT<pat::JetCollection> fatjetToken_;
        edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
        //MC particle info
        edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
        edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
        edm::EDGetTokenT<reco::GenJetCollection> genJetToken_;
        //trigger info
        edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
        edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
        edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
        //MET info
        edm::EDGetTokenT<pat::METCollection> metToken_;
        edm::EDGetTokenT<edm::TriggerResults> metFilterBits_;

        //EDM handles for each input object
        edm::Handle<edm::TriggerResults> triggerBits;
        edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
        edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
        edm::Handle<edm::TriggerResults> metFilterBits;
        edm::Handle<reco::VertexCollection> vertices;
        edm::Handle<pat::PackedCandidateCollection> pfs;
        edm::Handle<pat::MuonCollection> muons;
        edm::Handle<pat::ElectronCollection> electrons;
        edm::Handle<pat::PhotonCollection> photons;
        edm::Handle<pat::TauCollection> taus;
        edm::Handle<pat::JetCollection> jets;
        edm::Handle<pat::JetCollection> fatjets;
        edm::Handle<pat::METCollection> mets;
        edm::Handle<edm::View<reco::GenParticle> > pruned;
        edm::Handle<edm::View<pat::PackedGenParticle> > packed;
        edm::Handle<reco::GenJetCollection> genjets;

};

#endif
