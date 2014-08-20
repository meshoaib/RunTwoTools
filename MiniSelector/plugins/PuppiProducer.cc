#include <string>

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Common/interface/Association.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

#include "miniAODTest/METAnalyzer/plugins/puppiCleanContainer.hh"

typedef math::XYZTLorentzVector LorentzVector;

class PuppiProducer : public edm::EDProducer {
    public:
        explicit PuppiProducer(const edm::ParameterSet&);
        ~PuppiProducer();

        virtual void produce(edm::Event&, const edm::EventSetup&);

    private:
        edm::EDGetTokenT<pat::PackedCandidateCollection> Cands_;
};

PuppiProducer::PuppiProducer(const edm::ParameterSet& iConfig) :
  Cands_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("src")))
{
  produces< std::vector<pat::PackedCandidate> > ();
}

PuppiProducer::~PuppiProducer() {}

void PuppiProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<pat::PackedCandidateCollection> cands;
    iEvent.getByToken( Cands_, cands );

    //Get PUPPI rescaled momenta
    puppiCleanContainer puppiBuilder(*cands, 2.5, false, true); //the second boolean determines whether tuning is used
    std::vector<fastjet::PseudoJet> puppiParticles = puppiBuilder.puppiEvent(7, 0.5);//parameters taken from code found at https://twiki.cern.ch/twiki/bin/viewauth/CMS/PUPPI

    //vector for output
    std::auto_ptr< std::vector<pat::PackedCandidate> > outPtrP( new std::vector<pat::PackedCandidate> );

    //rescale PackedCandidate momenta by PUPPI weights, throwing out particles with weight 0
    for(unsigned int ic=0, nc = cands->size(); ic < nc; ++ic) {
        const pat::PackedCandidate &candDummy=(*cands)[ic];
        auto_ptr<pat::PackedCandidate> candPtr(new pat::PackedCandidate(candDummy));
        pat::PackedCandidate &cand = (*candPtr);

        if(puppiParticles[ic].e() == 0) continue; //weight 0 particle 
        LorentzVector newP;
        newP.SetXYZT(puppiParticles[ic].px(), puppiParticles[ic].py(), puppiParticles[ic].pz(), puppiParticles[ic].e());
        cand.setP4(newP);
        
        outPtrP->push_back(pat::PackedCandidate(cand));
    }

    iEvent.put( outPtrP );
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PuppiProducer);
