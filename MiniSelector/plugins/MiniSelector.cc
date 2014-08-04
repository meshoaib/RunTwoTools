// -*- C++ -*-
// Class:      MiniSelector
/*
Description: Base class for miniAOD analysis with CRAB
*/
// Original Author:  Dustin James Anderson
//         Created:  Thu, 17 Jul 2014 15:00:06 GMT

// ********************************************************************************//
// ********************************************************************************//
// ------ PLACE USEFUL HELPER FUNCTIONS IN THIS FILE FOR EVERYONE TO USE --------- //
// ********************************************************************************//
// ********************************************************************************//

#include "MiniSelector.h"

vector<TLorentzVector> MiniSelector::getHemispheres(vector<TLorentzVector> jets){
    int nJets = jets.size();
    vector<TLorentzVector> possibleHem1s; //holds possible hemisphere combinations
    vector<TLorentzVector> possibleHem2s;

    //stolen from https://github.com/pierinim/BSMatLHC/blob/master/BSMApp/src/CMS/CMSHemisphere.cc
    int nComb = pow(2, nJets);

    //step 1: store all possible partitions of the input jets
    int j_count;
    //this loop looks like witchcraft but it basically just counts in binary
    for(int i = 1; i < nComb-1; i++){ //note we omit the trivial hemisphere combinations (0 and nComb-1)
        TLorentzVector j_temp1, j_temp2;
        int itemp = i;
        j_count = nComb/2;
        int count = 0;
        while(j_count > 0){ //decompose i into binary '1's and '0's ; put the '1' jets into j_temp1 and the '0' jets into j_temp2
            if(itemp/j_count == 1){
                j_temp1 += jets[count];
            } else {
                j_temp2 += jets[count];
            }
            itemp -= j_count*(itemp/j_count); //note this is always (0 or 1)*j_count
            j_count /= 2;
            count++;
        }
        possibleHem1s.push_back(j_temp1);
        possibleHem2s.push_back(j_temp2);
    }
 
    //step 2: choose the partition that minimizes m1^2 + m2^2
    double mMin = -1;
    TLorentzVector myHem1;
    TLorentzVector myHem2;
    for(size_t i=0; i < possibleHem1s.size(); i++){
        double mTemp = possibleHem1s[i].M2() + possibleHem2s[i].M2();
        if(mMin < 0 || mTemp < mMin){
            mMin = mTemp;
            myHem1 = possibleHem1s[i];
            myHem2 = possibleHem2s[i];
        }
    }

    //return the hemispheres in decreasing order of pt
    vector<TLorentzVector> hemsOut;
    if(myHem1.Pt() > myHem2.Pt()){
        hemsOut.push_back(myHem1);
        hemsOut.push_back(myHem2);
    } else {
        hemsOut.push_back(myHem2);
        hemsOut.push_back(myHem1);
    }

    return hemsOut;
}

double MiniSelector::computeMR(TLorentzVector hem1, TLorentzVector hem2){
    return sqrt(pow(hem1.P() + hem2.P(), 2) - pow(hem1.Pz() + hem2.Pz(), 2));
}

double MiniSelector::computeR2(TLorentzVector hem1, TLorentzVector hem2, TLorentzVector pfMet){
    double mR = computeMR(hem1, hem2);
    double term1 = pfMet.Pt()/2*(hem1.Pt() + hem2.Pt());
    double term2 = pfMet.Px()/2*(hem1.Px() + hem2.Px()) + pfMet.Py()/2*(hem1.Py() + hem2.Py()); //dot product of MET with (p1T + p2T)
    double mTR = sqrt(term1 - term2);
    return (mTR / mR) * (mTR / mR);
}

bool MiniSelector::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle){
    //particle is already the ancestor
    if(ancestor == particle ) return true;

    //otherwise loop on mothers, if any and return true if the ancestor is found
    for(size_t i=0;i< particle->numberOfMothers();i++)
    {
        if(isAncestor(ancestor,particle->mother(i))) return true;
    }
    //if we did not return yet, then particle and ancestor are not relatives
    return false;
}

// constructors and destructor
MiniSelector::MiniSelector(const edm::ParameterSet& iConfig): 
    //get inputs from config file
    vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
    electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
    tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
    photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
    jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
    fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
    pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
    prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
    packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed"))),
    genJetToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genjets"))),
    triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
    triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
    triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),     
    metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
    metFilterBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("metBits")))
{
}

MiniSelector::~MiniSelector()
{
}

// member functions

void MiniSelector::loadEvent(const edm::Event& iEvent){
    iEvent.getByToken(triggerBits_, triggerBits);
    iEvent.getByToken(triggerObjects_, triggerObjects);
    iEvent.getByToken(triggerPrescales_, triggerPrescales);
    iEvent.getByToken(metFilterBits_, metFilterBits);
    iEvent.getByToken(vtxToken_, vertices);
    iEvent.getByToken(pfToken_, pfs);
    iEvent.getByToken(muonToken_, muons);
    iEvent.getByToken(electronToken_, electrons);
    iEvent.getByToken(photonToken_, photons);
    iEvent.getByToken(tauToken_, taus);
    iEvent.getByToken(jetToken_, jets);
    iEvent.getByToken(fatjetToken_, fatjets);
    iEvent.getByToken(metToken_, mets);
    iEvent.getByToken(prunedGenToken_,pruned);
    iEvent.getByToken(packedGenToken_,packed);
    iEvent.getByToken(genJetToken_,genjets);
}   

//------ Method called for each event ------//

void MiniSelector::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    using namespace edm;

    //initialize
    loadEvent(iEvent);

}


//------ Method called once each job just before starting event loop ------//
void MiniSelector::beginJob(){
}

//------ Method called once each job just after ending the event loop ------//
void MiniSelector::endJob(){
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniSelector);
