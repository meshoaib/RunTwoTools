// Class:      RazorJetAnalyzer
// Description: template for a class extending MiniSelector

// Author:  Dustin James Anderson

//------ Class declaration ------//

#include "MiniSelector.h"

//analysis-specific includes
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

using namespace std;

class RazorJetAnalyzer : public MiniSelector {
    public:
        //analyzer constructor and destructor
        explicit RazorJetAnalyzer(const edm::ParameterSet&);
        ~RazorJetAnalyzer();

        //helper functions
        void loadEvent(const edm::Event& iEvent); //call at the start of each event

    private:
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

        //----- Member data ------//

        //non-default inputs 
        edm::EDGetTokenT<reco::PFJetCollection> puppiJetToken_;
        edm::EDGetTokenT<reco::PFJetCollection> puppiJetAk8Token_;
        edm::EDGetTokenT<reco::GenJetCollection> genJetAk8Token_;
        edm::Handle<reco::PFJetCollection> puppiJets;
        edm::Handle<reco::PFJetCollection> puppiJets8;
        edm::Handle<reco::GenJetCollection> genJets8;

        //output tree
        TTree *outputTree;

        //------ Variables for tree ------//
        
        //muon variables
        int nMuons;
        vector<double> muPt, muEta, muPhi;

        //electron variables
        int nElectrons;
        vector<double> elePt, eleEta, elePhi;
        
        //jet variables
        int nJets4;
        vector<double> jet4Pt, jet4Eta, jet4Phi, jet4Mass;

        //fatjet variables
        int nJets8;
        vector<double> jet8Pt, jet8Eta, jet8Phi, jet8Mass;

        //genjet variables
        int nGenJets4;
        vector<double> genJet4Pt, genJet4Eta, genJet4Phi, genJet4Mass;
        int nGenJets8;
        vector<double> genJet8Pt, genJet8Eta, genJet8Phi, genJet8Mass;

        //MET variables
        double met;
        double phiMet;
        double genMet;

        //PUPPI variables
        int nPuppiJets4;
        vector<double> puppiJet4Pt, puppiJet4Eta, puppiJet4Phi, puppiJet4Mass;
        int nPuppiJets8;
        vector<double> puppiJet8Pt, puppiJet8Eta, puppiJet8Phi, puppiJet8Mass;

        //razor variables
        double MR4, RSquared4;
        double MRPuppi, RSquaredPuppi;
        double MRGen4, RSquaredGen4;

        double MR8, RSquared8;
        double MRPuppi8, RSquaredPuppi8;
        double MRGen8, RSquaredGen8;

        //vertex
        int nPV;
};

// constructors and destructor
RazorJetAnalyzer::RazorJetAnalyzer(const edm::ParameterSet& iConfig): MiniSelector(iConfig),
    puppiJetToken_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("puppiJets"))),
    puppiJetAk8Token_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("puppiJetsAk8"))),  
    genJetAk8Token_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("fatGenJets")))
{
    //declare the TFileService for output
    edm::Service<TFileService> fs;
    //set up output tree
    outputTree = fs->make<TTree>("outputTree", "selected miniAOD information");

    //initialize tree branches
    outputTree->Branch("nJets4", &nJets4, "nJets4/I");
    outputTree->Branch("jet4Pt", "std::vector<double>", &jet4Pt);
    outputTree->Branch("jet4Eta", "std::vector<double>", &jet4Eta);
    outputTree->Branch("jet4Phi", "std::vector<double>", &jet4Phi);
    outputTree->Branch("jet4Mass", "std::vector<double>", &jet4Mass);

    outputTree->Branch("nJets8", &nJets8, "nJets8/I");
    outputTree->Branch("jet8Pt", "std::vector<double>", &jet8Pt);
    outputTree->Branch("jet8Eta", "std::vector<double>", &jet8Eta);
    outputTree->Branch("jet8Phi", "std::vector<double>", &jet8Phi);
    outputTree->Branch("jet8Mass", "std::vector<double>", &jet8Mass);

    outputTree->Branch("nPuppiJets4", &nPuppiJets4, "nPuppiJets4/I");
    outputTree->Branch("puppiJet4Pt", "std::vector<double>", &puppiJet4Pt);
    outputTree->Branch("puppiJet4Eta", "std::vector<double>", &puppiJet4Eta);
    outputTree->Branch("puppiJet4Phi", "std::vector<double>", &puppiJet4Phi);
    outputTree->Branch("puppiJet4Mass", "std::vector<double>", &puppiJet4Mass);

    outputTree->Branch("nPuppiJets8", &nPuppiJets8, "nPuppiJets8/I");
    outputTree->Branch("puppiJet8Pt", "std::vector<double>", &puppiJet8Pt);
    outputTree->Branch("puppiJet8Eta", "std::vector<double>", &puppiJet8Eta);
    outputTree->Branch("puppiJet8Phi", "std::vector<double>", &puppiJet8Phi);
    outputTree->Branch("puppiJet8Mass", "std::vector<double>", &puppiJet8Mass);

    outputTree->Branch("nGenJets4", &nGenJets4, "nGenJets4/I");
    outputTree->Branch("genJet4Pt", "std::vector<double>", &genJet4Pt);
    outputTree->Branch("genJet4Eta", "std::vector<double>", &genJet4Eta);
    outputTree->Branch("genJet4Phi", "std::vector<double>", &genJet4Phi);
    outputTree->Branch("genJet4Mass", "std::vector<double>", &genJet4Mass);
    
    outputTree->Branch("nGenJets8", &nGenJets8, "nGenJets8/I");
    outputTree->Branch("genJet8Pt", "std::vector<double>", &genJet8Pt);
    outputTree->Branch("genJet8Eta", "std::vector<double>", &genJet8Eta);
    outputTree->Branch("genJet8Phi", "std::vector<double>", &genJet8Phi);
    outputTree->Branch("genJet8Mass", "std::vector<double>", &genJet8Mass);
    
    outputTree->Branch("nMuons", &nMuons, "nMuons/I");
    outputTree->Branch("muPt", "std::vector<double>", &muPt);
    outputTree->Branch("muEta", "std::vector<double>", &muEta);
    outputTree->Branch("muPhi", "std::vector<double>", &muPhi);
    
    outputTree->Branch("nElectrons", &nElectrons, "nElectrons/I");
    outputTree->Branch("elePt", "std::vector<double>", &elePt);
    outputTree->Branch("eleEta", "std::vector<double>", &eleEta);
    outputTree->Branch("elePhi", "std::vector<double>", &elePhi);

    outputTree->Branch("met", &met, "met/D");
    outputTree->Branch("phiMet", &phiMet, "phiMet/D");
    outputTree->Branch("genMet", &genMet, "genMet/D");

    outputTree->Branch("MR4", &MR4, "MR4/D");
    outputTree->Branch("RSquared4", &RSquared4, "RSquared4/D");
    outputTree->Branch("MRGen4", &MRGen4, "MRGen4/D");
    outputTree->Branch("RSquaredGen4", &RSquaredGen4, "RSquaredGen4/D");
    outputTree->Branch("MRPuppi", &MRPuppi, "MRPuppi/D");
    outputTree->Branch("RSquaredPuppi", &RSquaredPuppi, "RSquaredPuppi/D");

    outputTree->Branch("MR8", &MR8, "MR8/D");
    outputTree->Branch("RSquared8", &RSquared8, "RSquared8/D");
    outputTree->Branch("MRGen8", &MRGen8, "MRGen8/D");
    outputTree->Branch("RSquaredGen8", &RSquaredGen8, "RSquaredGen8/D");
    outputTree->Branch("MRPuppi8", &MRPuppi8, "MRPuppi8/D");
    outputTree->Branch("RSquaredPuppi8", &RSquaredPuppi8, "RSquaredPuppi8/D");

    outputTree->Branch("nPV", &nPV, "nPV/I");

}


RazorJetAnalyzer::~RazorJetAnalyzer()
{
}


// member functions

void RazorJetAnalyzer::loadEvent(const edm::Event& iEvent){

    //load default inputs
    MiniSelector::loadEvent(iEvent);
    
    //read extra inputs
    iEvent.getByToken(puppiJetToken_, puppiJets);
    iEvent.getByToken(puppiJetAk8Token_, puppiJets8);
    iEvent.getByToken(genJetAk8Token_, genJets8);
    
    //reset tree variables
    nJets4 = 0;
    jet4Pt.clear();
    jet4Eta.clear();
    jet4Phi.clear();
    jet4Mass.clear();

    nJets8 = 0;
    jet8Pt.clear();
    jet8Eta.clear();
    jet8Phi.clear();
    jet8Mass.clear();

    nGenJets4 = 0;
    genJet4Pt.clear();
    genJet4Eta.clear();
    genJet4Phi.clear();
    genJet4Mass.clear();

    nGenJets8 = 0;
    genJet8Pt.clear();
    genJet8Eta.clear();
    genJet8Phi.clear();
    genJet8Mass.clear();

    nPuppiJets4 = 0;
    puppiJet4Pt.clear();
    puppiJet4Eta.clear();
    puppiJet4Phi.clear();
    puppiJet4Mass.clear();

    nPuppiJets8 = 0;
    puppiJet8Pt.clear();
    puppiJet8Eta.clear();
    puppiJet8Phi.clear();
    puppiJet8Mass.clear();

    nMuons = 0;
    muPt.clear();
    muEta.clear();
    muPhi.clear();

    nElectrons = 0;
    elePt.clear();
    eleEta.clear();
    elePhi.clear();

    met = -999;
    phiMet = -999;
    genMet = -999;

    MR4 = -999;
    RSquared4 = -999;
    MRGen4 = -999;
    RSquaredGen4 = -999;
    MRPuppi = -999;
    RSquaredPuppi = -999;
    MR8 = -999;
    RSquared8 = -999;
    MRGen8 = -999;
    RSquaredGen8 = -999;
    MRPuppi8 = -999;
    RSquaredPuppi8 = -999;

    nPV = -1;
}

//------ Method called for each event ------//

void RazorJetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    //initialize
    loadEvent(iEvent); 

    //select the primary vertex, if any
    if (vertices->empty()) return; // skip the event if no PV found
    //const reco::Vertex &PV = vertices->front();
    nPV = vertices->size();

    //muons
    for(const pat::Muon &mu : *muons){
        if(mu.pt() < 5 || !mu.isLooseMuon()) continue;
        nMuons++;
        muPt.push_back(mu.pt());
        muEta.push_back(mu.eta());
        muPhi.push_back(mu.phi());
    }

    //electrons
    for(const pat::Electron &ele : *electrons){
        if(ele.pt() < 5) continue;
        nElectrons++;
        elePt.push_back(ele.pt());
        eleEta.push_back(ele.eta());
        elePhi.push_back(ele.phi());
    }   

    //AK4 jets
    vector<TLorentzVector> goodJets;
    for (const pat::Jet &j : *jets) {
        if(j.pt() < 40) continue;
        if(fabs(j.eta()) > 2.4) continue;  
        nJets4++;
        TLorentzVector newJet(j.px(), j.py(), j.pz(), j.energy());
        goodJets.push_back(newJet);
        jet4Pt.push_back(j.pt());
        jet4Eta.push_back(j.eta());
        jet4Phi.push_back(j.phi());
        jet4Mass.push_back(j.mass());
    }
    
    //AK4 PUPPI jets
    vector<TLorentzVector> goodPuppiJets;
    for(const reco::PFJet &pJ : *puppiJets){
        if(pJ.pt() < 40) continue;
        if(fabs(pJ.eta()) > 2.4) continue;
        nPuppiJets4++;
        TLorentzVector newJet(pJ.px(), pJ.py(), pJ.pz(), pJ.energy());
        goodPuppiJets.push_back(newJet);
        puppiJet4Pt.push_back(pJ.pt());
        puppiJet4Eta.push_back(pJ.eta());
        puppiJet4Phi.push_back(pJ.phi());
        puppiJet4Mass.push_back(pJ.mass());
    }

    //AK4 genJets
    vector<TLorentzVector> goodGenJets;
    for(const reco::GenJet &j : *genjets){
        if(j.pt() < 40) continue;
        if(fabs(j.eta()) > 2.4) continue;
        nGenJets4++;
        TLorentzVector newJet(j.px(), j.py(), j.pz(), j.energy());
        goodGenJets.push_back(newJet);
        genJet4Pt.push_back(j.pt());
        genJet4Eta.push_back(j.eta());
        genJet4Phi.push_back(j.phi());
        genJet4Mass.push_back(j.mass());
    }

    //AK8 jets
    vector<TLorentzVector> goodJets8;
    for (const pat::Jet &j : *fatjets) {
        if(j.pt() < 40) continue;
        if(fabs(j.eta()) > 2.4) continue;  
        nJets8++;
        TLorentzVector newJet(j.px(), j.py(), j.pz(), j.energy());
        goodJets8.push_back(newJet);
        jet8Pt.push_back(j.pt());
        jet8Eta.push_back(j.eta());
        jet8Phi.push_back(j.phi());
        jet8Mass.push_back(j.mass());
    }
    
    //AK8 PUPPI jets
    vector<TLorentzVector> goodPuppiJets8;
    for(const reco::PFJet &pJ : *puppiJets8){
        if(pJ.pt() < 40) continue;
        if(fabs(pJ.eta()) > 2.4) continue;
        nPuppiJets8++;
        TLorentzVector newJet(pJ.px(), pJ.py(), pJ.pz(), pJ.energy());
        goodPuppiJets8.push_back(newJet);
        puppiJet8Pt.push_back(pJ.pt());
        puppiJet8Eta.push_back(pJ.eta());
        puppiJet8Phi.push_back(pJ.phi());
        puppiJet8Mass.push_back(pJ.mass());
    }

    //AK8 genJets
    vector<TLorentzVector> goodGenJets8;
    for(const reco::GenJet &j : *genJets8){
        if(j.pt() < 40) continue;
        if(fabs(j.eta()) > 2.4) continue;
        TLorentzVector newJet(j.px(), j.py(), j.pz(), j.energy());
        goodGenJets8.push_back(newJet);
        nGenJets8++;
        genJet8Pt.push_back(j.pt());
        genJet8Eta.push_back(j.eta());
        genJet8Phi.push_back(j.phi());
        genJet8Mass.push_back(j.mass());
    }
   
    //process MET
    const pat::MET &Met = mets->front();
    //store MET variables in tree
    met = Met.pt();
    phiMet = Met.phi();
    genMet = Met.genMET()->pt();
   
    //compute the razor variables using the selected jets
    if(goodJets.size() > 1){
        vector<TLorentzVector> hemispheres = getHemispheres(goodJets);
        TLorentzVector pfMet(Met.px(), Met.py(), 0.0, 0.0);
        MR4 = computeMR(hemispheres[0], hemispheres[1]);
        RSquared4 = computeR2(hemispheres[0], hemispheres[1], pfMet);
    }

    //compute the razor variables using the selected genJets
    if(goodGenJets.size() > 1){
        vector<TLorentzVector> hemispheres = getHemispheres(goodGenJets);
        TLorentzVector pfMet(Met.genMET()->px(), Met.genMET()->py(), 0.0, 0.0);
        MRGen4 = computeMR(hemispheres[0], hemispheres[1]);
        RSquaredGen4 = computeR2(hemispheres[0], hemispheres[1], pfMet);
    }

    //compute the razor variables using the selected AK4 Puppi Jets
    if(goodPuppiJets.size() > 1){
        vector<TLorentzVector> hemispheres = getHemispheres(goodPuppiJets);
        //TODO: recompute MET using the Puppi particles (might not be very good)
        TLorentzVector pfMet(Met.px(), Met.py(), 0.0, 0.0);
        MRPuppi = computeMR(hemispheres[0], hemispheres[1]);
        RSquaredPuppi = computeR2(hemispheres[0], hemispheres[1], pfMet);
    }

    //compute the razor variables using the selected fat jets
    if(goodJets8.size() > 1){
        vector<TLorentzVector> hemispheres = getHemispheres(goodJets8);
        TLorentzVector pfMet(Met.px(), Met.py(), 0.0, 0.0);
        MR8 = computeMR(hemispheres[0], hemispheres[1]);
        RSquared8 = computeR2(hemispheres[0], hemispheres[1], pfMet);
    }

    //compute the razor variables using the selected fat genJets
    if(goodGenJets8.size() > 1){
        vector<TLorentzVector> hemispheres = getHemispheres(goodGenJets8);
        TLorentzVector pfMet(Met.genMET()->px(), Met.genMET()->py(), 0.0, 0.0);
        MRGen8 = computeMR(hemispheres[0], hemispheres[1]);
        RSquaredGen8 = computeR2(hemispheres[0], hemispheres[1], pfMet);
    }

    //compute the razor variables using the selected AK8 Puppi Jets
    if(goodPuppiJets8.size() > 1){
        vector<TLorentzVector> hemispheres = getHemispheres(goodPuppiJets8);
        //TODO: recompute MET using the Puppi particles (might not be very good)
        TLorentzVector pfMet(Met.px(), Met.py(), 0.0, 0.0);
        MRPuppi8 = computeMR(hemispheres[0], hemispheres[1]);
        RSquaredPuppi8 = computeR2(hemispheres[0], hemispheres[1], pfMet);
    }

    //fill the tree
    outputTree->Fill();
}

//define this as a plug-in
DEFINE_FWK_MODULE(RazorJetAnalyzer);
