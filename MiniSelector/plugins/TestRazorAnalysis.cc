// Class:      TestRazorAnalysis
// Description: template for a class extending MiniSelector

// Author:  Dustin James Anderson

//------ Class declaration ------//

#include "MiniSelector.h"

//analysis-specific includes

using namespace std;

class TestRazorAnalysis : public MiniSelector {
    public:
        //analyzer constructor and destructor
        explicit TestRazorAnalysis(const edm::ParameterSet&);
        ~TestRazorAnalysis();

        //helper functions
        void loadEvent(const edm::Event& iEvent); //call at the start of each event

    private:
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

        //----- Member data ------//

        //non-default inputs 

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
        int nJets;
        vector<double> jetPt, jetEta, jetPhi, jetMass;

        //genjet variables
        int nGenJets;
        vector<double> genJetPt, genJetEta, genJetPhi, genJetMass;

        //MET variables
        double met;
        double phiMet;
        double genMet;

        //razor variables
        double MR, RSq;

        //vertex
        int nPV;
};

// constructors and destructor
TestRazorAnalysis::TestRazorAnalysis(const edm::ParameterSet& iConfig): MiniSelector(iConfig)
{
    //declare the TFileService for output
    edm::Service<TFileService> fs;
    //set up output tree
    outputTree = fs->make<TTree>("outputTree", "selected miniAOD information");

    //initialize tree branches
    outputTree->Branch("nJets", &nJets, "nJets/I");
    outputTree->Branch("jetPt", "std::vector<double>", &jetPt);
    outputTree->Branch("jetEta", "std::vector<double>", &jetEta);
    outputTree->Branch("jetPhi", "std::vector<double>", &jetPhi);
    outputTree->Branch("jetMass", "std::vector<double>", &jetMass);

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

    outputTree->Branch("MR", &MR, "MR/D");
    outputTree->Branch("RSq", &RSq, "RSq/D");

    outputTree->Branch("nPV", &nPV, "nPV/I");
}

TestRazorAnalysis::~TestRazorAnalysis()
{
}

// member functions

void TestRazorAnalysis::loadEvent(const edm::Event& iEvent){

    //load default inputs
    MiniSelector::loadEvent(iEvent);
    
    //read extra inputs
    
    //reset tree variables
    nJets = 0;
    jetPt.clear();
    jetEta.clear();
    jetPhi.clear();
    jetMass.clear();

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

    MR = -999;
    RSq = -999;

    nPV = -1;
}

//------ Method called for each event ------//

void TestRazorAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    //initialize
    loadEvent(iEvent); 

    //select the primary vertex, if any
    if (vertices->empty()) return; // skip the event if no PV found
    //const reco::Vertex &PV = vertices->front();
    nPV = vertices->size();

    //muons
    vector<TLorentzVector> goodMuons;
    for(const pat::Muon &mu : *muons){
        if(mu.pt() < 5 || !mu.isLooseMuon()) continue;
        nMuons++;
        muPt.push_back(mu.pt());
        muEta.push_back(mu.eta());
        muPhi.push_back(mu.phi());
        TLorentzVector newMuon(mu.px(), mu.py(), mu.pz(), mu.energy());
        goodMuons.push_back(newMuon);
    }

    //electrons
    vector<TLorentzVector> goodElectrons;
    for(const pat::Electron &ele : *electrons){
        if(ele.pt() < 5) continue;
        nElectrons++;
        elePt.push_back(ele.pt());
        eleEta.push_back(ele.eta());
        elePhi.push_back(ele.phi());
        TLorentzVector newElectron(ele.px(), ele.py(), ele.pz(), ele.energy());
        goodElectrons.push_back(newElectron);
    }   

    //AK4 jets
    vector<TLorentzVector> goodJets;
    for (const pat::Jet &j : *jets) {
        if(j.pt() < 20) continue;
        if(fabs(j.eta()) > 5.0) continue;  
        TLorentzVector newJet(j.px(), j.py(), j.pz(), j.energy());
        nJets++;
        goodJets.push_back(newJet);
        jetPt.push_back(j.pt());
        jetEta.push_back(j.eta());
        jetPhi.push_back(j.phi());
        jetMass.push_back(j.mass());
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
        MR = computeMR(hemispheres[0], hemispheres[1]);
        RSq = computeR2(hemispheres[0], hemispheres[1], pfMet);
    }

    //fill the tree
    outputTree->Fill();
}

//define this as a plug-in
DEFINE_FWK_MODULE(TestRazorAnalysis);
