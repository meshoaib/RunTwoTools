// Class:      ExampleAnalyzer
// Description: template for a class extending MiniSelector

// Author:  Dustin James Anderson

//------ Class declaration ------//

#include "MiniSelector.h"

//analysis-specific includes
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

class ExampleAnalyzer : public MiniSelector {
    public:
        //analyzer constructor and destructor
        explicit ExampleAnalyzer(const edm::ParameterSet&);
        ~ExampleAnalyzer();

        //helper functions
        void loadEvent(const edm::Event& iEvent); //call at the start of each event
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

        //virtual void beginJob() override;
        //virtual void endJob() override;
        //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
        //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
        //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
        //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;


        //----- Member data ------//

        //non-default inputs (ex: uncorrected PF MET)
        edm::EDGetTokenT<reco::PFMETCollection> rawPfMetToken_;
        edm::EDGetTokenT<reco::PFJetCollection> ak5PfToken_;
        
        edm::Handle<reco::PFMETCollection> rawPfMet;
        edm::Handle<reco::PFJetCollection> ak5Pf;

        //output tree
        TTree *outputTree;

        //------ Variables for tree ------//

        //MET variables
        double rawMet;
        int nMuon;
        int nAk5Jets;
};

// constants, enums and typedefs

// static data member definitions

// constructors and destructor
ExampleAnalyzer::ExampleAnalyzer(const edm::ParameterSet& iConfig): MiniSelector(iConfig),
    rawPfMetToken_(consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("rawPfMet"))),
    ak5PfToken_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("ak5PFJetsCHS")))
{
    //declare the TFileService for output
    edm::Service<TFileService> fs;
    //set up output tree
    outputTree = fs->make<TTree>("outputTree", "selected miniAOD information");

    //initialize tree branches
    outputTree->Branch("rawMet", &rawMet, "rawMet/D");
    outputTree->Branch("nMuon", &nMuon, "nMuon/I");
    outputTree->Branch("nAk5Jets", &nAk5Jets, "nAk5Jets/I");
}


ExampleAnalyzer::~ExampleAnalyzer()
{
}


// member functions

void ExampleAnalyzer::loadEvent(const edm::Event& iEvent){

    //load default inputs
    MiniSelector::loadEvent(iEvent);
    
    //read extra inputs
    iEvent.getByToken(rawPfMetToken_, rawPfMet);
    iEvent.getByToken(ak5PfToken_, ak5Pf);
    
    //reset tree variables
    rawMet = -999;
    nMuon = 0;
    nAk5Jets = 0;
}

//------ Method called for each event ------//

void ExampleAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    //initialize
    loadEvent(iEvent); 

    //print trigger information
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    std::cout << "\n === TRIGGER PATHS === " << std::endl;
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
        std::cout << "Trigger " << names.triggerName(i) <<
                ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
                ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)")
                << std::endl;
    }

    //select the primary vertex, if any
    if (vertices->empty()) return; // skip the event if no PV found
    const reco::Vertex &PV = vertices->front();

    for(const pat::Muon &mu : *muons){
        if(mu.pt() < 5 || !mu.isTightMuon(PV)) continue;
        nMuon++;
    }

    for(const reco::Jet &jet : *ak5Pf){
        if(jet.pt() < 40 || fabs(jet.eta()) > 2.4) continue;
        nAk5Jets++;
    }

    //store raw MET
    rawMet = rawPfMet->front().pt();

    //fill the tree
    outputTree->Fill();
}


//------ Method called once each job just before starting event loop ------//
/*
void ExampleAnalyzer::beginJob(){
}
*/

//------ Method called once each job just after ending the event loop ------//
/*
void ExampleAnalyzer::endJob(){
}
*/

//------ Method called when starting to processes a run ------//
/*
void ExampleAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&){
}
*/

//------ Method called when ending the processing of a run ------//
/*
void ExampleAnalyzer::endRun(edm::Run const&, edm::EventSetup const&){
}
*/

//------ Method called when starting to processes a luminosity block ------//
/*
void ExampleAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){
}
*/

//------ Method called when ending the processing of a luminosity block ------//
/*
void ExampleAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){
}
*/

//------ Method fills 'descriptions' with the allowed parameters for the module ------//
void ExampleAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ExampleAnalyzer);
