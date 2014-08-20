#include "NoTrees.hh"
#include "fastjet/internal/base.hh"
#include "fastjet/PseudoJet.hh"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

using namespace std;
using namespace fastjet;

//FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

enum ParticleType {
    X=0,     // undefined
    h,       // charged hadron
    e,       // electron 
    mu,      // muon 
    photon,   // photon
    h0,      // neutral hadron
    h_HF,        // HF tower identified as a hadron
    egamma_HF    // HF tower identified as an EM particle
};

//......................
class puppiCleanContainer{
    public:
        // default ctor
        puppiCleanContainer(std::vector<pat::PackedCandidate> inParticles,double iTracker=2.5, bool iExperiment=false,bool iTuned=true);
        ~puppiCleanContainer(); 
        std::vector<fastjet::PseudoJet> genParticles(){ return _genParticles; }
        std::vector<fastjet::PseudoJet> pfParticles(){ return _pfParticles; }    
        std::vector<fastjet::PseudoJet> pvParticles(){ return _chargedPV; }        
        std::vector<fastjet::PseudoJet> puParticles(){ return _chargedNoPV; }    
        std::vector<fastjet::PseudoJet> pfchsParticles(){ return _pfchsParticles; }    
        std::vector<fastjet::PseudoJet> puppiEvent     (int iOpt,double iQuant);

        std::vector<float> getPuppiWeights_chLV(){ return puppiWeights_chLV; };
        std::vector<float> getPuppiWeights_all(){ return puppiWeights_all; };
        std::vector<float> getPuppiAlphas_chLV(){ return alphas_chLV; };
        std::vector<float> getPuppiAlphas_all(){ return alphas_all; };

    protected:

        //Helper Functions
        double  goodVar  (fastjet::PseudoJet &iPart,std::vector<fastjet::PseudoJet> &iParts, int iOpt);    
        void    getRMSAvg(int iOpt,std::vector<fastjet::PseudoJet> &iConstits,std::vector<fastjet::PseudoJet> &iParticles,double iQuant,double iPtRMS);
        double  compute  (int iOpt,double iVal,double iMed,double iRMS,double iChi2Exp);
        double  getChi2FromdZ(double iDZ);
        std::vector<pat::PackedCandidate> _recoParticles;
        std::vector<PseudoJet> _pfParticles;
        std::vector<PseudoJet> _pfchsParticles;    
        std::vector<PseudoJet> _genParticles;
        std::vector<PseudoJet> _chargedPV;
        std::vector<PseudoJet> _chargedNoPV;
        std::vector<double> _vals;
        double fMed;
        double fRMS;
        double fMedHEta;
        double fRMSHEta;
        double fNeutralMinE;
        bool _isTuned;
        bool _isExperiment;
        double fTrackerEta; 

        std::vector<float> puppiWeights_chLV;
        std::vector<float> puppiWeights_all;
        std::vector<float> alphas_chLV;
        std::vector<float> alphas_all;

};

//FASTJET_END_NAMESPACE

