// Class:      L1TkMuonHPSPFTauAnalyzer
//
// Original Author:  Sandeep Bhowmik
//         Created:  Tue, 12 Mar 2019 18:38:39 GMT
//
#include "FWCore/Framework/interface/one/EDAnalyzer.h" 
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"     
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include <FWCore/Framework/interface/Frameworkfwd.h>
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TTree.h>

#include "DataFormats/Phase2L1Taus/interface/L1HPSPFTau.h"
#include "DataFormats/Phase2L1Taus/interface/L1HPSPFTauFwd.h"
#include "DataFormats/Phase2L1Taus/interface/L1TkMuonHPSPFTau.h"
#include "DataFormats/Phase2L1Taus/interface/L1TkMuonHPSPFTauFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"  
#include "DataFormats/L1Trigger/interface/Tau.h" 
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include <DataFormats/PatCandidates/interface/Tau.h> 
#include "L1TauAnalyzerPhase2/L1TauAnalyzerPhase2/plugins/GenVertexProducer.h"
//#include "TMVA/Reader.h"
#include "L1TauAnalyzerPhase2/L1TauAnalyzerPhase2/interface/TMVAInterface.h" // TMVAInterface
#include "L1TauAnalyzerPhase2/L1TauAnalyzerPhase2/interface/XGBInterface.h" // XGBInterface

class L1TkMuonHPSPFTauAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit L1TkMuonHPSPFTauAnalyzer(const edm::ParameterSet&);
  ~L1TkMuonHPSPFTauAnalyzer();
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  void Initialize();
  
  TTree *tree_;
  std::string treeName_;
  
  ULong64_t       indexevents_;
  Int_t           runNumber_;
  Int_t           lumi_;
  float MC_weight_;
  std::vector<float> genMuonPt_;
  std::vector<float> genMuonEta_;
  std::vector<float> genMuonPhi_;
  std::vector<int> genMuonCharge_;
  std::vector<float> genTauPt_;
  std::vector<float> genTauEta_;
  std::vector<float> genTauPhi_;
  std::vector<int> genTauCharge_;
  std::vector<Bool_t> isGenMatched_;
  std::vector<int> recoTauDecayMode_;
  std::vector<float> recoTauPt_;
  std::vector<float> recoTauEta_;
  std::vector<float> recoTauPhi_;
  std::vector<int> recoTauCharge_;
  std::vector<Bool_t> isRecoMatched_;
  std::vector<int> recoGMTauDecayMode_;
  std::vector<float> recoGMTauPt_;
  std::vector<float> recoGMTauEta_;
  std::vector<float> recoGMTauPhi_;
  std::vector<int> recoGMTauCharge_;
  std::vector<Bool_t> isRecoGMMatched_;
  std::vector<float> l1TkMuonPt_;
  std::vector<float> l1TkMuonEta_;
  std::vector<float> l1TkMuonPhi_;
  std::vector<int> l1TkMuonCharge_;
  std::vector<int> l1PFTauType_;
  std::vector<float> l1PFTauPt_;
  std::vector<float> l1PFTauEta_;
  std::vector<float> l1PFTauPhi_;
  std::vector<int> l1PFTauCharge_;
  std::vector<float> l1PFTauIso_;
  std::vector<float> l1PFTauIsoOverPt_;
  std::vector<Bool_t> l1PFTauTightIso_;
  std::vector<Bool_t> l1PFTauMediumIso_;
  std::vector<Bool_t> l1PFTauLooseIso_;
  std::vector<Bool_t> l1PFTauVLooseIso_;
  std::vector<Bool_t> l1PFTauTightRelIso_;
  std::vector<Bool_t> l1PFTauMediumRelIso_;
  std::vector<Bool_t> l1PFTauLooseRelIso_;
  std::vector<Bool_t> l1PFTauVLooseRelIso_;
  std::vector<float> l1PFTauZ_;
  std::vector<float> single_l1TkMuonPt_;
  std::vector<float> single_l1TkMuonEta_;
  std::vector<float> single_l1TkMuonPhi_;
  std::vector<int> single_l1TkMuonCharge_;
  std::vector<int> single_l1PFTauType_;
  std::vector<float> single_l1PFTauPt_;
  std::vector<float> single_l1PFTauEta_;
  std::vector<float> single_l1PFTauPhi_;
  std::vector<int> single_l1PFTauCharge_;
  std::vector<float> single_l1PFTauIso_;
  std::vector<Bool_t> single_l1PFTauTightIso_;
  std::vector<Bool_t> single_l1PFTauMediumIso_;
  std::vector<Bool_t> single_l1PFTauLooseIso_;
  std::vector<Bool_t> single_l1PFTauVLooseIso_;
  std::vector<Bool_t> single_l1PFTauTightRelIso_;
  std::vector<Bool_t> single_l1PFTauMediumRelIso_;
  std::vector<Bool_t> single_l1PFTauLooseRelIso_;
  std::vector<Bool_t> single_l1PFTauVLooseRelIso_;
  std::vector<float> single_l1PFTauZ_;

  std::vector<float> l1PFTauLeadTrackPtOverTauPt_;
  std::vector<float> l1PFTauChargedIso_;
  std::vector<float> l1PFTauNeutralIso_;
  std::vector<float> l1PFTauChargedIsoPileup_;
  std::vector<float> l1PFTauRho_;
  std::vector<float> l1PFTauNSignalChargedHadrons_;
  std::vector<float> l1PFTauNSignalElectrons_;
  std::vector<float> l1PFTauNSignalPhotons_;
  std::vector<float> l1PFTauNSignalChargedPFCands_;
  std::vector<float> l1PFTauSignalChargeSum_;
  std::vector<float> l1PFTauStripPtOverTauPt_;
  std::vector<float> l1PFTauStripMassOverTauPt_;
  std::vector<float> l1PFTauStripMassOverStripPt_;

  std::vector<float> l1PFTauStripPt_;
  std::vector<float> l1PFTauLeadTrackPt_;
  std::vector<int>   l1PFTauVtxIndex_;
  std::vector<float> l1PFTaudz_;
  std::vector<float> l1PFTauSumTrackPtOfVtx_;
  std::vector<float> l1PFTauLeadTrackHoverE_;
  std::vector<float> l1PFTauHoverE_;
  std::vector<float> l1PFTauSignalTrackMass_;
  std::vector<float> l1PFTauNStripElectrons_;
  std::vector<float> l1PFTauNStripPhotons_;
  std::vector<float> l1PFTauDeltaRLeadTrackStrip_;
  double genVertex_;
  int l1VertexN_;
  int recoVertexN_;

  std::vector<float> l1PFTauBDT_;

  bool fillBDT_;
  TTree *treeBDT_;
  std::string treeBDTName_;
  std::string bdtRootFileName_;
  TFile* bdtRootFile_;

  bool applyBDT_;
  std::string bdtInputFileName_;
  /*
  TMVA::Reader *tmva_reader_;
  float bdt_l1PFTauPt_;
  float bdt_l1PFTauEta_;
  float bdt_l1PFTauIso_;
  float bdt_l1PFTauNeutralIso_;
  float bdt_l1PFTauChargedIsoPileup_;
  float bdt_l1PFTauNSignalPhotons_;
  float bdt_l1PFTauSignalChargeSum_;
  float bdt_l1PFTauStripPt_;
  float bdt_l1PFTauLeadTrackPt_;
  float bdt_l1PFTaudz_;
  float bdt_l1PFTauLeadTrackHoverE_;
  float bdt_l1PFTauHoverE_;
  float bdt_l1PFTauSignalTrackMass_;
  */
  TMVAInterface *bdt_tmva_L1HPSPFTau;
  XGBInterface *bdt_xgb_L1HPSPFTau;

  bool createHistRoorFile_;
  std::string histRootFileName_;
  TFile* histRootFile_;
  TH1F* hist_genMuonPt_;
  TH1F* hist_genMuonEta_;
  TH1F* hist_genMuonPhi_;
  TH1F* hist_genTauPt_;
  TH1F* hist_genTauEta_;
  TH1F* hist_genTauPhi_;
  TH1F* hist_isGenMatched_;
  TH1F* hist_recoTauPt_;
  TH1F* hist_recoTauEta_;
  TH1F* hist_recoTauPhi_;
  TH1F* hist_isRecoMatched_;
  TH1F* hist_recoGMTauPt_;
  TH1F* hist_recoGMTauEta_;
  TH1F* hist_recoGMTauPhi_;
  TH1F* hist_isRecoGMMatched_;
  TH1F* hist_l1TkMuonPt_;
  TH1F* hist_l1TkMuonEta_;
  TH1F* hist_l1TkMuonPhi_;
  TH1F* hist_l1PFTauPt_;
  TH1F* hist_l1PFTauEta_;
  TH1F* hist_l1PFTauPhi_;
  TH1F* hist_l1PFTauReso_vs_Gen_;
  TH1F* hist_l1PFTauReso_vs_Reco_;
  TH1F* hist_l1PFTauReso_vs_RecoGM_;

  // ----------member data ---------------------------

  bool debug_;
  bool isReco_;
  double min_muon_pt_;
  double min_tau_pt_;
  double max_eta_;
  typedef std::vector<int> vint;
  vint muonPDGIds_;
  edm::EDGetTokenT<GenEventInfoProduct>            genTagToken_;
  edm::EDGetTokenT<double>                         genVertexToken_;
  edm::EDGetTokenT<std::vector<l1t::Vertex>>       l1VertexToken_;
  edm::EDGetTokenT<reco::GenParticleCollection>    genMuonToken_;
  edm::EDGetTokenT<std::vector<reco::GenJet>>      genTauToken_;
  edm::EDGetTokenT<l1t::TauBxCollection>           l1TauToken_;
  edm::EDGetTokenT<l1t::L1HPSPFTauCollection>      l1PFTauToken_;
  edm::EDGetTokenT<l1t::L1TkMuonHPSPFTauCollection>l1TkMuonPFTauToken_;
  edm::EDGetTokenT<l1t::L1TkMuonParticleCollection>l1TkMuonToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>>      recoVertexToken_;
  edm::EDGetTokenT<std::vector<pat::Tau>>          recoTauToken_;
  edm::EDGetTokenT<pat::TauRefVector>              recoGMTauToken_;

};


L1TkMuonHPSPFTauAnalyzer::L1TkMuonHPSPFTauAnalyzer(const edm::ParameterSet& iConfig)
  : debug_          (iConfig.getUntrackedParameter<bool>("debug", false))
  , isReco_         (iConfig.getUntrackedParameter<bool>("isReco", false))
  , min_muon_pt_         (iConfig.getUntrackedParameter<double>("min_muon_pt", 10))
  , min_tau_pt_         (iConfig.getUntrackedParameter<double>("min_tau_pt", 20))
  , max_eta_        (iConfig.getUntrackedParameter<double>("max_eta", 2.4))
  , muonPDGIds_     (iConfig.getUntrackedParameter<vint>("muonPDGIds"))
  , genTagToken_    (consumes<GenEventInfoProduct>          (iConfig.getParameter<edm::InputTag>("genTagToken")))
  , genVertexToken_ (consumes<double>                       (iConfig.getParameter<edm::InputTag>("genVertexToken")))
  , l1VertexToken_  (consumes<std::vector<l1t::Vertex>>     (iConfig.getParameter<edm::InputTag>("l1VertexToken")))
  , genMuonToken_    (consumes<reco::GenParticleCollection> (iConfig.getParameter<edm::InputTag>("genMuonToken")))
  , genTauToken_    (consumes<std::vector<reco::GenJet>>    (iConfig.getParameter<edm::InputTag>("genTauToken")))
  , l1TauToken_     (consumes<l1t::TauBxCollection>         (iConfig.getParameter<edm::InputTag>("l1TauToken")))
  , l1PFTauToken_   (consumes<l1t::L1HPSPFTauCollection>    (iConfig.getParameter<edm::InputTag>("l1PFTauToken")))
  , l1TkMuonPFTauToken_(consumes<l1t::L1TkMuonHPSPFTauCollection>(iConfig.getParameter<edm::InputTag>("l1TkMuonPFTauToken")))
  , l1TkMuonToken_(consumes<l1t::L1TkMuonParticleCollection>(iConfig.getParameter<edm::InputTag>("l1TkMuonToken")))
  , recoVertexToken_(consumes<std::vector<reco::Vertex>>    (iConfig.getParameter<edm::InputTag>("recoVertexToken")))
  , recoTauToken_   (consumes<std::vector<pat::Tau>>        (iConfig.getParameter<edm::InputTag>("recoTauToken")))
  , recoGMTauToken_ (consumes<pat::TauRefVector>            (iConfig.getParameter<edm::InputTag>("recoGMTauToken")))

{
   //now do what ever initialization is needed
  treeName_             = iConfig.getParameter<std::string>("treeName");
  edm::Service<TFileService> fs;
  tree_                 = fs -> make<TTree>(treeName_.c_str(), treeName_.c_str());
  fillBDT_              = iConfig.getUntrackedParameter<bool>("fillBDT", false);
  bdtRootFileName_      = iConfig.getParameter<std::string>("bdtRootFileName");
  bdtRootFile_          = new TFile(bdtRootFileName_.c_str(), "RECREATE");
  treeBDTName_          = iConfig.getParameter<std::string>("treeBDTName");
  treeBDT_              = new TTree(treeBDTName_.c_str(), treeBDTName_.c_str());
  createHistRoorFile_   = iConfig.getUntrackedParameter<bool>("createHistRoorFile", false);
  histRootFileName_     = iConfig.getParameter<std::string>("histRootFileName");
  histRootFile_         = new TFile(histRootFileName_.c_str(), "RECREATE");
  applyBDT_             = iConfig.getUntrackedParameter<bool>("applyBDT", false);
  bdtInputFileName_     = iConfig.getParameter<std::string>("bdtInputFileName");
  return;
}


L1TkMuonHPSPFTauAnalyzer::~L1TkMuonHPSPFTauAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
L1TkMuonHPSPFTauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  Initialize();

  if(debug_){
    std::cout<<" Starting L1HPSPFTau Analyzer ............     "<< std::endl;  
  }
    using namespace edm;

   indexevents_ = iEvent.id().event();
   runNumber_ = iEvent.id().run();
   lumi_ = iEvent.luminosityBlock();

   edm::Handle<GenEventInfoProduct>        genEvt;
   try {iEvent.getByToken(genTagToken_,    genEvt);}
   catch (...) {;}
   if(genEvt.isValid()) MC_weight_ = genEvt->weight();

   edm::Handle<double>                     genVertexHandle;
   iEvent.getByToken(genVertexToken_,      genVertexHandle);
   genVertex_ = *genVertexHandle;

   edm::Handle<std::vector<l1t::Vertex> >  l1VertexHandle;
   iEvent.getByToken(l1VertexToken_,       l1VertexHandle);
   l1VertexN_ = l1VertexHandle->size();

   edm::Handle<reco::GenParticleCollection> genMuonHandle;
   iEvent.getByToken(genMuonToken_,         genMuonHandle);

   edm::Handle<std::vector<reco::GenJet>>  genTauHandle;
   iEvent.getByToken(genTauToken_,         genTauHandle);

   edm::Handle< BXVector<l1t::Tau> >       l1TauHandle;
   iEvent.getByToken(l1TauToken_,          l1TauHandle);

   edm::Handle<l1t::L1HPSPFTauCollection>  l1PFTauHandle;
   iEvent.getByToken(l1PFTauToken_,        l1PFTauHandle);

   edm::Handle<l1t::L1TkMuonHPSPFTauCollection>  l1TkMuonPFTauHandle;
   iEvent.getByToken(l1TkMuonPFTauToken_,        l1TkMuonPFTauHandle);

   edm::Handle<l1t::L1TkMuonParticleCollection>  l1TkMuonHandle;
   iEvent.getByToken(l1TkMuonToken_,        l1TkMuonHandle);

   for(auto genMuon : *genMuonHandle){
     if ( muonPDGIds_.size() > 0 ){
       bool isSelected = false;
       for ( auto PDGId : muonPDGIds_ ){
	 if ( genMuon.pdgId() == PDGId ) isSelected = true;
       }
       if ( !isSelected ) continue;
     }
     if (fabs(genMuon.eta())>max_eta_)
       continue;
     if (fabs(genMuon.pt())<min_muon_pt_)
       continue;

     genMuonPt_.push_back(genMuon.pt());
     genMuonEta_.push_back(genMuon.eta());
     genMuonPhi_.push_back(genMuon.phi());
     genMuonCharge_.push_back(genMuon.charge());

     hist_genMuonPt_->Fill(genMuon.pt());
     hist_genMuonEta_->Fill(genMuon.eta());
     hist_genMuonPhi_->Fill(genMuon.phi());

     if(debug_){
       std::cout<<" GenMuon pt "<<genMuon.pt()<<" eta "<< genMuon.eta()<<" phi "<< genMuon.phi()<<" charge "<< genMuon.charge()<<std::endl;
     }
   }

   for(auto genTau : *genTauHandle){
     if (fabs(genTau.eta())>max_eta_)
       continue;
     //if (fabs(genTau.pt())<min_muon_pt_)
     //continue;
     genTauPt_.push_back(genTau.pt());
     genTauEta_.push_back(genTau.eta());
     genTauPhi_.push_back(genTau.phi());
     genTauCharge_.push_back(genTau.charge());
     bool isMatched = false;
     for(auto l1PFTau : *l1PFTauHandle){
       double deltaEta = l1PFTau.eta() - genTau.eta();
       double deltaPhi = l1PFTau.phi() - genTau.phi();
       if ( (deltaEta*deltaEta + deltaPhi*deltaPhi) < 0.25 ){
         isMatched = true;
         hist_l1PFTauReso_vs_Gen_->Fill(l1PFTau.pt() / genTau.pt());
         break;
       }
     }
     isGenMatched_.push_back(isMatched);

     hist_isGenMatched_->Fill(isMatched);
     hist_genTauPt_->Fill(genTau.pt());
     hist_genTauEta_->Fill(genTau.eta());
     hist_genTauPhi_->Fill(genTau.phi());

     if(debug_){
       std::cout<<" GenTau pt "<<genTau.pt()<<" eta "<< genTau.eta()<<" phi "<< genTau.phi()<<" charge "<< genTau.charge()<<std::endl;
       std::cout<<" GenTau Z " <<genTau.vertex().z() << std::endl;
     }
   }

   for(auto l1PFTau : *l1PFTauHandle){
     if(fabs(l1PFTau.pt())<min_tau_pt_)
       continue;
     if(fabs(l1PFTau.eta())>max_eta_)
       continue;

     single_l1PFTauPt_.push_back(l1PFTau.pt());
     single_l1PFTauEta_.push_back(l1PFTau.eta());
     single_l1PFTauPhi_.push_back(l1PFTau.phi());
     single_l1PFTauCharge_.push_back(l1PFTau.charge());
     single_l1PFTauType_.push_back(l1PFTau.tauType());
     single_l1PFTauIso_.push_back(l1PFTau.sumChargedIso());
     single_l1PFTauTightIso_.push_back(l1PFTau.passTightIso());
     single_l1PFTauMediumIso_.push_back(l1PFTau.passMediumIso());
     single_l1PFTauLooseIso_.push_back(l1PFTau.passLooseIso());
     single_l1PFTauVLooseIso_.push_back(l1PFTau.passVLooseIso());
     double l1PFTauZ = 1000;
     if ( l1PFTau.leadChargedPFCand().isNonnull() && l1PFTau.leadChargedPFCand()->pfTrack().isNonnull())
       {
         l1PFTauZ = l1PFTau.leadChargedPFCand()->pfTrack()->vertex().z();
       }
     single_l1PFTauZ_.push_back(l1PFTauZ);

     if(l1PFTau.pt()!=0){
       if(l1PFTau.sumChargedIso()/l1PFTau.pt() < 0.40){
	 single_l1PFTauVLooseRelIso_.push_back(true);
       }else{
	 single_l1PFTauVLooseRelIso_.push_back(false);
       }
       if(l1PFTau.sumChargedIso()/l1PFTau.pt() < 0.20){
	 single_l1PFTauLooseRelIso_.push_back(true);
       }else{
	 single_l1PFTauLooseRelIso_.push_back(false);
       }
       if(l1PFTau.sumChargedIso()/l1PFTau.pt() < 0.10){
	 single_l1PFTauMediumRelIso_.push_back(true);
       }else{
	 single_l1PFTauMediumRelIso_.push_back(false);
       }
       if(l1PFTau.sumChargedIso()/l1PFTau.pt() < 0.05){
	 single_l1PFTauTightRelIso_.push_back(true);
       }else{
	 single_l1PFTauTightRelIso_.push_back(false);
       }
     }

     if(debug_){
       std::cout<<" L1PFTau pt "<<l1PFTau.pt()<<" eta "<< l1PFTau.eta()<<" phi "<< l1PFTau.phi()<<" charge "<< l1PFTau.charge()<<std::endl;
     }
   }

   for(auto l1TkMuon : *l1TkMuonHandle){
     if(fabs(l1TkMuon.pt())<min_muon_pt_)
       continue;
     if(fabs(l1TkMuon.eta())>max_eta_)
       continue;

     single_l1TkMuonPt_.push_back(l1TkMuon.pt());
     single_l1TkMuonEta_.push_back(l1TkMuon.eta());
     single_l1TkMuonPhi_.push_back(l1TkMuon.phi());
     single_l1TkMuonCharge_.push_back(l1TkMuon.charge());

     if(debug_){
       std::cout<<" L1TkMuon pt "<<l1TkMuon.pt()<<" eta "<< l1TkMuon.eta()<<" phi "<< l1TkMuon.phi()<<" charge "<< l1TkMuon.charge()<<std::endl;
     }

   }


   for(auto l1TkMuonPFTau : *l1TkMuonPFTauHandle){

     if(fabs(l1TkMuonPFTau.L1TkMuonParticle()->pt())<min_muon_pt_)
       continue;
     if(fabs(l1TkMuonPFTau.L1TkMuonParticle()->eta())>max_eta_)
       continue;
     if(fabs(l1TkMuonPFTau.L1HPSPFTau()->pt())<min_tau_pt_)
       continue;
     if(fabs(l1TkMuonPFTau.L1HPSPFTau()->eta())>max_eta_)
       continue;

     l1TkMuonPt_.push_back(l1TkMuonPFTau.L1TkMuonParticle()->pt());
     l1TkMuonEta_.push_back(l1TkMuonPFTau.L1TkMuonParticle()->eta());
     l1TkMuonPhi_.push_back(l1TkMuonPFTau.L1TkMuonParticle()->phi());
     l1TkMuonCharge_.push_back(l1TkMuonPFTau.L1TkMuonParticle()->charge());

     l1PFTauPt_.push_back(l1TkMuonPFTau.L1HPSPFTau()->pt());
     l1PFTauEta_.push_back(l1TkMuonPFTau.L1HPSPFTau()->eta());
     l1PFTauPhi_.push_back(l1TkMuonPFTau.L1HPSPFTau()->phi());
     l1PFTauCharge_.push_back(l1TkMuonPFTau.L1HPSPFTau()->charge());
     l1PFTauType_.push_back(l1TkMuonPFTau.L1HPSPFTau()->tauType());
     l1PFTauIso_.push_back(l1TkMuonPFTau.L1HPSPFTau()->sumChargedIso());
     l1PFTauTightIso_.push_back(l1TkMuonPFTau.L1HPSPFTau()->passTightIso());
     l1PFTauMediumIso_.push_back(l1TkMuonPFTau.L1HPSPFTau()->passMediumIso());
     l1PFTauLooseIso_.push_back(l1TkMuonPFTau.L1HPSPFTau()->passLooseIso());
     l1PFTauVLooseIso_.push_back(l1TkMuonPFTau.L1HPSPFTau()->passVLooseIso());
     double l1PFTauZ = 1000;
     if ( l1TkMuonPFTau.L1HPSPFTau()->leadChargedPFCand().isNonnull() && l1TkMuonPFTau.L1HPSPFTau()->leadChargedPFCand()->pfTrack().isNonnull())
       {
         l1PFTauZ = l1TkMuonPFTau.L1HPSPFTau()->leadChargedPFCand()->pfTrack()->vertex().z();
       }
     l1PFTauZ_.push_back(l1PFTauZ);

     if(l1TkMuonPFTau.L1HPSPFTau()->pt()!=0){
       if(l1TkMuonPFTau.L1HPSPFTau()->sumChargedIso()/l1TkMuonPFTau.L1HPSPFTau()->pt() < 0.40){
	 l1PFTauVLooseRelIso_.push_back(true);
       }
       else{
	 l1PFTauVLooseRelIso_.push_back(false);
       }
       if(l1TkMuonPFTau.L1HPSPFTau()->sumChargedIso()/l1TkMuonPFTau.L1HPSPFTau()->pt() < 0.20){
	 l1PFTauLooseRelIso_.push_back(true);
       }
       else{
	 l1PFTauLooseRelIso_.push_back(false);
       }
       if(l1TkMuonPFTau.L1HPSPFTau()->sumChargedIso()/l1TkMuonPFTau.L1HPSPFTau()->pt() < 0.10){
	 l1PFTauMediumRelIso_.push_back(true);
       }
       else{
	 l1PFTauMediumRelIso_.push_back(false);
       }
       if(l1TkMuonPFTau.L1HPSPFTau()->sumChargedIso()/l1TkMuonPFTau.L1HPSPFTau()->pt() < 0.05){
	 l1PFTauTightRelIso_.push_back(true);
       }
       else{
	 l1PFTauTightRelIso_.push_back(false);
       }
     }





     if(debug_){
       std::cout<<"l1Tk Muon Pt "<<l1TkMuonPFTau.L1TkMuonParticle()->pt()<<std::endl;
       std::cout<<"l1PF Tau Pt "<<l1TkMuonPFTau.L1HPSPFTau()->pt()<<std::endl;
     }
   }


   
   if(isReco_){
     edm::Handle<std::vector<reco::Vertex> > recoVertexHandle;
     iEvent.getByToken(recoVertexToken_,     recoVertexHandle);
     recoVertexN_ = recoVertexHandle->size();

     edm::Handle<std::vector<pat::Tau>>      recoTauHandle;
     iEvent.getByToken(recoTauToken_,        recoTauHandle);

     edm::Handle<pat::TauRefVector>          recoGMTauHandle;
     iEvent.getByToken(recoGMTauToken_,      recoGMTauHandle);

     for(auto recoTau : *recoTauHandle){
       if (fabs(recoTau.eta())>max_eta_)
	 continue;
       recoTauPt_.push_back(recoTau.pt());
       recoTauEta_.push_back(recoTau.eta());
       recoTauPhi_.push_back(recoTau.phi());
       recoTauCharge_.push_back(recoTau.charge());
       recoTauDecayMode_.push_back(recoTau.decayMode());
       bool isMatched = false;
       for(auto l1PFTau : *l1PFTauHandle){
	 double deltaEta = l1PFTau.eta() - recoTau.eta();
	 double deltaPhi = l1PFTau.phi() - recoTau.phi();
	 if ( (deltaEta*deltaEta + deltaPhi*deltaPhi) < 0.25 ){
	   isMatched = true;
	   hist_l1PFTauReso_vs_Reco_->Fill(l1PFTau.pt() / recoTau.pt());
	   break;
	 }
       }
       isRecoMatched_.push_back(isMatched);
       
       hist_isRecoMatched_->Fill(isMatched);
       hist_recoTauPt_->Fill(recoTau.pt());
       hist_recoTauEta_->Fill(recoTau.eta());
       hist_recoTauPhi_->Fill(recoTau.phi());
       
       if(debug_){
	 std::cout<<" RecoTau pt "<<recoTau.pt()<<" eta "<< recoTau.eta()<<" phi "<< recoTau.phi()<<" charge "<< recoTau.charge()<<" DecayMode "<< recoTau.decayMode()<<std::endl;
	 std::cout<<" RecoTau Z " <<recoTau.vertex().z() << std::endl;
       }
     } //for(auto recoTau : *recoTauHandle){

     for(auto recoGMTau : *recoGMTauHandle){
       if (fabs(recoGMTau->eta())>max_eta_)
	 continue;
       recoGMTauPt_.push_back(recoGMTau->pt());
       recoGMTauEta_.push_back(recoGMTau->eta());
       recoGMTauPhi_.push_back(recoGMTau->phi());
       recoGMTauCharge_.push_back(recoGMTau->charge());
       recoGMTauDecayMode_.push_back(recoGMTau->decayMode());
       bool isMatched = false;
       for(auto l1PFTau : *l1PFTauHandle){
	 double deltaEta = l1PFTau.eta() - recoGMTau->eta();
	 double deltaPhi = l1PFTau.phi() - recoGMTau->phi();
	 if ( (deltaEta*deltaEta + deltaPhi*deltaPhi) < 0.25 ){
	   isMatched = true;
	   hist_l1PFTauReso_vs_RecoGM_->Fill(l1PFTau.pt() / recoGMTau->pt());
	   break;
	 }
       }
       isRecoGMMatched_.push_back(isMatched);

       hist_isRecoGMMatched_->Fill(isMatched);
       hist_recoGMTauPt_->Fill(recoGMTau->pt());
       hist_recoGMTauEta_->Fill(recoGMTau->eta());
       hist_recoGMTauPhi_->Fill(recoGMTau->phi());

       if(debug_){
	 std::cout<<" RecoGMTau pt "<<recoGMTau->pt()<<" eta "<< recoGMTau->eta()<<" phi "<< recoGMTau->phi()<<" charge "<< recoGMTau->charge()<<" DecayMode "<< recoGMTau->decayMode()<<std::endl;
       }
     } //for(auto recoGMTau : *recoGMTauHandle){
   } //if(isReco_){







   tree_ -> Fill();
   treeBDT_ -> Fill();
}

void L1TkMuonHPSPFTauAnalyzer::Initialize() {
  indexevents_ = 0;
  runNumber_ = 0;
  lumi_ = 0;
  MC_weight_ = 1;
  genMuonPt_ .clear();
  genMuonEta_ .clear();
  genMuonPhi_ .clear();
  genMuonCharge_ .clear();
  genTauPt_ .clear();
  genTauEta_ .clear();
  genTauPhi_ .clear();
  genTauCharge_ .clear();
  isGenMatched_ .clear();
  recoTauPt_ .clear();
  recoTauEta_ .clear();
  recoTauPhi_ .clear();
  recoTauCharge_ .clear();
  isRecoMatched_ .clear();
  recoTauDecayMode_ .clear();
  recoGMTauPt_ .clear();
  recoGMTauEta_ .clear();
  recoGMTauPhi_ .clear();
  recoGMTauCharge_ .clear();
  isRecoGMMatched_ .clear();
  recoGMTauDecayMode_ .clear();
  l1TkMuonPt_ .clear();
  l1TkMuonEta_ .clear();
  l1TkMuonPhi_ .clear();
  l1TkMuonCharge_ .clear();
  l1PFTauPt_ .clear();
  l1PFTauEta_ .clear();
  l1PFTauPhi_ .clear();
  l1PFTauCharge_ .clear();
  l1PFTauType_ .clear();
  l1PFTauIso_ .clear();
  l1PFTauIsoOverPt_ .clear();
  l1PFTauTightIso_ .clear();
  l1PFTauMediumIso_ .clear();
  l1PFTauLooseIso_ .clear();
  l1PFTauVLooseIso_ .clear();
  l1PFTauTightRelIso_ .clear();
  l1PFTauMediumRelIso_ .clear();
  l1PFTauLooseRelIso_ .clear();
  l1PFTauVLooseRelIso_ .clear();
  l1PFTauLeadTrackPtOverTauPt_ .clear();
  l1PFTauChargedIso_ .clear();
  l1PFTauNeutralIso_ .clear();
  l1PFTauChargedIsoPileup_ .clear();
  l1PFTauRho_ .clear();
  l1PFTauNSignalChargedHadrons_ .clear();
  l1PFTauNSignalElectrons_ .clear();
  l1PFTauNSignalPhotons_ .clear();
  l1PFTauNSignalChargedPFCands_ .clear();
  l1PFTauSignalChargeSum_ .clear();
  l1PFTauStripPtOverTauPt_ .clear();
  l1PFTauStripMassOverTauPt_ .clear();
  l1PFTauStripMassOverStripPt_ .clear();
  l1PFTauZ_ .clear();
  l1PFTauStripPt_.clear();
  l1PFTauLeadTrackPt_.clear();
  l1PFTauVtxIndex_.clear();
  l1PFTaudz_.clear();
  l1PFTauSumTrackPtOfVtx_.clear();
  l1PFTauLeadTrackHoverE_.clear();
  l1PFTauHoverE_.clear();
  l1PFTauSignalTrackMass_.clear();
  l1PFTauNStripElectrons_.clear();
  l1PFTauNStripPhotons_.clear();
  l1PFTauDeltaRLeadTrackStrip_.clear();
  l1PFTauBDT_.clear();
  genVertex_ = 0;
  l1VertexN_ = 0;
  recoVertexN_ = 0;
  single_l1TkMuonPt_ .clear();
  single_l1TkMuonEta_ .clear();
  single_l1TkMuonPhi_ .clear();
  single_l1TkMuonCharge_ .clear();
  single_l1PFTauPt_ .clear();
  single_l1PFTauEta_ .clear();
  single_l1PFTauPhi_ .clear();
  single_l1PFTauCharge_ .clear();
  single_l1PFTauType_ .clear();
  single_l1PFTauIso_ .clear();
  single_l1PFTauTightIso_ .clear();
  single_l1PFTauMediumIso_ .clear();
  single_l1PFTauLooseIso_ .clear();
  single_l1PFTauVLooseIso_ .clear();
  single_l1PFTauTightRelIso_ .clear();
  single_l1PFTauMediumRelIso_ .clear();
  single_l1PFTauLooseRelIso_ .clear();
  single_l1PFTauVLooseRelIso_ .clear();
  single_l1PFTauZ_ .clear();
}


// ------------ method called once each job just before starting event loop  ------------
void
L1TkMuonHPSPFTauAnalyzer::beginJob()
{
  tree_ -> Branch("EventNumber",&indexevents_,"EventNumber/l");
  tree_ -> Branch("RunNumber",&runNumber_,"RunNumber/I");
  tree_ -> Branch("lumi",&lumi_,"lumi/I");
  tree_ -> Branch("MC_weight",&MC_weight_,"MC_weight/F");
  tree_ -> Branch("genMuonPt",  &genMuonPt_);
  tree_ -> Branch("genMuonEta", &genMuonEta_);
  tree_ -> Branch("genMuonPhi", &genMuonPhi_);
  tree_ -> Branch("genMuonCharge", &genMuonCharge_);
  tree_ -> Branch("genTauPt",  &genTauPt_);
  tree_ -> Branch("genTauEta", &genTauEta_);
  tree_ -> Branch("genTauPhi", &genTauPhi_);
  tree_ -> Branch("genTauCharge", &genTauCharge_);
  tree_ -> Branch("isGenMatched", &isGenMatched_);
  tree_ -> Branch("recoTauPt",  &recoTauPt_);
  tree_ -> Branch("recoTauEta", &recoTauEta_);
  tree_ -> Branch("recoTauPhi", &recoTauPhi_);
  tree_ -> Branch("recoTauCharge", &recoTauCharge_);
  tree_ -> Branch("isRecoMatched", &isRecoMatched_);
  tree_ -> Branch("recoTauDecayMode", &recoTauDecayMode_);
  tree_ -> Branch("recoGMTauPt",  &recoGMTauPt_);
  tree_ -> Branch("recoGMTauEta", &recoGMTauEta_);
  tree_ -> Branch("recoGMTauPhi", &recoGMTauPhi_);
  tree_ -> Branch("recoGMTauCharge", &recoGMTauCharge_);
  tree_ -> Branch("isRecoGMMatched", &isRecoGMMatched_);
  tree_ -> Branch("recoGMTauDecayMode", &recoGMTauDecayMode_);
  tree_ -> Branch("l1TkMuonPt",  &l1TkMuonPt_);
  tree_ -> Branch("l1TkMuonEta", &l1TkMuonEta_);
  tree_ -> Branch("l1TkMuonPhi", &l1TkMuonPhi_);
  tree_ -> Branch("l1TkMuonCharge", &l1TkMuonCharge_);
  tree_ -> Branch("l1PFTauPt",  &l1PFTauPt_);
  tree_ -> Branch("l1PFTauEta", &l1PFTauEta_);
  tree_ -> Branch("l1PFTauPhi", &l1PFTauPhi_);
  tree_ -> Branch("l1PFTauCharge", &l1PFTauCharge_);
  tree_ -> Branch("l1PFTauType", &l1PFTauType_);
  tree_ -> Branch("l1PFTauIso", &l1PFTauIso_);
  tree_ -> Branch("l1PFTauIsoOverPt", &l1PFTauIsoOverPt_);
  tree_ -> Branch("l1PFTauTightIso", &l1PFTauTightIso_);
  tree_ -> Branch("l1PFTauMediumIso", &l1PFTauMediumIso_);
  tree_ -> Branch("l1PFTauLooseIso", &l1PFTauLooseIso_);
  tree_ -> Branch("l1PFTauVLooseIso", &l1PFTauVLooseIso_);
  tree_ -> Branch("l1PFTauTightRelIso", &l1PFTauTightRelIso_);
  tree_ -> Branch("l1PFTauMediumRelIso", &l1PFTauMediumRelIso_);
  tree_ -> Branch("l1PFTauLooseRelIso", &l1PFTauLooseRelIso_);
  tree_ -> Branch("l1PFTauVLooseRelIso", &l1PFTauVLooseRelIso_);
  //tree_ -> Branch("l1PFTauBDT", &l1PFTauBDT_);
  tree_ -> Branch("l1PFTauZ", &l1PFTauZ_);
  //tree_ -> Branch("l1PFTaudz", &l1PFTaudz_ );
  tree_ -> Branch("genVertex", &genVertex_, "genVertex/D");
  tree_ -> Branch("l1VertexN", &l1VertexN_, "l1VertexN/I");
  tree_ -> Branch("recoVertexN", &recoVertexN_, "recoVertexN/I");
  tree_ -> Branch("single_l1TkMuonPt",  &single_l1TkMuonPt_);
  tree_ -> Branch("single_l1TkMuonEta", &single_l1TkMuonEta_);
  tree_ -> Branch("single_l1TkMuonPhi", &single_l1TkMuonPhi_);
  tree_ -> Branch("single_l1TkMuonCharge", &single_l1TkMuonCharge_);
  tree_ -> Branch("single_l1PFTauPt",  &single_l1PFTauPt_);
  tree_ -> Branch("single_l1PFTauEta", &single_l1PFTauEta_);
  tree_ -> Branch("single_l1PFTauPhi", &single_l1PFTauPhi_);
  tree_ -> Branch("single_l1PFTauCharge", &single_l1PFTauCharge_);
  tree_ -> Branch("single_l1PFTauType", &single_l1PFTauType_);
  tree_ -> Branch("single_l1PFTauIso", &single_l1PFTauIso_);
  tree_ -> Branch("single_l1PFTauTightIso", &single_l1PFTauTightIso_);
  tree_ -> Branch("single_l1PFTauMediumIso", &single_l1PFTauMediumIso_);
  tree_ -> Branch("single_l1PFTauLooseIso", &single_l1PFTauLooseIso_);
  tree_ -> Branch("single_l1PFTauVLooseIso", &single_l1PFTauVLooseIso_);
  tree_ -> Branch("single_l1PFTauTightRelIso", &single_l1PFTauTightRelIso_);
  tree_ -> Branch("single_l1PFTauMediumRelIso", &single_l1PFTauMediumRelIso_);
  tree_ -> Branch("single_l1PFTauLooseRelIso", &single_l1PFTauLooseRelIso_);
  tree_ -> Branch("single_l1PFTauVLooseRelIso", &single_l1PFTauVLooseRelIso_);
  tree_ -> Branch("single_l1PFTauZ", &single_l1PFTauZ_);

  hist_genMuonPt_ = new TH1F("genMuonPt","genMuonPt", 100, 0., 1000.);
  hist_genMuonEta_ = new TH1F("genMuonEta","genMuonEta",50, -3., 3.);
  hist_genMuonPhi_ = new TH1F("genMuonPhi","genMuonPhi",50, -3., 3.);
  hist_genTauPt_ = new TH1F("genTauPt","genTauPt", 100, 0., 1000.);
  hist_genTauEta_ = new TH1F("genTauEta","genTauEta",50, -3., 3.);
  hist_genTauPhi_ = new TH1F("genTauPhi","genTauPhi",50, -3., 3.);
  hist_isGenMatched_ = new TH1F("isGenMatched","isGenMatched", 3, -1., 2.);
  hist_recoTauPt_ = new TH1F("recoTauPt","recoTauPt", 100, 0., 1000.);
  hist_recoTauEta_ = new TH1F("recoTauEta","recoTauEta",50, -3., 3.);
  hist_recoTauPhi_ = new TH1F("recoTauPhi","recoTauPhi",50, -3., 3.);
  hist_isRecoMatched_ = new TH1F("isRecoMatched","isRecoMatched", 3, -1., 2.);
  hist_recoGMTauPt_ = new TH1F("recoGMTauPt","recoGMTauPt", 100, 0., 1000.);
  hist_recoGMTauEta_ = new TH1F("recoGMTauEta","recoGMTauEta",50, -3., 3.);
  hist_recoGMTauPhi_ = new TH1F("recoGMTauPhi","recoGMTauPhi",50, -3., 3.);
  hist_isRecoGMMatched_ = new TH1F("isRecoGMMatched","isRecoGMMatched", 3, -1., 2.);
  hist_l1TkMuonPt_ = new TH1F("l1TkMuonPt","l1TkMuonPt", 100, 0., 1000.);
  hist_l1TkMuonEta_ = new TH1F("l1TkMuonEta","l1TkMuonEta",50, -3., 3.);
  hist_l1TkMuonPhi_ = new TH1F("l1TkMuonPhi","l1TkMuonPhi",50, -3., 3.);
  hist_l1PFTauPt_ = new TH1F("l1PFTauPt","l1PFTauPt", 100, 0., 1000.);
  hist_l1PFTauEta_ = new TH1F("l1PFTauEta","l1PFTauEta",50, -3., 3.);
  hist_l1PFTauPhi_ = new TH1F("l1PFTauPhi","l1PFTauPhi",50, -3., 3.);
  hist_l1PFTauReso_vs_Gen_ = new TH1F("l1PFTauReso_vs_Gen","l1PFTauReso_vs_Gen", 60, 0., 3.);
  hist_l1PFTauReso_vs_Reco_ = new TH1F("l1PFTauReso_vs_Reco","l1PFTauReso_vs_Reco", 60, 0., 3.);
  hist_l1PFTauReso_vs_RecoGM_ = new TH1F("l1PFTauReso_vs_RecoGM","l1PFTauReso_vs_RecoGM", 60, 0., 3.);

  if(fillBDT_)
    {
      treeBDT_ -> Branch("MC_weight",&MC_weight_,"MC_weight/F");
      treeBDT_ -> Branch("genMuonPt",  &genMuonPt_);
      treeBDT_ -> Branch("genMuonEta", &genMuonEta_);
      treeBDT_ -> Branch("genMuonPhi", &genMuonPhi_);
      treeBDT_ -> Branch("genMuonCharge", &genMuonCharge_);
      treeBDT_ -> Branch("genTauPt",  &genTauPt_);
      treeBDT_ -> Branch("genTauEta", &genTauEta_);
      treeBDT_ -> Branch("genTauPhi", &genTauPhi_);
      treeBDT_ -> Branch("genTauCharge", &genTauCharge_);
      treeBDT_ -> Branch("l1TkMuonPt",   &l1TkMuonPt_);
      treeBDT_ -> Branch("l1TkMuonEta",  &l1TkMuonEta_);
      treeBDT_ -> Branch("l1TkMuonPhi", &l1TkMuonPhi_);
      treeBDT_ -> Branch("l1TkMuonCharge", &l1TkMuonCharge_);
      treeBDT_ -> Branch("l1PFTauPt",   &l1PFTauPt_);
      treeBDT_ -> Branch("l1PFTauEta",  &l1PFTauEta_);
      treeBDT_ -> Branch("l1PFTauPhi", &l1PFTauPhi_);
      treeBDT_ -> Branch("l1PFTauCharge", &l1PFTauCharge_);
      treeBDT_ -> Branch("l1PFTauType", &l1PFTauType_);
      treeBDT_ -> Branch("l1PFTauIso", &l1PFTauIso_);
      treeBDT_ -> Branch("l1PFTauZ", &l1PFTauZ_);
      treeBDT_ -> Branch("l1PFTauLeadTrackPtOverTauPt",  &l1PFTauLeadTrackPtOverTauPt_ );
      treeBDT_ -> Branch("l1PFTauChargedIso",  &l1PFTauChargedIso_ );
      treeBDT_ -> Branch("l1PFTauNeutralIso",  &l1PFTauNeutralIso_ );
      treeBDT_ -> Branch("l1PFTauChargedIsoPileup",  &l1PFTauChargedIsoPileup_ );
      treeBDT_ -> Branch("l1PFTauRho",  &l1PFTauRho_ );
      treeBDT_ -> Branch("l1PFTauNSignalChargedHadrons",  &l1PFTauNSignalChargedHadrons_ );
      treeBDT_ -> Branch("l1PFTauNSignalElectrons",  &l1PFTauNSignalElectrons_ );
      treeBDT_ -> Branch("l1PFTauNSignalPhotons",  &l1PFTauNSignalPhotons_ );
      treeBDT_ -> Branch("l1PFTauNSignalChargedPFCands",  &l1PFTauNSignalChargedPFCands_ );
      treeBDT_ -> Branch("l1PFTauSignalChargeSum",  &l1PFTauSignalChargeSum_ );
      treeBDT_ -> Branch("l1PFTauStripPtOverTauPt",  &l1PFTauStripPtOverTauPt_ );
      treeBDT_ -> Branch("l1PFTauStripMassOverTauPt",  &l1PFTauStripMassOverTauPt_ );
      treeBDT_ -> Branch("l1PFTauStripMassOverStripPt",  &l1PFTauStripMassOverStripPt_ );
      treeBDT_ -> Branch("l1PFTauStripPt", &l1PFTauStripPt_ );
      treeBDT_ -> Branch("l1PFTauLeadTrackPt", &l1PFTauLeadTrackPt_ );
      treeBDT_ -> Branch("l1PFTauVtxIndex", &l1PFTauVtxIndex_ );
      treeBDT_ -> Branch("l1PFTaudz", &l1PFTaudz_ );
      treeBDT_ -> Branch("l1PFTauSumTrackPtOfVtx", &l1PFTauSumTrackPtOfVtx_ );
      treeBDT_ -> Branch("l1PFTauLeadTrackHoverE", &l1PFTauLeadTrackHoverE_ );
      treeBDT_ -> Branch("l1PFTauHoverE", &l1PFTauHoverE_ );
      treeBDT_ -> Branch("l1PFTauSignalTrackMass", &l1PFTauSignalTrackMass_ );
      treeBDT_ -> Branch("l1PFTauNStripElectrons", &l1PFTauNStripElectrons_ );
      treeBDT_ -> Branch("l1PFTauNStripPhotons", &l1PFTauNStripPhotons_ );
      treeBDT_ -> Branch("l1PFTauDeltaRLeadTrackStrip", &l1PFTauDeltaRLeadTrackStrip_ );
    }

  if(applyBDT_){
    /*
    bdt_l1PFTauPt_ = 0.;
    bdt_l1PFTauEta_ = 0.;
    bdt_l1PFTauIso_ = 0.;
    bdt_l1PFTauNeutralIso_ = 0.;
    bdt_l1PFTauChargedIsoPileup_ = 0.;
    bdt_l1PFTauNSignalPhotons_ = 0.;
    bdt_l1PFTauSignalChargeSum_ = 0.;
    bdt_l1PFTauStripPt_ = 0.;
    bdt_l1PFTauLeadTrackPt_ = 0.;
    bdt_l1PFTaudz_ = 0.;
    bdt_l1PFTauLeadTrackHoverE_ = 0.;
    bdt_l1PFTauHoverE_ = 0.;
    bdt_l1PFTauSignalTrackMass_ = 0.;
    tmva_reader_ = new TMVA::Reader( "V:Color:!Silent" );

    tmva_reader_->AddVariable("l1PFTauPt", &bdt_l1PFTauPt_);
    tmva_reader_->AddVariable("l1PFTauEta", &bdt_l1PFTauEta_);
    tmva_reader_->AddVariable("l1PFTauIso", &bdt_l1PFTauIso_);
    tmva_reader_->AddVariable("l1PFTauNeutralIso", &bdt_l1PFTauNeutralIso_);
    tmva_reader_->AddVariable("l1PFTauChargedIsoPileup", &bdt_l1PFTauChargedIsoPileup_);
    tmva_reader_->AddVariable("l1PFTauNSignalPhotons", &bdt_l1PFTauNSignalPhotons_);
    tmva_reader_->AddVariable("l1PFTauSignalChargeSum", &bdt_l1PFTauSignalChargeSum_);
    tmva_reader_->AddVariable("l1PFTauStripPt", &bdt_l1PFTauStripPt_);
    tmva_reader_->AddVariable("l1PFTauLeadTrackPt", &bdt_l1PFTauLeadTrackPt_);
    tmva_reader_->AddVariable("l1PFTaudz", &bdt_l1PFTaudz_);
    tmva_reader_->AddVariable("l1PFTauLeadTrackHoverE", &bdt_l1PFTauLeadTrackHoverE_);
    //tmva_reader_->AddVariable("l1PFTauHoverE", &bdt_l1PFTauHoverE_);
    tmva_reader_->AddVariable("l1PFTauSignalTrackMass", &bdt_l1PFTauSignalTrackMass_);
    tmva_reader_->BookMVA("BDTG", "/home/sbhowmik/BDT/CMSSW_9_4_9/src/BDTTraining/L1HPS/test/L1HPSPFTau/L1HPSPFTau_XGB_testVars_default_12Var.xml");
    */
    std::vector<std::string> bdtInputVariables = {
      //"l1PFTauPt",
      "l1PFTauPhi",
      "l1PFTauIso",
      //"l1PFTauNeutralIso",
      "l1PFTauStripPt",
      "l1PFTauLeadTrackPt",
      "l1PFTaudz",
      "l1PFTauSignalTrackMass"
    };
    std::string fileName_tmva = bdtInputFileName_;
    bdt_tmva_L1HPSPFTau = new TMVAInterface(fileName_tmva, bdtInputVariables, { "iF_Recl[0]", "iF_Recl[1]", "iF_Recl[2]" });
    bdt_tmva_L1HPSPFTau->enableBDTTransform();
    std::string fileName_xgb = "L1TauAnalyzerPhase2/L1TauAnalyzerPhase2/data/L1HPSPFTau_XGB_testVars_default_12Var.pkl";
    //bdt_xgb_L1HPSPFTau = new XGBInterface(fileName_xgb, bdtInputVariables);
  }


  return;
}

// ------------ method called once each job just after ending the event loop  ------------
void
L1TkMuonHPSPFTauAnalyzer::endJob()
{
  if(fillBDT_)
    {
      bdtRootFile_->cd();
      treeBDT_->SaveAs();
      bdtRootFile_->Write();
      bdtRootFile_->Close();
    }
  if(createHistRoorFile_){
    histRootFile_->cd();

    hist_genMuonPt_->Write();
    hist_genMuonEta_->Write();
    hist_genMuonPhi_->Write();
    hist_genTauPt_->Write();
    hist_genTauEta_->Write();
    hist_genTauPhi_->Write();
    hist_isGenMatched_->Write();
    hist_recoTauPt_->Write();
    hist_recoTauEta_->Write();
    hist_recoTauPhi_->Write();
    hist_isRecoMatched_->Write();
    hist_recoGMTauPt_->Write();
    hist_recoGMTauEta_->Write();
    hist_recoGMTauPhi_->Write();
    hist_isRecoGMMatched_->Write();
    hist_l1TkMuonPt_->Write();
    hist_l1TkMuonEta_->Write();
    hist_l1TkMuonPhi_->Write();
    hist_l1PFTauPt_->Write();  
    hist_l1PFTauEta_->Write();
    hist_l1PFTauPhi_->Write();
    hist_l1PFTauReso_vs_Gen_->Write();
    hist_l1PFTauReso_vs_Reco_->Write();
    hist_l1PFTauReso_vs_RecoGM_->Write();

  //  histRootFile_->Write();
    histRootFile_->Close();
    delete bdt_tmva_L1HPSPFTau;
    delete bdt_xgb_L1HPSPFTau;
  }


  return;
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TkMuonHPSPFTauAnalyzer);
