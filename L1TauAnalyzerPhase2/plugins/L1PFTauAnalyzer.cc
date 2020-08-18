// Class:      L1PFTauAnalyzer
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
#include <TH2F.h>

#include "DataFormats/Phase2L1Taus/interface/L1HPSPFTau.h"
#include "DataFormats/Phase2L1Taus/interface/L1HPSPFTauFwd.h"
#include "DataFormats/L1Trigger/interface/L1PFTau.h"
#include "DataFormats/VertexReco/interface/Vertex.h"  
#include "DataFormats/L1Trigger/interface/Tau.h" 
#include "DataFormats/Phase2L1ParticleFlow/interface/PFTau.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include <DataFormats/PatCandidates/interface/Tau.h> 
#include "L1TauAnalyzerPhase2/L1TauAnalyzerPhase2/plugins/GenVertexProducer.h"
#include "TMVA/Reader.h"

class L1PFTauAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit L1PFTauAnalyzer(const edm::ParameterSet&);
  ~L1PFTauAnalyzer();
  
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
  std::vector<int> l1PFTauType_;
  std::vector<float> l1PFTauPt_;
  std::vector<float> l1PFTauEta_;
  std::vector<float> l1PFTauPhi_;
  std::vector<int> l1PFTauCharge_;
  std::vector<float> l1PFTauIso_;
  std::vector<Bool_t> l1PFTauTightIso_;
  std::vector<Bool_t> l1PFTauMediumIso_;
  std::vector<Bool_t> l1PFTauLooseIso_;
  std::vector<Bool_t> l1PFTauVLooseIso_;
  std::vector<Bool_t> l1PFTauTightRelIso_;
  std::vector<Bool_t> l1PFTauMediumRelIso_;
  std::vector<Bool_t> l1PFTauLooseRelIso_;
  std::vector<Bool_t> l1PFTauVLooseRelIso_;
  std::vector<float> l1PFTauZ_;

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
  TMVA::Reader *tmva_reader;

  bool createHistRoorFile_;
  std::string histRootFileName_;
  TFile* histRootFile_;
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
  TH1F* hist_l1PFTauPt_;
  TH1F* hist_l1PFTauEta_;
  TH1F* hist_l1PFTauPhi_;
  TH1F* hist_l1PFTauReso_vs_Gen_;
  TH1F* hist_l1PFTauReso_vs_Reco_;
  TH1F* hist_l1PFTauReso_vs_RecoGM_;
  TH1F* hist_l1PFTauIso_;
  TH2F* hist_l1PFTauIso_vs_l1PFTauPt_;

  // ----------member data ---------------------------

  bool debug_;
  bool isReco_;
  double min_pt_;
  double max_eta_;
  edm::EDGetTokenT<GenEventInfoProduct>            genTagToken_;
  edm::EDGetTokenT<double>                         genVertexToken_;
  edm::EDGetTokenT<std::vector<l1t::Vertex>>       l1VertexToken_;
  edm::EDGetTokenT<std::vector<reco::GenJet>>      genTauToken_;
  edm::EDGetTokenT<l1t::TauBxCollection>           l1TauToken_;
  //edm::EDGetTokenT<l1t::L1PFTauCollection>         l1PFTauToken_;
  edm::EDGetTokenT<l1t::PFTauCollection>           l1PFTauToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>>      recoVertexToken_;
  edm::EDGetTokenT<std::vector<pat::Tau>>          recoTauToken_;
  edm::EDGetTokenT<pat::TauRefVector>              recoGMTauToken_;

};


L1PFTauAnalyzer::L1PFTauAnalyzer(const edm::ParameterSet& iConfig)
  : debug_          (iConfig.getUntrackedParameter<bool>("debug", false))
  , isReco_         (iConfig.getUntrackedParameter<bool>("isReco", false))
  , min_pt_         (iConfig.getUntrackedParameter<double>("min_pt", 20))
  , max_eta_        (iConfig.getUntrackedParameter<double>("max_eta", 2.4))
  , genTagToken_    (consumes<GenEventInfoProduct>          (iConfig.getParameter<edm::InputTag>("genTagToken")))
  , genVertexToken_ (consumes<double>                       (iConfig.getParameter<edm::InputTag>("genVertexToken")))
  , l1VertexToken_  (consumes<std::vector<l1t::Vertex>>     (iConfig.getParameter<edm::InputTag>("l1VertexToken")))
  , genTauToken_    (consumes<std::vector<reco::GenJet>>    (iConfig.getParameter<edm::InputTag>("genTauToken")))
  , l1TauToken_     (consumes<l1t::TauBxCollection>         (iConfig.getParameter<edm::InputTag>("l1TauToken")))
  //, l1PFTauToken_   (consumes<l1t::L1PFTauCollection>       (iConfig.getParameter<edm::InputTag>("l1PFTauToken")))
  , l1PFTauToken_   (consumes<l1t::PFTauCollection>         (iConfig.getParameter<edm::InputTag>("l1PFTauToken")))
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
  return;
}


L1PFTauAnalyzer::~L1PFTauAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
L1PFTauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  Initialize();

  if(debug_){
    std::cout<<" Starting L1PFTau Analyzer ............     "<< std::endl;  
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

   edm::Handle<std::vector<reco::GenJet>>  genTauHandle;
   iEvent.getByToken(genTauToken_,         genTauHandle);

   edm::Handle< BXVector<l1t::Tau> >       l1TauHandle;
   iEvent.getByToken(l1TauToken_,          l1TauHandle);

   //edm::Handle<l1t::L1PFTauCollection>     l1PFTauHandle;
   edm::Handle<l1t::PFTauCollection>       l1PFTauHandle;
   iEvent.getByToken(l1PFTauToken_,        l1PFTauHandle);

   for(auto genTau : *genTauHandle){
     if (fabs(genTau.eta())>max_eta_)
       continue;
     //if (fabs(genTau.pt())<min_pt_)
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
     l1PFTauPt_.push_back(l1PFTau.pt());
     l1PFTauEta_.push_back(l1PFTau.eta());
     l1PFTauPhi_.push_back(l1PFTau.phi());
     l1PFTauCharge_.push_back(l1PFTau.charge());
     /*
     l1PFTauType_.push_back(l1PFTau.tauType());
     l1PFTauIso_.push_back(l1PFTau.chargedIso());
     l1PFTauTightIso_.push_back(l1PFTau.passTightIso());
     l1PFTauMediumIso_.push_back(l1PFTau.passMediumIso());
     l1PFTauLooseIso_.push_back(l1PFTau.passLooseIso());
     l1PFTauVLooseIso_.push_back(l1PFTau.passVLooseIso());
     */
     l1PFTauType_.push_back(1);
     l1PFTauIso_.push_back(l1PFTau.chargedIso());
     l1PFTauTightIso_.push_back(l1PFTau.passTightNN());
     l1PFTauMediumIso_.push_back(l1PFTau.passTightPF());
     l1PFTauLooseIso_.push_back(l1PFTau.passLooseNN());
     l1PFTauVLooseIso_.push_back(l1PFTau.passLoosePF());

     if(l1PFTau.pt()!=0)
       {
	 if(l1PFTau.chargedIso()/l1PFTau.pt() < 0.40)
	   {
	     l1PFTauVLooseRelIso_.push_back(true);
	   }
	 else
	   {
	     l1PFTauVLooseRelIso_.push_back(false);
	   }
	 if(l1PFTau.chargedIso()/l1PFTau.pt() < 0.20)
           {
             l1PFTauLooseRelIso_.push_back(true);
           }
         else
           {
             l1PFTauLooseRelIso_.push_back(false);
           }
	 if(l1PFTau.chargedIso()/l1PFTau.pt() < 0.10)
           {
             l1PFTauMediumRelIso_.push_back(true);
           }
         else
           {
             l1PFTauMediumRelIso_.push_back(false);
           }
	 if(l1PFTau.chargedIso()/l1PFTau.pt() < 0.05)
           {
             l1PFTauTightRelIso_.push_back(true);
           }
         else
           {
             l1PFTauTightRelIso_.push_back(false);
           }
       }
     double l1PFTauZ = 1000;
     ///l1PFTauZ = l1PFTau.z0();
     /*
     double leadTrackPt = 0;
     double leadTrackEta = 0;
     double leadTrackPhi = 0;
     if ( l1PFTau.leadChargedPFCand().isNonnull() && l1PFTau.leadChargedPFCand()->pfTrack().isNonnull())
       {
         l1PFTauZ = l1PFTau.leadChargedPFCand()->pfTrack()->vertex().z();
	 leadTrackPt = l1PFTau.leadChargedPFCand()->pfTrack()->pt();
	 leadTrackEta = l1PFTau.leadChargedPFCand()->pfTrack()->eta();
	 leadTrackPhi = l1PFTau.leadChargedPFCand()->pfTrack()->phi();
       }
     */
     l1PFTauZ_.push_back(l1PFTauZ);


     // define variables for BDT training
     /*
     l1PFTauLeadTrackPt_.push_back(leadTrackPt);
     l1PFTauLeadTrackPtOverTauPt_.push_back(leadTrackPt/l1PFTau.pt());
     double DeltaRLeadTrackStrip = sqrt(pow((leadTrackEta - l1PFTau.strip_p4().eta()),2) + pow((leadTrackPhi - l1PFTau.strip_p4().phi()),2));
     l1PFTauDeltaRLeadTrackStrip_.push_back(DeltaRLeadTrackStrip);
     double l1PFTauLeadTrackHoverE = 0.;
     if(l1PFTau.leadChargedPFCand().isNonnull() && l1PFTau.leadChargedPFCand()->pfCluster().isNonnull()){
       l1PFTauLeadTrackHoverE = l1PFTau.leadChargedPFCand()->pfCluster()->hOverE();
     }
     l1PFTauLeadTrackHoverE_.push_back(l1PFTauLeadTrackHoverE);
     int l1PFTauVtxIndex = -1;
     double l1PFTaudz = 10000.;
     for(unsigned int i=0; i<l1VertexHandle->size(); i++){
       double temp_dz = abs(l1PFTauZ - l1VertexHandle->at(i).z0());
       if(temp_dz < l1PFTaudz){
	 l1PFTaudz = temp_dz;
	 l1PFTauVtxIndex = i;
       }
     }
     l1PFTauVtxIndex_.push_back(l1PFTauVtxIndex);
     l1PFTaudz_.push_back(l1PFTaudz);
     l1PFTauChargedIso_.push_back(l1PFTau.chargedIso());
     l1PFTauNeutralIso_.push_back(l1PFTau.sumNeutralIso());
     l1PFTauChargedIsoPileup_.push_back(l1PFTau.chargedIsoPileup());
     l1PFTauRho_.push_back(l1PFTau.rhoCorr());
     l1PFTauNSignalChargedHadrons_.push_back(l1PFTau.signalChargedHadrons().size());
     l1PFTauNSignalElectrons_.push_back(l1PFTau.signalElectrons().size());
     l1PFTauNSignalPhotons_.push_back(l1PFTau.signalPhotons().size());
     int nSignalCan = l1PFTau.signalAllL1PFCandidates().size();
     int NSignalChargedPFCands = 0;
     int SumChargeSignalChargedPFCands = 0;
     reco::Candidate::LorentzVector l1PFTauSignalTrack_p4 ;  
     double SignalPFCandsTrack_Pt = 0;
     double sum_E = 0;
     double sum_H = 0;
     for (int i=0; i<nSignalCan; i++){
       if(l1PFTau.signalAllL1PFCandidates().at(i)->charge() != 0){
	 NSignalChargedPFCands++;
	 SumChargeSignalChargedPFCands += l1PFTau.signalAllL1PFCandidates().at(i)->charge();
       }
       if(l1PFTau.signalAllL1PFCandidates().at(i)->charge() != 0 && l1PFTau.signalAllL1PFCandidates().at(i)->pfTrack().isNonnull()){
	 l1PFTauSignalTrack_p4 += l1PFTau.signalAllL1PFCandidates().at(i)->pfTrack()->p4();
	 if((l1VertexHandle->size()!=0) && abs(l1PFTau.signalAllL1PFCandidates().at(i)->pfTrack()->vertex().z() - l1VertexHandle->at(l1PFTauVtxIndex).z0()) <= l1PFTaudz){
	   SignalPFCandsTrack_Pt += l1PFTau.signalAllL1PFCandidates().at(i)->pfTrack()->pt();
	 }
       }
       if(l1PFTau.signalAllL1PFCandidates().at(i)->pfCluster().isNonnull()){
	 double totalEnergy = l1PFTau.signalAllL1PFCandidates().at(i)->pfCluster()->energy();
	 double HoverE = l1PFTau.signalAllL1PFCandidates().at(i)->pfCluster()->hOverE();
	 double E = totalEnergy / (1+HoverE);
	 double H = HoverE * E;
	 sum_E += E;
	 sum_H += H;
       }
     }
     double l1PFTauHoverE = 0;
     if(sum_E!=0){
       l1PFTauHoverE = sum_H / sum_E ;
     }
     l1PFTauHoverE_.push_back(l1PFTauHoverE);
     l1PFTauSignalTrackMass_.push_back(l1PFTauSignalTrack_p4.mass());
     l1PFTauSumTrackPtOfVtx_.push_back(SignalPFCandsTrack_Pt);
     l1PFTauNSignalChargedPFCands_.push_back(NSignalChargedPFCands);
     l1PFTauSignalChargeSum_.push_back(abs(SumChargeSignalChargedPFCands));
     l1PFTauStripPt_.push_back(l1PFTau.strip_p4().pt());
     l1PFTauNStripElectrons_.push_back(l1PFTau.stripElectrons().size());
     l1PFTauNStripPhotons_.push_back(l1PFTau.stripPhotons().size());
     l1PFTauStripPtOverTauPt_.push_back(l1PFTau.strip_p4().pt()/l1PFTau.pt());
     l1PFTauStripMassOverTauPt_.push_back(l1PFTau.strip_p4().mass()/sqrt(l1PFTau.pt()));
     if(l1PFTau.strip_p4().pt()!=0)
       {
	 l1PFTauStripMassOverStripPt_.push_back(l1PFTau.strip_p4().mass()/sqrt(l1PFTau.strip_p4().pt()));
       }
     else
       {
	 l1PFTauStripMassOverStripPt_.push_back(0);
       }

     hist_l1PFTauPt_->Fill(l1PFTau.pt());
     hist_l1PFTauEta_->Fill(l1PFTau.eta());
     hist_l1PFTauPhi_->Fill(l1PFTau.phi());
     */
     hist_l1PFTauIso_->Fill(l1PFTau.chargedIso());
     hist_l1PFTauIso_vs_l1PFTauPt_->Fill(l1PFTau.pt(), l1PFTau.chargedIso());

     if(debug_){
       std::cout<<" L1PFTau pt "<<l1PFTau.pt()<<" eta "<< l1PFTau.eta()<<" phi "<< l1PFTau.phi()<<" charge "<< l1PFTau.charge()<<" Type "<< 1 <<" chargedIso "<< l1PFTau.chargedIso()<<std::endl;
       std::cout<<" L1PFTau Z " << l1PFTauZ  << std::endl;
     }
     /*
     // --- compute output of BDT 
     if(applyBDT_){
       float bdt_l1PFTauPt = 0.;
       float bdt_l1PFTauEta = 0.;
       float bdt_l1PFTauIso = 0.;
       float bdt_l1PFTauNeutralIso = 0.;
       float bdt_l1PFTauChargedIsoPileup = 0.;
       float bdt_l1PFTauNSignalPhotons = 0.;
       float bdt_l1PFTauSignalChargeSum = 0.;
       float bdt_l1PFTauStripPt = 0.;
       float bdt_l1PFTauLeadTrackPt = 0.;
       float bdt_l1PFTaudz = 0.;
       float bdt_l1PFTauLeadTrackHoverE = 0.;
       float bdt_l1PFTauHoverE = 0.;
       float bdt_l1PFTauSignalTrackMass = 0.;
       
       tmva_reader = new TMVA::Reader( "V:Color:!Silent" );

       tmva_reader->AddVariable("l1PFTauPt", &bdt_l1PFTauPt);
       tmva_reader->AddVariable("l1PFTauEta", &bdt_l1PFTauEta);
       tmva_reader->AddVariable("l1PFTauIso", &bdt_l1PFTauIso);
       tmva_reader->AddVariable("l1PFTauNeutralIso", &bdt_l1PFTauNeutralIso);
       tmva_reader->AddVariable("l1PFTauChargedIsoPileup", &bdt_l1PFTauChargedIsoPileup);
       tmva_reader->AddVariable("l1PFTauNSignalPhotons", &bdt_l1PFTauNSignalPhotons);
       tmva_reader->AddVariable("l1PFTauSignalChargeSum", &bdt_l1PFTauSignalChargeSum);
       tmva_reader->AddVariable("l1PFTauStripPt", &bdt_l1PFTauStripPt);
       tmva_reader->AddVariable("l1PFTauLeadTrackPt", &bdt_l1PFTauLeadTrackPt);
       tmva_reader->AddVariable("l1PFTaudz", &bdt_l1PFTaudz);
       tmva_reader->AddVariable("l1PFTauLeadTrackHoverE", &bdt_l1PFTauLeadTrackHoverE);
       tmva_reader->AddVariable("l1PFTauHoverE", &bdt_l1PFTauHoverE);
       tmva_reader->AddVariable("l1PFTauSignalTrackMass", &bdt_l1PFTauSignalTrackMass);
       
       tmva_reader->BookMVA("BDTG", "/home/sbhowmik/BDT/CMSSW_9_4_9/src/BDTTraining/L1HPS/test/L1HPSPFTau/L1HPSPFTau_XGB_testVars_default_13Var.xml");
       
       bdt_l1PFTauPt = l1PFTau.pt();
       bdt_l1PFTauEta = l1PFTau.eta();
       bdt_l1PFTauIso = l1PFTau.sumChargedIso();
       bdt_l1PFTauNeutralIso = l1PFTau.sumNeutralIso();
       bdt_l1PFTauChargedIsoPileup = l1PFTau.sumChargedIsoPileup();
       bdt_l1PFTauNSignalPhotons = l1PFTau.signalPhotons().size();
       bdt_l1PFTauSignalChargeSum = abs(SumChargeSignalChargedPFCands);
       bdt_l1PFTauStripPt = l1PFTau.strip_p4().pt();
       bdt_l1PFTauLeadTrackPt = leadTrackPt;
       bdt_l1PFTaudz = l1PFTaudz;
       bdt_l1PFTauLeadTrackHoverE = l1PFTauLeadTrackHoverE;
       bdt_l1PFTauHoverE = l1PFTauHoverE;
       bdt_l1PFTauSignalTrackMass = l1PFTauSignalTrack_p4.mass();
       
       Float_t BDT = tmva_reader->EvaluateMVA("BDTG");
       l1PFTauBDT_.push_back(BDT);
       if(debug_){
	 std::cout<<"BDT "<<BDT<<std::endl;
       }
     } //if(applyBDT){
     */
   } //for(auto l1PFTau : *l1PFTauHandle){


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

void L1PFTauAnalyzer::Initialize() {
  indexevents_ = 0;
  runNumber_ = 0;
  lumi_ = 0;
  MC_weight_ = 1;
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
  l1PFTauPt_ .clear();
  l1PFTauEta_ .clear();
  l1PFTauPhi_ .clear();
  l1PFTauCharge_ .clear();
  l1PFTauType_ .clear();
  l1PFTauIso_ .clear();
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
}


// ------------ method called once each job just before starting event loop  ------------
void
L1PFTauAnalyzer::beginJob()
{
  tree_ -> Branch("EventNumber",&indexevents_,"EventNumber/l");
  tree_ -> Branch("RunNumber",&runNumber_,"RunNumber/I");
  tree_ -> Branch("lumi",&lumi_,"lumi/I");
  tree_ -> Branch("MC_weight",&MC_weight_,"MC_weight/F");
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
  tree_ -> Branch("l1PFTauPt",  &l1PFTauPt_);
  tree_ -> Branch("l1PFTauEta", &l1PFTauEta_);
  tree_ -> Branch("l1PFTauPhi", &l1PFTauPhi_);
  tree_ -> Branch("l1PFTauCharge", &l1PFTauCharge_);
  tree_ -> Branch("l1PFTauType", &l1PFTauType_);
  tree_ -> Branch("l1PFTauIso", &l1PFTauIso_);
  tree_ -> Branch("l1PFTauTightIso", &l1PFTauTightIso_);
  tree_ -> Branch("l1PFTauMediumIso", &l1PFTauMediumIso_);
  tree_ -> Branch("l1PFTauLooseIso", &l1PFTauLooseIso_);
  tree_ -> Branch("l1PFTauVLooseIso", &l1PFTauVLooseIso_);
  tree_ -> Branch("l1PFTauTightRelIso", &l1PFTauTightRelIso_);
  tree_ -> Branch("l1PFTauMediumRelIso", &l1PFTauMediumRelIso_);
  tree_ -> Branch("l1PFTauLooseRelIso", &l1PFTauLooseRelIso_);
  tree_ -> Branch("l1PFTauVLooseRelIso", &l1PFTauVLooseRelIso_);
  tree_ -> Branch("l1PFTauBDT", &l1PFTauBDT_);
  tree_ -> Branch("l1PFTauZ", &l1PFTauZ_);
  tree_ -> Branch("l1PFTaudz", &l1PFTaudz_ );
  tree_ -> Branch("genVertex", &genVertex_, "genVertex/D");
  tree_ -> Branch("l1VertexN", &l1VertexN_, "l1VertexN/I");
  tree_ -> Branch("recoVertexN", &recoVertexN_, "recoVertexN/I");

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
  hist_l1PFTauPt_ = new TH1F("l1PFTauPt","l1PFTauPt", 1000, 0., 1000.);
  hist_l1PFTauEta_ = new TH1F("l1PFTauEta","l1PFTauEta",50, -3., 3.);
  hist_l1PFTauPhi_ = new TH1F("l1PFTauPhi","l1PFTauPhi",50, -3., 3.);
  hist_l1PFTauReso_vs_Gen_ = new TH1F("l1PFTauReso_vs_Gen","l1PFTauReso_vs_Gen", 60, 0., 3.);
  hist_l1PFTauReso_vs_Reco_ = new TH1F("l1PFTauReso_vs_Reco","l1PFTauReso_vs_Reco", 60, 0., 3.);
  hist_l1PFTauReso_vs_RecoGM_ = new TH1F("l1PFTauReso_vs_RecoGM","l1PFTauReso_vs_RecoGM", 60, 0., 3.);
  hist_l1PFTauIso_ = new TH1F("l1PFTauIso","l1PFTauIso", 100, 0., 1.);
  hist_l1PFTauIso_vs_l1PFTauPt_ = new TH2F("l1PFTauIso_vs_l1PFTauPt","l1PFTauIso_vs_l1PFTauPt", 100, 0., 500., 100, 0., 1.);


  if(fillBDT_)
    {
      treeBDT_ -> Branch("MC_weight",&MC_weight_,"MC_weight/F");
      treeBDT_ -> Branch("genTauPt",  &genTauPt_);
      treeBDT_ -> Branch("genTauEta", &genTauEta_);
      treeBDT_ -> Branch("genTauPhi", &genTauPhi_);
      treeBDT_ -> Branch("genTauCharge", &genTauCharge_);
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

  return;
}

// ------------ method called once each job just after ending the event loop  ------------
void
L1PFTauAnalyzer::endJob()
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
    hist_l1PFTauPt_->Write();  
    hist_l1PFTauEta_->Write();
    hist_l1PFTauPhi_->Write();
    hist_l1PFTauReso_vs_Gen_->Write();
    hist_l1PFTauReso_vs_Reco_->Write();
    hist_l1PFTauReso_vs_RecoGM_->Write();
    hist_l1PFTauIso_->Write();
    hist_l1PFTauIso_vs_l1PFTauPt_->Write();
  //  histRootFile_->Write();
    histRootFile_->Close();
  }


  return;
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1PFTauAnalyzer);
