#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
//#include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h"
#include "SusyAnaTools/Skims/plugins/ElectronEffectiveArea.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

#include "DataFormats/METReco/interface/MET.h"

#include "TLorentzVector.h"

typedef std::vector< edm::Handle< edm::ValueMap<reco::IsoDeposit> > >   IsoDepositMaps;
typedef std::vector< edm::Handle< edm::ValueMap<double> > >             IsoDepositVals;

class prodElectrons : public edm::EDFilter {

  public:

    explicit prodElectrons(const edm::ParameterSet & iConfig);
    ~prodElectrons();

  private:

    virtual bool filter(edm::Event & iEvent, const edm::EventSetup & iSetup);

    edm::InputTag electronSrc_;
    edm::InputTag conversionsSrc_;
    edm::InputTag vtxSrc_;
    edm::InputTag metSrc_;
    edm::InputTag beamSpotSrc_;
    bool   doEleVeto_, doEleIso_;
    double minElePt_, maxEleEta_;
    bool debug_;

};


typedef std::vector< edm::Handle< edm::ValueMap<reco::IsoDeposit> > >   IsoDepositMaps;
typedef std::vector< edm::Handle< edm::ValueMap<double> > >             IsoDepositVals;


prodElectrons::prodElectrons(const edm::ParameterSet & iConfig) {
  electronSrc_   = iConfig.getParameter<edm::InputTag>("ElectronSource");
  conversionsSrc_= iConfig.getParameter<edm::InputTag>("ConversionsSource");
  vtxSrc_        = iConfig.getParameter<edm::InputTag>("VertexSource");
  metSrc_        = iConfig.getParameter<edm::InputTag>("metSource");
  beamSpotSrc_   = iConfig.getParameter<edm::InputTag>("BeamSpotSource");
  minElePt_      = iConfig.getParameter<double>("MinElePt");
  maxEleEta_     = iConfig.getParameter<double>("MaxEleEta");
  doEleVeto_     = iConfig.getParameter<bool>("DoElectronVeto");
  doEleIso_      = iConfig.getParameter<bool>("DoElectronIsolation");
  debug_         = iConfig.getParameter<bool>("Debug");

  produces<std::vector<pat::Electron> >("");
  produces<std::vector<TLorentzVector> >("elesLVec");
  produces<std::vector<double> >("elesCharge");
  produces<std::vector<double> >("elesMtw");
  produces<std::vector<double> >("elesRelIso");
  produces<std::vector<bool> >("elesisEB");
  produces<int>("nElectrons");
}


prodElectrons::~prodElectrons() {
}


bool prodElectrons::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // electrons
  edm::Handle< edm::View<pat::Electron> > electrons;   
  iEvent.getByLabel(electronSrc_, electrons);

  // conversions
  edm::Handle< std::vector<reco::Conversion> > conversions;
  iEvent.getByLabel(conversionsSrc_, conversions);

  // beam spot
  edm::Handle<reco::BeamSpot> beamspot;
  iEvent.getByLabel(beamSpotSrc_, beamspot);
//  const reco::BeamSpot &beamSpot = *(beamspot.product());
  
  // vertices
  edm::Handle< std::vector<reco::Vertex> > vertices;
  iEvent.getByLabel(vtxSrc_, vertices);

  edm::Handle<edm::View<reco::MET> > met;
  iEvent.getByLabel(metSrc_, met);

  float cut_sigmaIEtaIEta[2]  = {999.9, 999.9};
  float cut_dEtaIn[2]         = {999.9, 999.9};
  float cut_dPhiIn[2]         = {999.9, 999.9};
  float cut_hoe[2]            = {999.9, 999.9};
  float cut_iso[2]            = {999.9, 999.9};
  float cut_ooemoop[2]        = {999.9, 999.9};
  float cut_d0vtx[2]          = {999.9, 999.9};
  float cut_dzvtx[2]          = {999.9, 999.9};
  unsigned int cut_mHits[2]   = {999, 999};
  bool cut_convVeto[2]          = {false, false};
    
  cut_sigmaIEtaIEta[0] = 0.011100; cut_sigmaIEtaIEta[1] = 0.033987;
  cut_dEtaIn[0]        = 0.016315; cut_dEtaIn[1]        = 0.010671;
  cut_dPhiIn[0]        = 0.252044; cut_dPhiIn[1]        = 0.245263;
  cut_hoe[0]           = 0.345843; cut_hoe[1]           = 0.134691;
  cut_iso[0]           = 0.164369; cut_iso[1]           = 0.212604;
  cut_ooemoop[0]       = 0.248070; cut_ooemoop[1]       = 0.157160;
  cut_d0vtx[0]         = 0.060279; cut_d0vtx[1]         = 0.273097;
  cut_dzvtx[0]         = 0.800538; cut_dzvtx[1]         = 0.885860;
  cut_mHits[0]         = 2;        cut_mHits[1]         = 3;
  cut_convVeto[0]      = true;     cut_convVeto[1]      = true;

  // check which ones to keep
  std::auto_ptr<std::vector<pat::Electron> > prod(new std::vector<pat::Electron>());
  std::auto_ptr<std::vector<TLorentzVector> > elesLVec(new std::vector<TLorentzVector>());
  std::auto_ptr<std::vector<double> > elesCharge(new std::vector<double>());
  std::auto_ptr<std::vector<double> > elesMtw(new std::vector<double>());
  std::auto_ptr<std::vector<double> > elesRelIso(new std::vector<double>());
  std::auto_ptr<std::vector<bool> > elesisEB(new std::vector<bool>());

  // loop on electrons
  for( edm::View<pat::Electron>::const_iterator ele = electrons->begin(); ele != electrons->end(); ele++ ){

    double pt = ele->pt();
    if (ele->pt() < minElePt_) continue;

    // get the ID variables from the electron object
    // kinematic variables
    bool isEB           = ele->isEB() ? true : false;
//    float eta           = ele->superCluster()->eta();

    // id variables
    float sigmaIEtaIEta = ele->full5x5_sigmaIetaIeta();
    float dEtaIn        = ele->deltaEtaSuperClusterTrackAtVtx();
    float dPhiIn        = ele->deltaPhiSuperClusterTrackAtVtx();
    float hoe           = ele->hadronicOverEm();
    float ooemoop       = 1e30;
    if( ele->ecalEnergy() !=0 && std::isfinite(ele->ecalEnergy()) ){
       ooemoop = std::abs(1.0/ele->ecalEnergy() - ele->eSuperClusterOverP()/ele->ecalEnergy());
    }

    // impact parameter variables
    float d0vtx         = 0.0;
    float dzvtx         = 0.0;
    if (vertices->size() > 0) {
        reco::VertexRef vtx(vertices, 0);    
        d0vtx = ele->gsfTrack()->dxy(vtx->position());
        dzvtx = ele->gsfTrack()->dz(vtx->position());
    } else {
        d0vtx = ele->gsfTrack()->dxy();
        dzvtx = ele->gsfTrack()->dz();
    }

    // conversion rejection variables
    bool convVeto = ele->passConversionVeto();
    float mHits = ele->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
    
    // choose cut if barrel or endcap
    unsigned int idx = isEB ? 0 : 1;

    // test cuts
    if (sigmaIEtaIEta >= cut_sigmaIEtaIEta[idx])     continue;
    if (fabs(dEtaIn) >= cut_dEtaIn[idx])             continue;
    if (fabs(dPhiIn) >= cut_dPhiIn[idx])             continue;
    if (hoe >= cut_hoe[idx])                         continue;
    if (fabs(ooemoop) >= cut_ooemoop[idx])           continue;
    if (fabs(d0vtx) >= cut_d0vtx[idx])               continue;
    if (fabs(dzvtx) >= cut_dzvtx[idx])               continue;
    if (mHits > cut_mHits[idx])                      continue;
    if (convVeto != cut_convVeto[idx])               continue;

    if(debug_) {
      reco::VertexRef vtx(vertices, 0);
      std::cout << "iEle " << ele - electrons->begin()  << ": "
		<< " (pt,eta,phi) "<<ele->pt()<<", "<<ele->eta()<<", "<<ele->phi() << " "
		<< ", isEB " << ele->isEB() << ", isEE " << ele->isEE() << "\n"
		<< ", dEtaIn " << ele->deltaEtaSuperClusterTrackAtVtx()
		<< ", dPhiIn " << ele->deltaPhiSuperClusterTrackAtVtx()
		<< ", sigmaIEtaIEta "<< ele->sigmaIetaIeta()
		<< ", hoe " << ele->hadronicOverEm()
		<< ", d0vtx " << ele->gsfTrack()->dxy(vtx->position())
		<< ", dzvtx " << ele->gsfTrack()->dz(vtx->position())
		<< ", passSelection " 
		<< std::endl;
    }


    // isolation cuts                                                                                                                                        
    reco::GsfElectron::PflowIsolationVariables pfIso = ele->pfIsolationVariables();
    float absiso = pfIso.sumChargedHadronPt + std::max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt );

    // compute final isolation
    double iso = absiso/pt;

    if (doEleIso_) {
      if(iso >= cut_iso[idx]) continue;
    }

    // electron is ID'd and isolated! - only accept if vertex present
    if (vertices->size()>0){
       prod->push_back(pat::Electron(*ele));
       TLorentzVector perLVec; perLVec.SetPtEtaPhiE(ele->pt(), ele->eta(), ele->phi(), ele->energy());
       elesLVec->push_back(perLVec);

       double mtw = sqrt( 2*( (*met)[0].pt()*ele->pt() -( (*met)[0].px()*ele->px() + (*met)[0].py()*ele->py() ) ) );

       elesCharge->push_back(ele->charge());
       elesMtw->push_back(mtw);
       elesRelIso->push_back(iso);
       elesisEB->push_back(isEB);
    }
  }


  // determine result before losing ownership of the pointer
  bool result = (doEleVeto_ ? (prod->size() == 0) : true);

  std::auto_ptr<int> nElectrons (new int);

  *nElectrons = prod->size();

  // store in the event
  iEvent.put(prod);
  iEvent.put(elesLVec, "elesLVec");
  iEvent.put(elesCharge, "elesCharge");
  iEvent.put(elesMtw, "elesMtw");
  iEvent.put(elesRelIso, "elesRelIso");
  iEvent.put(elesisEB, "elesisEB");
  iEvent.put(nElectrons, "nElectrons");

  return result;
}


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(prodElectrons);
