#include <memory>
#include <algorithm>
#include <assert.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/METReco/interface/MET.h"

#include "TLorentzVector.h"

#include "StopTupleMaker/SkimsAUX/plugins/common.h"

class prodMuons : public edm::EDFilter
{
 public:
  explicit prodMuons(const edm::ParameterSet & iConfig);
  ~prodMuons();
 private:
  virtual bool filter(edm::Event & iEvent, const edm::EventSetup & iSetup);

  edm::InputTag muonSrc_;
  edm::InputTag vtxSrc_;
  edm::InputTag metSrc_;
  edm::InputTag pfCandsSrc_;
  edm::InputTag rhoSrc_;
  edm::EDGetTokenT<std::vector<pat::Muon>> MuonTok_;
  edm::EDGetTokenT< std::vector<reco::Vertex> > VtxTok_;
  edm::EDGetTokenT<edm::View<reco::MET> > MetTok_;
  edm::EDGetTokenT<pat::PackedCandidateCollection>  PfcandTok_;
  edm::EDGetTokenT<double> RhoTok_;
  bool debug_;
  //bool doMuonVeto_;//, doMuonID_, doMuonVtx_;
//  int doMuonIso_; // 0: don't do any isolation; 1: relIso;  2: miniIso
  double minMuPt_, maxMuEta_, maxMuD0_, maxMuDz_, maxMuRelIso_, maxMuMiniIso_, minMuNumHit_;
//  double minMuPtForMuon2Clean_;
//  bool specialFix_;
//  edm::InputTag badGlobalMuonTaggerSrc_, cloneGlobalMuonTaggerSrc_;
//  edm::EDGetTokenT<edm::PtrVector<reco::Muon>> badGlobalMuonTok_, cloneGlobalMuonTok_;

  bool isLooseMuon(const pat::Muon & muon);
  bool isMediumMuon(const pat::Muon & muon );
  bool isMediumPlusMuon(const pat::Muon & muon, const reco::Vertex::Point & vtxpos);
  bool isTightMuon(const pat::Muon & muon, const reco::Vertex::Point & vtxpos);
//  bool isTightMuonOld(const pat::Muon & muon, const reco::Vertex::Point & vtxpos);
};


prodMuons::prodMuons(const edm::ParameterSet & iConfig) 
{
  muonSrc_      = iConfig.getParameter<edm::InputTag>("MuonSource");
  vtxSrc_       = iConfig.getParameter<edm::InputTag>("VertexSource");
  metSrc_       = iConfig.getParameter<edm::InputTag>("metSource");
  pfCandsSrc_   = iConfig.getParameter<edm::InputTag>("PFCandSource");
  rhoSrc_       = iConfig.getParameter<edm::InputTag>("RhoSource");
  minMuPt_      = iConfig.getParameter<double>("MinMuPt");
  maxMuEta_     = iConfig.getParameter<double>("MaxMuEta");
  maxMuD0_      = iConfig.getParameter<double>("MaxMuD0");
  maxMuDz_      = iConfig.getParameter<double>("MaxMuDz");
  maxMuRelIso_  = iConfig.getParameter<double>("MaxMuRelIso");
  maxMuMiniIso_ = iConfig.getParameter<double>("MaxMuMiniIso");
  minMuNumHit_  = iConfig.getParameter<double>("MinMuNumHit");
  doMuonVeto_   = iConfig.getParameter<bool>("DoMuonVeto"); //?
//  doMuonID_     = iConfig.getParameter<bool>("DoMuonID");
//  doMuonVtx_    = iConfig.getParameter<bool>("DoMuonVtxAssociation");
//  doMuonIso_    = iConfig.getParameter<int>("DoMuonIsolation");
  debug_        = iConfig.getParameter<bool>("Debug");
//  specialFix_   = iConfig.getParameter<bool>("specialFix");

  MuonTok_ = consumes<std::vector<pat::Muon>>(muonSrc_);
  VtxTok_ = consumes<std::vector<reco::Vertex>>(vtxSrc_);
  MetTok_ = consumes<edm::View<reco::MET>>(metSrc_);
  PfcandTok_ = consumes<pat::PackedCandidateCollection>(pfCandsSrc_);
  RhoTok_ = consumes<double>(rhoSrc_);

  produces<std::vector<pat::Muon> >("");
//  produces<std::vector<pat::Muon> >("mu2Clean");
  //0: loose;  1: medium;  1.5: mediumPlus 2: tight
  produces<std::vector<float> >("muonsIDtype");
  produces<std::vector<int> >("muonsFlagLoose");
  produces<std::vector<int> >("muonsFlagMedium");
  produces<std::vector<int> >("muonsFlagMediumPlus");
  produces<std::vector<int> >("muonsFlagTight");
  produces<std::vector<TLorentzVector> >("muonsLVec");
  produces<std::vector<float> >("muonsCharge");
  produces<std::vector<float> >("muonsMtw");
  produces<std::vector<float> >("muonsRelIso");
  produces<std::vector<float> >("muonsMiniIso");
  produces<std::vector<float> >("muonspfActivity");
  //here add muon vtx info /////////////////////////////////////////////??TODO  ADD THE IMPACT STUFF HERE!!!!!

  produces<std::vector<float> >("muonsDB_PV2D");
  produces<std::vector<float> >("muonsDB_PV3D");
  produces<std::vector<float> >("muonsDB_PVDZ");
  produces<std::vector<float> >("muonsDB_BS2D");
  produces<std::vector<float> >("muonsDB_BS3D");
  produces<std::vector<float> >("muonseDB_PV2D");
  produces<std::vector<float> >("muonseDB_PV3D");
  produces<std::vector<float> >("muonseDB_PVDZ");
  produces<std::vector<float> >("muonseDB_BS2D");
  produces<std::vector<float> >("muonseDB_BS3D");
  produces<std::vector<float> >("muonsBT_Dz");
  produces<std::vector<float> >("muonsBT_Dxy");
  produces<std::vector<float> >("muonsBTe_Dxy");
  produces<std::vector<float> >("muonsBTe_Dz");
  produces<int>("nMuons");
}//end constructor

prodMuons::~prodMuons()
{
}//end destructor


//workhorse
bool prodMuons::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	 //read in the objects
  
  edm::Handle<std::vector<pat::Muon> > muons;
  iEvent.getByToken(MuonTok_, muons);
  
  edm::Handle< std::vector<reco::Vertex> > vertices;
  iEvent.getByToken(VtxTok_, vertices);

  reco::Vertex::Point vtxpos = (vertices->size() > 0 ? (*vertices)[0].position() : reco::Vertex::Point());
  edm::Handle<edm::View<reco::MET> > met;
  iEvent.getByToken(MetTok_, met);

  edm::Handle<pat::PackedCandidateCollection> pfcands;
  iEvent.getByToken(PfcandTok_, pfcands);

  edm::Handle< double > rho_;
  iEvent.getByToken(RhoTok_,rho_); 
  //iEvent.getByToken("fixedGridRhoFastjetCentralNeutral", rho_); // Central rho recommended for SUSY
  double rho = *rho_;

	 //check which ones to keep
  std::unique_ptr<std::vector<pat::Muon> > prod(new std::vector<pat::Muon>()); //?
 // std::unique_ptr<std::vector<pat::Muon> > mu2Clean(new std::vector<pat::Muon>());
  std::unique_ptr<std::vector<TLorentzVector> > muonsLVec(new std::vector<TLorentzVector>());
  std::unique_ptr<std::vector<float> > muonsCharge(new std::vector<float>());
  std::unique_ptr<std::vector<float> > muonsMtw(new std::vector<float>());
  std::unique_ptr<std::vector<float> > muonsRelIso(new std::vector<float>());
  std::unique_ptr<std::vector<float> > muonsMiniIso(new std::vector<float>());
  std::unique_ptr<std::vector<float> > muonspfActivity(new std::vector<float>());

  std::unique_ptr<std::vector<float> > muonsIDtype(new std::vector<float>());
  std::unique_ptr<std::vector<int> > muonsFlagLoose(new std::vector<int>());
  std::unique_ptr<std::vector<int> > muonsFlagMedium(new std::vector<int>());
  std::unique_ptr<std::vector<int> > muonsFlagMediumPlus( new std::vector<int>());
  std::unique_ptr<std::vector<int> > muonsFlagTight(new std::vector<int>());
  std::unique_ptr<int> nMuons(new int);

  //vertex info?
  std::unique_ptr<std::vector<float> > muonsDB_PV2D( new std::vector<float>());
  std::unique_ptr<std::vector<float> > muonsDB_PV3D( new std::vector<float>());
  std::unique_ptr<std::vector<float> > muonsDB_PVDZ( new std::vector<float>());
  std::unique_ptr<std::vector<float> > muonsDB_BS2D( new std::vector<float>());
  std::unique_ptr<std::vector<float> > muosnDB_BS3D( new std::vector<float>());
  std::unique_ptr<std::vector<float> > muonseDB_PV2D( new std::vector<float>());
  std::unique_ptr<std::vector<float> > muonseDB_PV3D( new std::vector<float>());
  std::unique_ptr<std::vector<float> > muonseDB_PVDZ( new std::vector<float>());
  std::unique_ptr<std::vector<float> > muonseDB_BS2D( new std::vector<float>());
  std::unique_ptr<std::vector<float> > muonseDB_BS3D( new std::vector<float>());

  std::unique_ptr<std::vector<float> > muonsBT_Dz( new std::vector<float>());
  std::unique_ptr<std::vector<float> > muonsBT_Dxy( new std::vector<float>());
  std::unique_ptr<std::vector<float> > muonsBTe_Dxy( new std::vector<float>());
  std::unique_ptr<std::vector<float> > muonsBTe_Dz( new std::vector<float>());

   if (vertices->size() > 0) 
  {
    //for (std::vector<int> m = muons->begin(); m != muons->end(); ++m) {
    for (std::vector<pat::Muon>::const_iterator m = muons->begin(); m != muons->end(); ++m)
    {
      //acceptance cuts from config
      if (m->pt() < minMuPt_) continue;
      if (std::abs(m->eta()) >= maxMuEta_) continue;

      bool isLooseID = isLooseMuon((*m));
	  bool isMediumID = isMediumMuon((*m));
      bool isMediumPlusID = isMediumPlusMuon((*m), vtxpos);
      bool isTightID = isTightMuon((*m), vtxpos);

      // only store muons passing medium or tight ID
     // if ( !(isMediumID || isTightID) ) continue;
	 //store ALL Muons!

      // isolation cuts
      double muRelIso = getRelIso(*m); 
      double miniIso = commonFunctions::GetMiniIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&(*m)), "muon", rho);
      double pfActivity = commonFunctions::GetMiniIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&(*m)), "muon", rho, true);

  
	  //isolation cuts from config
      if (muRelIso >= maxMuRelIso_) continue;
      if(debug_) {std::cout << "PassedMuon Isolation" << std::endl;}
      if (miniIso >= maxMuMiniIso_ ) continue;
   

      // muon is ID'd and isolated! - only accept if vertex present
   //   if (vertices->size() > 0)  //this redundant?
   //   {
        prod->push_back(pat::Muon(*m));
        TLorentzVector perLVec; 
		perLVec.SetPtEtaPhiE(m->pt(), m->eta(), m->phi(), m->energy());
        muonsLVec->push_back(perLVec);

        double mtw = getMtw(*m, *met);

        muonsCharge->push_back(m->charge());
        muonsMtw->push_back(mtw);
        muonsRelIso->push_back(muRelIso);
        muonsMiniIso->push_back(miniIso);
		muonspfActivity->push_back(pfActivity);

        float idtype =-1;
		if( isLooseId ) idtype = 0;
        if( isMediumID ) idtype = 1;
		if( isMediumPlusID ) idtype = 1.5;
        if( isTightID ) idtype = 2;

        muonsIDtype->push_back(idtype); 

        if( isLooseID ) muonsFlagLoose->push_back(1); else muonsFlagLoose->push_back(0);
        if( isMediumID ) muonsFlagMedium->push_back(1); else muonsFlagMedium->push_back(0);
		if( isMediumPlusID ) muonsFlagMediumPlus->push_back(1); else muonsFlagMediumPlus->push_back(0);
        if( isTightID ) muonsFlagTight->push_back(1); else muonsFlagTight->push_back(0);



		

		//need to fill vertex info?
		muonsDB_PV2D->push_back( m->dB(pat::Muon::PV2D) );
		muonsDB_PV3D->push_back( m->dB(pat::Muon::PV3D) );
		muonsDB_PVDZ->push_back( m->dB(pat::Muon::PVDZ) );
		muonsDB_BS2D->push_back( m->dB(pat::Muon::BS2D) );
		muonsDB_BS3D->push_back( m->dB(pat::Muon::BS3D) );
		
		muonseDB_PV2D->push_back( m->edB(pat::Muon::PV2D) );
		muonseDB_PV3D->push_back( m->edB(pat::Muon::PV3D) );
		muonseDB_PVDZ->push_back( m->edB(pat::Muon::PVDZ) );
		muonseDB_BS2D->push_back( m->edB(pat::Muon::BS2D) );
		muonseDB_BS3D->push_back( m->edB(pat::Muon::BS3D) );
		
		muonsBT_Dz->push_back( m->muonsBestTrack()->dz(vtxPos));
		muonsBT_Dxy->push_back( m->muonsBestTrack()->dxy(vtxPos));
		muonsBTe_Dz->push_back(m->muonBestTrack()->dzError());//error not calculated wrt to beamspot?????
		muonsBTe_Dxy->push_back(m->muonBestTrack()->dxyError());
		//Best track is calculated WRT vtxpos i.e. fitted beamspot/PV

		
        
  //    }
      //add muons to clean from jets 
      //if(isMediumID && miniIso < maxMuMiniIso_ && m->pt() > minMuPtForMuon2Clean_) mu2Clean->push_back(*m);

    }// end loop over pat muons
  }//end if vert size > 0
	//not sure if this is necessary (doMuonVeto is false in config)
	 // determine result before losing ownership of the pointer
  	//bool result = (doMuonVeto_ ? (prod->size() == 0) : true);
	
  *nMuons = prod->size();
	// store in the event
  iEvent.put(std::move(prod));
  iEvent.put(std::move(muonsIDtype), "muonsIDtype");
  iEvent.put(std::move(muonsFlagLoose), "muonsFlagLoose");
  iEvent.put(std::move(muonsFlagMedium), "muonsFlagMedium");
  iEvent.put(std::move(muonsFlagMediumPlus), "muonsFlagMediumPlus");
  iEvent.put(std::move(muonsFlagTight), "muonsFlagTight");
  iEvent.put(std::move(muonsLVec), "muonsLVec");
  iEvent.put(std::move(muonsCharge), "muonsCharge");
  iEvent.put(std::move(muonsMtw), "muonsMtw");
  iEvent.put(std::move(muonsRelIso), "muonsRelIso");
  iEvent.put(std::move(muonsMiniIso), "muonsMiniIso");
  iEvent.put(std::move(muonspfActivity), "muonspfActivity");
  iEvent.put(std::move(nMuons), "nMuons");
  iEvent.put(std::move(muonsDB_PV2D), "muonsDB_PV2D");
  iEvent.put(std::move(muonsDB_PV3D), "muonsDB_PV3D");
  iEvent.put(std::move(muonsDB_PVDZ), "muonsDB_PVDZ");
  iEvent.put(std::move(muonsDB_BS2D), "muonsDB_BS2D");
  iEvent.put(std::move(muonsDB_BS3D), "muonsDB_BS3D");
  iEvent.put(std::move(muonseDB_PV2D), "muonseDB_PV2D");
  iEvent.put(std::move(muonseDB_PV3D), "muonseDB_PV3D");
  iEvent.put(std::move(muonseDB_PVDZ), "muonseDB_PVDZ");
  iEvent.put(std::move(muonseDB_BS2D), "muonseDB_BS2D");
  iEvent.put(std::move(muonseDB_BS3D), "muonseDB_BS3D");
  iEvent.put(std::move(muonsBT_Dz), "muonsBT_Dz");
  iEvent.put(std::move(muonsBT_Dxy), "muonsBT_Dxy");
  iEvent.put(std::move(muonsBTe_Dz), "muonsBTe_Dz");
  iEvent.put(std::move(muonsBTe_Dxy),"muonsBTe_Dxy");


}//end filter
//kinematic variable fn's
double prodMuons::getMtw(const pat::Muon & muon, const reco::MET & met)
{
	return sqrt( 2*( (met)[0].pt()*m.pt() -( (met)[0].px()*m.px() + (met)[0].py()*m.py() ) ) );

}
///Isolation functions
double prodMuons::getRelIso(const pat::Muon & muon)
{
	return (m.pfIsolationR04().sumChargedHadronPt + std::max(0., m.pfIsolationR04().sumNeutralHadronEt + m.pfIsolationR04().sumPhotonEt - 0.5*m.pfIsolationR04().sumPUPt) ) / m.pt();
}
///POG classifications functions
bool prodMuons::isLooseMuon(const pat::Muon & muon)
{
  bool isLoose = true;
  if(doMuonID_)
  {
    isLoose = muon.isLooseMuon();
  }
  return isLoose;
}
bool prodMuons::isMediumMuon(const pat::Muon & muon)
{
	 bool isMedium = true;
  //medium WP + dz/dxy cuts
  bool goodGlob = muon.isGlobalMuon() && 
                  muon.globalTrack()->normalizedChi2() < 3 && 
                  muon.combinedQuality().chi2LocalPosition < 12 && 
                  muon.combinedQuality().trkKink < 20; 

		
    isMedium = muon.isLooseMuon() && 
               muon.innerTrack()->validFraction() > 0.8 && //revert back to standard before we introduce the HIPs and mitigation, see email in : https://hypernews.cern.ch/HyperNews/CMS/get/susy/2250.html
               //muon.innerTrack()->validFraction() > 0.49 && // Short term ID tunning to mitigate the HIP impact on lower tracker hit efficiencies
               muon.segmentCompatibility() > (goodGlob ? 0.303 : 0.451);
	
	return isMedium;
}
bool prodMuons::isMediumPlusMuon(const pat::Muon & muon, const reco::Vertex::Point & vtxPos )
{
  //Always default to true. If don't do muon ID (i.e., doMuonID_ = false) or don't do muon vtx (i.e., doMuonVtx_ = false), then the muon passes!
  bool isMedium = true, isMediumVtx = true;
  //medium WP + dz/dxy cuts
  bool goodGlob = muon.isGlobalMuon() && 
                  muon.globalTrack()->normalizedChi2() < 3 && 
                  muon.combinedQuality().chi2LocalPosition < 12 && 
                  muon.combinedQuality().trkKink < 20; 

    isMedium = muon.isLooseMuon() && 
               muon.innerTrack()->validFraction() > 0.8 && //revert back to standard before we introduce the HIPs and mitigation, see email in : https://hypernews.cern.ch/HyperNews/CMS/get/susy/2250.html
               //muon.innerTrack()->validFraction() > 0.49 && // Short term ID tunning to mitigate the HIP impact on lower tracker hit efficiencies
               muon.segmentCompatibility() > (goodGlob ? 0.303 : 0.451);



    isMediumVtx = muon.dB() < maxMuD0_ && fabs(muon.muonBestTrack()->dz(vtxPos)) < maxMuDz_;


  bool isMediumPlus = isMedium && isMediumVtx;
  return isMediumPlus; 
}
bool prodMuons::isTightMuon(const pat::Muon & muon, const reco::Vertex::Point & vtxpos)
{
  bool isTight = true;
  if(debug_ && (muon.muonID("AllGlobalMuons") != 0) ) 
  {
    std::cout << " (pt,eta,phi) "<<muon.pt()<<", "<<muon.eta()<<", "<<muon.phi() << " "
              << " isPFMuon " << muon.isPFMuon() << "\n"
              << " NormChi2 "<<muon.globalTrack()->normalizedChi2() 
              << ", NValidMuonHits " << muon.globalTrack()->hitPattern().numberOfValidMuonHits()
              << ", MMatchedStations " << muon.numberOfMatchedStations()
              << ", NValidPixelHits " << muon.innerTrack()->hitPattern().numberOfValidPixelHits()
              << ", trackerLayersWithMeasurement "  << muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() 
              << "\n"
              << ", dxy " << std::abs(muon.innerTrack()->dxy(vtxpos))
              << ", dz " << std::abs(muon.vz() - vtxpos.z() )
              << std::endl;

    std::cout << " sumChargedHadronPt " << muon.pfIsolationR04().sumChargedHadronPt
              << ", sumNeutralHadronEt " << muon.pfIsolationR04().sumNeutralHadronEt
              << ", sumPhotonEt " << muon.pfIsolationR04().sumPhotonEt
              << ", sumPUPt " << muon.pfIsolationR04().sumPUPt 
              << ", relIso " <<  (muon.pfIsolationR04().sumChargedHadronPt + std::max(0., muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5*muon.pfIsolationR04().sumPUPt) )/ muon.pt()
              << std::endl;
  }

  // ID cuts - always ask isGlobalMuon()
  if (muon.muonID("AllGlobalMuons") == 0){ isTight = false; return isTight; }
  if (doMuonID_) 
  {
    if(!muon.isPFMuon() ) isTight = false; 
    if( muon.globalTrack()->normalizedChi2() >= 10. ) isTight = false;
    if( muon.globalTrack()->hitPattern().numberOfValidMuonHits() <=0 ) isTight = false;
    if( muon.numberOfMatchedStations() <=1 ) isTight = false;
    if( muon.innerTrack()->hitPattern().numberOfValidPixelHits() == 0) isTight = false;
    if( muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() <=5 ) isTight = false;
    if(debug_) {std::cout << "PassedMuon ID" << std::endl;}
  }

  // vertex association cuts - ignore if no vertex (see further)
  if (doMuonVtx_) 
  {
    if (std::abs(muon.innerTrack()->dxy(vtxpos)) >= maxMuD0_) isTight = false;
    if (std::abs(muon.innerTrack()->dz(vtxpos))  >= maxMuDz_) isTight = false;
    if(debug_) {std::cout << "PassedMuon Vtx Association" << std::endl;}
  }
  return isTight; 
}
