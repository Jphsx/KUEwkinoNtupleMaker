#ifndef __myGenJetsPartSelector__
#define __myGenJetsPartSelector__

/* \class myGenJetsPartSelector
 * 
 * Selects particles for the GenJet input.
 * Deselect specified particles, also radiation from resoances.
 * Or only select partonic final state.
 * The algorithm is based on code of Christophe Saout. 
 *
 * \author: Andreas Oehler
 * 
 * 
 */

#include "DataFormats/Provenance/interface/Provenance.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/Common/interface/EDProductfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


class myGenJetsPartSelector : public edm::EDProducer {
  // collection type
 public:
  typedef std::vector<bool>     ParticleBitmap;
  typedef std::vector<const reco::GenParticle*> ParticleVector;
      
  myGenJetsPartSelector(const edm::ParameterSet & ); 
  ~myGenJetsPartSelector(); 
  // select object from a collection and 
  // possibly event content
  virtual void produce (edm::Event &evt, const edm::EventSetup &evtSetup);
      
  bool getPartonicFinalState() const { return partonicFinalState; }
  bool getExcludeResonances() const { return excludeResonances; }
  bool getTausAndJets() const { return tausAsJets; }
  double getPtMin() const { return ptMin; }
  const std::vector<unsigned int> &getIgnoredParticles() const
    { return ignoreParticleIDs; }
  
  void setPartonicFinalState(bool flag = true)
  { partonicFinalState = flag; }
  void setExcludeResonances(bool flag = true)
  { excludeResonances = flag; }
  void setTausAsJets(bool flag = true) { tausAsJets = flag; }
  void setPtMin(double ptMin) { this->ptMin = ptMin; }
  bool isParton(int pdgId) const;
  static bool isHadron(int pdgId);
  static bool isResonance(int pdgId);

  bool isIgnored(int pdgId) const;
  bool hasPartonChildren(ParticleBitmap &invalid,
			 const ParticleVector &p,
			 const reco::GenParticle *particle) const;
  
  enum ResonanceState {
    kNo = 0,
    kDirect,
    kIndirect
  };
  ResonanceState fromResonance(ParticleBitmap &invalid,
			       const ParticleVector &p,
			       const reco::GenParticle *particle) const;
  
  
  // iterators over selected objects: collection begin
  
 private:
  //container selected_;  //container required by selector
  myGenJetsPartSelector(){} //should not be used!
  
  edm::InputTag inTag;
  edm::EDGetTokenT<std::vector<reco::GenParticleCollection>> InTagTok_;
  int testPartonChildren(ParticleBitmap &invalid,
			 const ParticleVector &p,
			 const reco::GenParticle *particle) const;
  
  std::vector<unsigned int>	ignoreParticleIDs;
  std::vector<unsigned int> excludeFromResonancePids;
  void setExcludeFromResonancePids(const std::vector<unsigned int> &particleIDs);
  void setIgnoredParticles(const std::vector<unsigned int> &particleIDs);
  bool isExcludedFromResonance(int pdgId) const;
  
  bool			partonicFinalState;
  bool			excludeResonances;
  bool			tausAsJets;
  double			ptMin;

  bool find_mother(const reco::Candidate* genp, const int motherID);
  std::vector<int> ignoreParticleDecayIDs;
  bool invertignoreParticleDecay;

  int find_idx(const reco::Candidate* target);

  edm::Handle<reco::GenParticleCollection> genParticles;
  
};

#endif
