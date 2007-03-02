#ifndef StubAnalyzer_h
#define StubAnalyzer_h

/** \class StubAnalyzer
 *
 * StubAnalyzer is the EDProducer subclass which finds seeds
 *
 * \authors Jason Keller, Tony Kelly (UNL)
 *
 * \version   First Version 6/23/2006

 *
 ************************************************************/

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Common/interface/EDProduct.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"

#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Track/interface/EmbdSimTrack.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "SimDataFormats/Vertex/interface/EmbdSimVertex.h"
#include "SimDataFormats/Track/interface/EmbdSimTrackContainer.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TTree.h"
#include "TFile.h"

class StubAnalyzer : public edm::EDAnalyzer
{
 public:
  
  explicit StubAnalyzer(const edm::ParameterSet& conf);
  
  virtual ~StubAnalyzer();
  
  virtual void beginJob(const edm::EventSetup& es);

  virtual void analyze(const edm::Event& e, const edm::EventSetup& es);

  virtual void endJob();

 protected:
  void fillLocal(const TrajectorySeed&, int);
  void fillGlobal(const PixelGeomDetUnit*, int);
  void fillSim(std::vector<PSimHit>::const_iterator, DetId, int);
  void fillSimTrack(const EmbdSimTrack&, int);
  void fillRecoTrack(const reco::Track&, int);
  void fillRecoError(const reco::Track&, int);
  void init();
  
 private:
  edm::ParameterSet conf_;
  TTree *stubtree_;
  TFile *stubfile_;

  static const int maxelem = 5000;
  static const int maxtrack = 500;
  struct tracknseed
  {
    int tracknum;
    int seednum;
  }tracknseed_;
 
  struct local
  {
    int localelem;
    double x1[maxelem];
    double y1[maxelem];
    double z1[maxelem];
  } local_;

  struct global
  {
    int globalelem;
    double gx1[maxelem];
    double gy1[maxelem];
    double gz1[maxelem];
  } global_;
  
 struct sim
 {
   int simelem;
   int subdetid[maxelem];
   unsigned int TID[maxelem];
   int PID[maxelem];
 } sim_;
 
  struct simtrack
  {
    int simtrackelem;
    double eta[maxtrack];
    double phi[maxtrack];
    double theta[maxtrack];
    double pT[maxtrack];
    int PID[maxtrack];
  } simtrack_;

  struct recotrack
  {
    int recotrackelem;
    double eta[maxtrack];
    double phi[maxtrack];
    double theta[maxtrack];
    int charge[maxtrack]; 
    double pT[maxtrack];
    double d0[maxtrack];
    double dz[maxtrack];
    double curveT[maxtrack];
    
  } recotrack_;

  struct recoerror
  {
    int recoerrorelem;
    double thetaerror[maxtrack];
    double d0error[maxtrack];
    double dzerror[maxtrack];
    double curveTerror[maxtrack];
  } recoerror_;
 
};


#endif
