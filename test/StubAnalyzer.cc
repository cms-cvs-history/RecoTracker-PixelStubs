// File: StubAnalyzer.cc
// Description:  see StubAnalyzer.h
// Author: Aaron Dominguez, Jason Keller, Tony Kelly (UNL)
// Creation Date:  OGU Jun. 23 2006 Revised version.
//--------------------------------------------
#include <memory>
#include <string>
#include <iostream>
#include "RecoTracker/PixelStubs/test/StubAnalyzer.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

/*
StubAnalyzer is a class, with data members conf, stub tree, 
and stubfile.
*/
StubAnalyzer::StubAnalyzer(edm::ParameterSet const& conf) : 
  conf_(conf), stubtree_(0), stubfile_(0)
  {
  }
StubAnalyzer::~StubAnalyzer()
  { 
  } 

void StubAnalyzer::endJob() 
{
  std::cout << " StubAnalyzer::endJob" << std::endl;
  stubfile_->Write();
  stubfile_->Close();
}//endJob
 
void StubAnalyzer::beginJob(const edm::EventSetup& es)
{
   //Make a New TTree
    std::string rootFile = conf_.getParameter<std::string>("RootFile");
    stubfile_ = new TFile(rootFile.c_str(),"RECREATE");
    stubtree_ = new TTree("stubtree_","Here it goes.");
 
    int buffer = 64000;

    //Branches
    std::cout << "Making Tracknum/Seednum Branch:" << std::endl;
    stubtree_->Branch("tracknseed", &tracknseed_.tracknum, "tracknum/I:seednum");

    std::cout << "Making Local Seed Position Branches:" << std::endl;
    stubtree_->Branch("numlocal",&local_.localelem, "localelem/I", buffer);
    stubtree_->Branch("localx",local_.x1,"x1[localelem]/D",buffer);
    stubtree_->Branch("localy",local_.y1,"y1[localelem]/D",buffer);
    stubtree_->Branch("localz",local_.z1,"z1[localelem]/D",buffer);

    std::cout << "Making Global Seed Position Branches:" << std::endl;
    stubtree_->Branch("numglobal",&global_.globalelem,"globalelem/I", buffer);
    stubtree_->Branch("globalx",global_.gx1,"gx1[globalelem]/D",buffer);
    stubtree_->Branch("globaly",global_.gy1,"gy1[globalelem]/D",buffer);
    stubtree_->Branch("globalz",global_.gz1,"gz1[globalelem]/D",buffer);

    std::cout << "Making sim Branches:" << std::endl;
    stubtree_->Branch("numsim",&sim_.simelem,"simelem/I", buffer);
    stubtree_->Branch("subdetId",sim_.subdetid,"subdetid[simelem]/I",buffer);
    stubtree_->Branch("TID", sim_.TID, "TID[simelem]/i", buffer);
    stubtree_->Branch("PID", sim_.PID, "PID[simelem]/I", buffer);

    std::cout << "Making simtrack Branches:" << std::endl;
    stubtree_->Branch("numsimtrack", &simtrack_.simtrackelem, "simtrackelem/I", buffer);
    stubtree_->Branch("simeta", simtrack_.eta, "eta[simtrackelem]/D", buffer);
    stubtree_->Branch("simphi", simtrack_.phi, "phi[simtrackelem]/D", buffer);
    stubtree_->Branch("simtheta", simtrack_.theta, "theta[simtrackelem]/D", buffer);
    stubtree_->Branch("simpT", simtrack_.pT, "pT[simtrackelem]/D", buffer);
    stubtree_->Branch("simPID", simtrack_.PID, "PID[simtrackelem]/I", buffer);

    std::cout << "Making recotrack Branches:" << std::endl;
    stubtree_->Branch("numrecotrack", &recotrack_.recotrackelem, "recotrackelem/I", buffer);
    stubtree_->Branch("recoeta", recotrack_.eta, "eta[recotrackelem]/D", buffer);
    stubtree_->Branch("recophi", recotrack_.phi, "phi[recotrackelem]/D", buffer);
    stubtree_->Branch("recotheta", recotrack_.theta, "theta[recotrackelem]/D", buffer);
    stubtree_->Branch("recocharge", recotrack_.charge, "charge[recotrackelem}/I", buffer);
    stubtree_->Branch("recopT", recotrack_.pT, "pT[recotrackelem]/D", buffer);
    stubtree_->Branch("recod0", recotrack_.d0, "d0[recotrackelem]/D", buffer);
    stubtree_->Branch("recodz", recotrack_.dz, "dz[recotrackelem]/D", buffer);
    stubtree_->Branch("recocurveT", recotrack_.curveT, "curveT[recotrackelem]/D", buffer);

    std::cout << "Making recoerror Branches:" << std::endl;
    stubtree_->Branch("numrecoerror", &recoerror_.recoerrorelem, "recoerrorelem/I", buffer);
    stubtree_->Branch("recothetaerror", recoerror_.thetaerror, "thetaerror[recoerrorelem]/D", buffer);
    stubtree_->Branch("recod0error", recoerror_.d0error, "d0error[recoerrorelem]/D", buffer);
    stubtree_->Branch("recodzerror", recoerror_.dzerror, "dzerror[recoerrorelem]/D", buffer);
    stubtree_->Branch("recocurveTerror", recoerror_.curveTerror, "curveTerror[recoerrorelem]/D", buffer);

    std::cout << "Done Making Branches!" << std::endl;
}//beginJob

void StubAnalyzer::analyze(const edm::Event& e, const edm::EventSetup& es)
  {
    using namespace edm;
    using namespace std;

    init();

    std::string trackFitter = conf_.getParameter<std::string>("TrackFitter");
    edm::Handle<reco::TrackCollection> recotracks; 
    e.getByLabel(trackFitter, recotracks);

    std::string seedGenerator = conf_.getParameter<std::string>("SeedGenerator");         
    edm::Handle<TrajectorySeedCollection> seeds;
    e.getByLabel(seedGenerator, seeds);

    //  Display on Monitor Number of seeds and tracks.
    LogInfo("Stubs") << "Seeds: " << seeds->size() << " " << "Tracks: " << recotracks->size();
    //std::cout << "Setting tracknum" << std::endl;
    tracknseed_.tracknum = recotracks->size();
    //std::cout << "Setting seednum" << std::endl;
    tracknseed_.seednum  = seeds->size();
 
    //std::cout << "Retrieving Handles" << std::endl;
    edm::ESHandle<TrackerGeometry> geom;
    es.get<TrackerDigiGeometryRecord>().get( geom );
    const TrackerGeometry& Tracker(*geom);

    //--- Get the simtracks for matching
    Handle<edm::SimTrackContainer> simtracks;
    e.getByLabel("g4SimHits",simtracks);

    TrackerHitAssociator associate(e);
    std::vector<PSimHit> pixhit;

    //std::cout << "Entering for loop to fill TrajectorySeed \"s\" class" << std::endl;

  TrajectorySeed s;
  DetId detId1;
  for(int i = 0; i != tracknseed_.seednum; ++i)
  { 
     s = (*seeds)[i];
     pixhit.clear();
     pixhit = associate.associateHit(*s.recHits().first);

     if(!pixhit.empty())
     {
       //std::cout << "Running fillLocal" << std::endl;
       fillLocal(s, i);

       detId1 = ((*s.recHits().first).geographicalId());
       const PixelGeomDetUnit *GeomDetUnit1 = dynamic_cast<const PixelGeomDetUnit*> (Tracker.idToDet(detId1));
       //std::cout << "Running fillGlobal" << std::endl;
       fillGlobal(GeomDetUnit1, i);

       //std::cout << "Running fillSim" << std::endl;
       for(std::vector<PSimHit>::const_iterator pixiter = pixhit.begin(); pixiter != pixhit.end(); ++pixiter)
       {
         fillSim(pixiter, detId1, i);
       }//for pixiter
     }//if pixhit
  }//for seednum

  SimTrack st;
  //std::cout << "Running fillSimTrack" << std::endl;
  for(int i = 0; i != (int)simtracks->size(); ++i)
  {
    st = (*simtracks)[i];
    fillSimTrack(st, i);
  }//for simtrack

  reco::Track rt;
  //std::cout << "Running fillRecoTrack" << std::endl;
  for(int i = 0; i != tracknseed_.tracknum; ++i)
  {
    rt = (*recotracks)[i];
    fillRecoTrack(rt, i);
    fillRecoError(rt, i);
  }//for recotrack
    
  //std::cout << "Filling tree" << std::endl;
  stubtree_->Fill();
}//analyze

void StubAnalyzer::fillSim(std::vector<PSimHit>::const_iterator ipix, DetId detid, int iter)
{
   sim_.simelem = iter;
   sim_.subdetid[iter] = detid.subdetId();
   sim_.TID[iter] = (*ipix).trackId();
   sim_.PID[iter] = (*ipix).particleType();
}

void StubAnalyzer::fillLocal(const TrajectorySeed& s, int iter)
{
  local_.localelem = iter;
  local_.x1[iter] = (*s.recHits().first).localPosition().x(); 
  local_.y1[iter] = (*s.recHits().first).localPosition().y(); 
  local_.z1[iter] = (*s.recHits().first).localPosition().z();
}

void StubAnalyzer::fillGlobal(const PixelGeomDetUnit *pixgeom, int iter)
{
  GlobalPoint GP1 = 
     pixgeom->surface().toGlobal(Local3DPoint(local_.x1[iter],local_.y1[iter],local_.z1[iter]));
  global_.globalelem= iter;
  global_.gx1[iter] = GP1.x();
  global_.gy1[iter] = GP1.y();
  global_.gz1[iter] = GP1.z();
}

void StubAnalyzer::fillSimTrack(const SimTrack& trks, int iter)
{
  simtrack_.simtrackelem = iter;
  simtrack_.eta[iter] = trks.momentum().eta();
  simtrack_.phi[iter] = trks.momentum().phi();
  simtrack_.theta[iter] = trks.momentum().theta();
  simtrack_.pT[iter] = trks.momentum().perp();
  simtrack_.PID[iter] = trks.type();
}

void StubAnalyzer::fillRecoTrack(const reco::Track& t, int iter)
{
  recotrack_.recotrackelem = iter;
  recotrack_.eta[iter] = t.eta();
  recotrack_.phi[iter] = t.phi();
  recotrack_.theta[iter] = t.theta();
  recotrack_.charge[iter] = t.charge();
  recotrack_.pT[iter] = t.pt() * t.charge();
  recotrack_.d0[iter] = t.d0();
  recotrack_.dz[iter] = t.dz();
  recotrack_.curveT[iter] = t.qoverp();
}

void StubAnalyzer::fillRecoError(const reco::Track& t, int iter)
{
  recoerror_.recoerrorelem = iter;
  recoerror_.thetaerror[iter] = t.thetaError();
  recoerror_.d0error[iter] = t.d0Error();
  recoerror_.dzerror[iter] = t.dzError();
  recoerror_.curveTerror[iter] = t.qoverpError();
}

void StubAnalyzer::init()
{
  tracknseed_.tracknum = 0;
  tracknseed_.seednum = 0;
  local_.localelem = 0;
  global_.globalelem = 0;
  sim_.simelem = 0;
  simtrack_.simtrackelem = 0;
  recotrack_.recotrackelem = 0;
  recoerror_.recoerrorelem = 0;
}

