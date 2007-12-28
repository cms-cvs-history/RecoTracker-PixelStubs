#include "RecoTracker/PixelStubs/interface/PixelStub.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerTopology/interface/RectangularPixelTopology.h"
#include "RecoLocalTracker/SiPixelRecHits/interface/PixelCPEBase.h"
#include "RecoLocalTracker/SiPixelRecHits/interface/PixelCPETemplateReco.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include <cmath>
#include <algorithm>
#include "TMath.h"

PixelStub::PixelStub(const SiPixelRecHit &hit, double factor, double tempfactor, int method, const edm::ParameterSet &cfg, const edm::EventSetup &es, const PixelCPEBase &cpe)
  : wyOffset_(13,0), betaRes_(13,0), hit_(hit), betaCutFactor_(factor), tempCutFactor_(tempfactor), alpha_(0), beta_(0), method_(method), cfg_(cfg), es_(es), cpe_(&cpe), RadtoDeg(180.0/TMath::Pi())

{
  wyOffset_[2]=-0.179757;
  wyOffset_[3]=-0.249349;
  wyOffset_[4]=-0.451085;
  wyOffset_[5]=-0.451085;
  wyOffset_[6]=-0.576694;
  wyOffset_[7]=-0.638618;
  wyOffset_[8]=-0.71524;
  wyOffset_[9]=-0.878925;
  wyOffset_[10]=-0.878925;
  wyOffset_[11]=-0.760377;
  wyOffset_[12]=-0.760377;

  betaRes_[1]=27;
  betaRes_[2]=12.3848;
  betaRes_[3]=5.10328;
  betaRes_[4]=3.86236;
  betaRes_[5]=2.18467;
  betaRes_[6]=1.62346;
  betaRes_[7]=1.04936;
  betaRes_[8]=0.83669;
  betaRes_[9]=0.828255;
  betaRes_[10]=0.709503;
  betaRes_[11]=0.5;
  betaRes_[12]=betaRes_[11];

	edm::ESHandle<TrackerGeometry> geom;
	es.get<TrackerDigiGeometryRecord>().get(geom);

	const GeomDetUnit *geoUnit = geom->idToDetUnit(hit_.geographicalId());
	theDet_ = geoUnit;
		
	edm::ESHandle<MagneticField> magFieldHandle;
	es.get<IdealMagneticFieldRecord>().get(magFieldHandle);
	magField_ = &(*magFieldHandle);

}
PixelStub::~PixelStub(){}

const SiPixelRecHit& PixelStub::hit() const {
  return hit_;
}

double PixelStub::wy(bool fancy) const { //method_ == 1
  const SiPixelCluster &c = (*hit_.cluster());
  int size_y = std::min(c.sizeY(),12);
  std::cout << "Size " << size_y << " Id " << hit_.geographicalId().rawId() << std::endl;
  double Wy = size_y - 1.0;
  if (fancy) {
    double qL = 0, qF = 0; // Charge on first & last col
    for (unsigned int i=0; i<c.pixels().size(); ++i) {
      std::cout << c.pixels()[i].y << " " << c.minPixelCol()+0.5 << " " << c.pixels()[i].adc << std::endl;
      if (std::abs(c.pixels()[i].y - (c.minPixelCol()+0.5))<0.1) qF += c.pixels()[i].adc;
      if (std::abs(c.pixels()[i].y - (c.maxPixelCol()+0.5))<0.1) qL += c.pixels()[i].adc;
    }
    double yhit = c.y();
    double yC = (size_y)/2.0 + c.minPixelCol();
    double Wi = (size_y-2.0);
    if ( std::abs(qL-qF) > 10 ) { // Use fancy method
      Wy = std::abs( ( yhit - yC) * 2.0 * (qL+qF) / (qL-qF) );
      if ( (Wy - Wi > 1.0) || (Wy - Wi) < 0.0 ) Wy = size_y - 1.0;
      else {
	Wy += Wi;
      }
    }
    Wy -= wyOffset_[size_y];
  }
  else 
    Wy -= wyOffset_[size_y];
  return Wy;  
}

double PixelStub::piOver2MinusBeta() const {
  double angle_calc= TMath::ATan2(wy()*0.0150,0.0285)*180/TMath::Pi();
  return angle_calc;
}

double PixelStub::beta() const {
	//PixelStubs - stubs
	if (method_==1) return 90.0 - piOver2MinusBeta();
	
	//Method similar to PixelStubs, using size and theThickness only - crude
	else if (method_==2) return crudeAngles();

	else return 0;
}

bool PixelStub::compatible(const PixelStub &otherstub) const {

	bool rtn = false;

	if (method_ < 3){
		double beta_diff = std::fabs(beta() - otherstub.beta());
		std::cout << "Cluster Size: " << hit_.cluster()->sizeY() << " Beta Diff:  " << beta_diff << std::endl;
		std::cout << "The cut factor is set at " << betaCutFactor_ << std::endl;
		if (beta_diff < betaCutFactor_*std::sqrt( res()*res() + otherstub.res()*otherstub.res())) {
			std::cout << "Compatible!" << std::endl;
			std::cout << "First stub beta: " << beta() << "\tSecond stub beta: " << otherstub.beta() << std::endl;
			rtn = true;
		}
	}
	else {
		GlobalPoint psGP = gpMaker();
		GlobalPoint opsGP = otherstub.gpMaker();
		GlobalVector gv( (psGP.x() - opsGP.x()), (psGP.y() - opsGP.y()), (psGP.z() - opsGP.z()));
		
		templateProbs(gv);
		otherstub.templateProbs(gv);
			
		/*		std::cout << "First stub alpha: " << alpha_ * RadtoDeg << "\tSecond stub alpha: " << otherstub.alpha_ * RadtoDeg << std::endl;
		std::cout << "First stub beta: " << beta_ * RadtoDeg << "\tSecond stub beta: " << otherstub.beta_ * RadtoDeg << std::endl;
		std::cout << "First stub prob in X: " << probX_ << "\t\tSecond stub prob in X: " << otherstub.probX_ << std::endl;
		std::cout << "First stub prob in Y: " << probY_ << "\t\tSecond stub prob in Y: " << otherstub.probY_ << std::endl;
		std::cout << "First stub qBin: " << qBin_ << "\t\tSecond stub qBin: " << otherstub.qBin_ << std::endl;*/

		if (method_ == 3 ) {
			if (fabs(fabs(beta_) - fabs(otherstub.beta_)) < betaCutFactor_ / RadtoDeg) {
				//std::cout << "Compatible!" << std::endl;
				rtn = true;
			}
		}
		else {
			if (probY_ > tempCutFactor_ && otherstub.probY_ > tempCutFactor_
				&& qBin_ < 4 && otherstub.qBin_ < 4 ) {//&& probX_ > tempCutFactor_ && otherstub.probX_ > tempCutFactor_) {
				//	std::cout << "Compatible!" << std::endl;
				rtn = true;
			}
		}
	}
	return rtn;
}
  
double PixelStub::res() const {
  const SiPixelCluster &c = (*hit_.cluster());
  int size_y = std::min(c.sizeY(),12);
  return betaRes_[size_y];
}

double PixelStub::crudeAngles() const { //method_ == 2
	const SiPixelCluster &clu = (*hit_.cluster());
	double theThickness = theDet_->surface().bounds().thickness();
	return TMath::ATan2(clu.sizeY()-1,theThickness)*180/TMath::Pi();
}

void PixelStub::templateProbs(const GlobalVector &gv) const { //method_ == 4
	const SiPixelCluster &clu = (*hit_.cluster());

	LocalVector lv = theDet_->surface().toLocal( gv );
	alpha_ = TMath::ATan2(lv.z(),lv.x());
	beta_ = TMath::ATan2(lv.z(),lv.y());

	std::pair<LocalPoint, LocalError> locValue = 
		cpe_->localParameters( clu, *theDet_, alpha_, beta_ );

	probX_ = cpe_->probabilityX();
	probY_ = cpe_->probabilityY();
	qBin_  = cpe_->qBin();
}

GlobalPoint PixelStub::gpMaker() const {
	//	const SiPixelCluster &clu = (*hit_.cluster());
	
	// get cluster center of gravity (of charge)
	//float xcenter = clu.x();
	//float ycenter = clu.y();
	
	// get the cluster position in local coordinates (cm) and covert to global coordinates (cm)
	//const PixelGeomDetUnit *pixDet_ = dynamic_cast<const PixelGeomDetUnit*>( theDet_ );
	
	//const RectangularPixelTopology *theTopol
	//	= dynamic_cast<const RectangularPixelTopology*>( & (pixDet_->specificTopology()) );
	
	//LocalPoint lp = theTopol->localPosition( MeasurementPoint(xcenter, ycenter) );
	LocalPoint lp = hit_.localPosition();
	
	GlobalPoint gp = theDet_->surface().toGlobal( lp );

	return gp;
}
