#ifndef PixelStubs_PixelStub_h
#define PixelStubs_PixelStub_h
/** \class PixelStub PixelStub.h RecoTracker/PixelStubs/PixelStub.h 
 * Use shape in local y of pixel cluster to estimate the track incidence angle.
 * Can make comparisons for compatibilities.
 *  $Date: 2007/05/12 13:53:13 $
 *  $Revision: 2.0 $
 *  \author Aaron Dominguez (UNL) edited by Dave Fehling
 */

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
//#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include <vector>
#include "TMath.h"

class MagneticField;
class PixelCPEBase;

class PixelStub {
 public:
  PixelStub(const SiPixelRecHit &hit, double factor, double tempfactor, int method, const edm::ParameterSet &cfg, const edm::EventSetup &es, const PixelCPEBase &cpe);
  ~PixelStub();

  const SiPixelRecHit &hit() const;
  double wy(bool fancy=false) const;
  double piOver2MinusBeta() const;
  double beta() const;
  bool compatible(const PixelStub &otherstub) const;

	double crudeAngles() const;
	void templateProbs(const GlobalVector &gv) const;
	GlobalPoint gpMaker() const;

private:
  double res() const;

  std::vector<double> wyOffset_;
  std::vector<double> betaRes_;
  const SiPixelRecHit &hit_;
  double betaCutFactor_;
	double tempCutFactor_;
	mutable double alpha_;
	mutable double beta_;
	int method_;
	const edm::ParameterSet cfg_;
	const edm::EventSetup &es_;
	const GeomDetUnit *theDet_;
	const PixelGeomDetUnit *pixDet_;
	const MagneticField *magField_;
 
	mutable const PixelCPEBase * cpe_;
	mutable float probX_;
	mutable float probY_;
	mutable int qBin_;
	
	double RadtoDeg; //we should make this static const;
};
#endif
