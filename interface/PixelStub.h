#ifndef PixelStubs_PixelStub_h
#define PixelStubs_PixelStub_h
/** \class PixelStub PixelStub.h RecoTracker/PixelStubs/PixelStub.h 
 * Use shape in local y of pixel cluster to estimate the track incidence angle.
 * Can make comparisons for compatibilities.
 *  $Date: 2007/03/01 23:50:38 $
 *  $Revision: 1.1 $
 *  \author Aaron Dominguez (UNL)
 */
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include <vector>
class PixelStub {
 public:
  PixelStub(const SiPixelRecHit &hit, double factor);
  ~PixelStub();

  const SiPixelRecHit &hit() const;
  double wy(bool fancy=false) const;
  double piOver2MinusBeta() const;
  double beta() const;
  bool compatible(const PixelStub &otherstub) const;

 private:
  double res() const;
  double betaCutFactor;

  std::vector<double> wyOffset_;
  std::vector<double> betaRes_;
  const SiPixelRecHit &hit_;
};
#endif
