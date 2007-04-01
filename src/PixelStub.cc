#include "RecoTracker/PixelStubs/interface/PixelStub.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include <cmath>
#include <algorithm>
#include "TMath.h"

PixelStub::PixelStub(const SiPixelRecHit &hit, double factor)
  : wyOffset_(13,0), betaRes_(13,0), hit_(hit), betaCutFactor(factor)
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
}
PixelStub::~PixelStub(){}

const SiPixelRecHit& PixelStub::hit() const {
  return hit_;
}

double PixelStub::wy(bool fancy) const {
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
  return 90.0 - piOver2MinusBeta();
}

bool PixelStub::compatible(const PixelStub &otherstub) const {
  double beta_diff = std::abs(piOver2MinusBeta() - otherstub.piOver2MinusBeta());
  std::cout << hit_.cluster()->sizeY() << " " << beta_diff << std::endl;
  std::cout << "The cut factor is set at " << betaCutFactor << std::endl;
  if (beta_diff < betaCutFactor*std::sqrt( res()*res() + otherstub.res()*otherstub.res()) ) {
    std::cout << "Compatible!" << std::endl;
    return true;
  }
  else {
    return false;
  }
}
  
double PixelStub::res() const {
  const SiPixelCluster &c = (*hit_.cluster());
  int size_y = std::min(c.sizeY(),12);
  return betaRes_[size_y];
}
