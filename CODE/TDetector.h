#ifndef TDETECTOR_HH
#define TDETECTOR_HH

#include "TRandom.h"
#include "TGasMixture.h"
#include "TGem.h"
#include "TPixMap.h"
#include "TTrack.h"
#include <vector>
#include <algorithm>
//#include <pair.h>
#include "TReadout.h"
#include "TDimension.h"
#include <list>
#include "TView.h"
#include "TGaxis.h"
struct ADC
{
  int Channel;
  int Charge;
  int DigiCharge;
};

class TDetector
{
 public:
  TDetector(TGasMixture *MIXTURE, TRandom *RND, Bool_t VERBOSE,TDimension *Dimension);
  void TransferToGEM(TTrack Track);
  //  void GEMSampling();
  //  void TranferToRadOutPlane();
  
  //  void DrawSignal(std::vector<ADC> Digi);
  
  std::vector<ADC>        mySampling(std::vector<std::pair<double,double> >  XYSec,TReadout Readout);
  int                     CountChannels(std::vector<int> chvect, int ch);
  //  std::vector<ADC>        ReadOutSampling(TTrack track);
  void                    myDrawSignal(std::vector<ADC> Digi,TReadout Readout);
  
 private:
  
  int  Digitization(int Charge);
  TRandom *rng;
  //  TPixMap *PixMap;
  TGem    *Gem;
  //TReadout *Readout;
  double Gem_Pitch;
  double Gem_Gain;
  double Gem_Radius;
  double Pitch;
  std::pair<double,double> Pi;
  double NoiseRMS;
  int    Pedestal;
  double ElectronsPerADC;
  double PixelX;
  double PixelY;
 
};
#endif
