#ifndef TDETECTOR_HH
#define TDETECTOR_HH

#include "MC.h"
#include "TGasMixture.h"
#include "TGem.h"
#include "TPixMap.h"
#include "TTrack.h"
#include "TReadout.h"
#include "TDimension.h"

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
  std::vector<ADC>        mySampling(std::vector<std::pair<double,double> >  XYSec,TReadout Readout);
  int                     CountChannels(std::vector<int> chvect, int ch);
  void                    myDrawSignal(std::vector<ADC> Digi,TReadout Readout);
  
 private:
  
  int  Digitization(int Charge);
  TRandom *rng;
  TGem    *Gem;
  double Gem_Pitch;
  double Gem_Gain;
  double Gem_Radius;
  double Pitch;
  std::pair<int,int> Pi;
  double NoiseRMS;
  int    Pedestal;
  double ElectronsPerADC;
  double PixelX;
  double PixelY;
 
};
#endif
