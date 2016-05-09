#ifndef TGEM_HH
#define TGEM_HH

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "TRandom.h"
#include "TMath.h"
#include <vector>
#include "TDimension.h"
#include <set>
class TGem{

 public:
  
  TGem(TRandom *RNG,TDimension *Dimension);
  std::pair<double,double>                   Sampling(double x, double y);
  std::pair<double,double>                   GemSampling(double xin, double yin);
  inline void SetGain(double gain)  {Gain  = gain;}
  inline void SetPitch(double pitch){Pitch = pitch;}
  inline double GetPitch()  {return Pitch;}
  int GetSecondaryElectrons();
  std::vector<std::pair<double,double> >  DiffusionofSecondaryElectrons(std::vector<std::pair<double,double> > VectorIn ) ;
  inline int GetPac() {return Pac;}
  inline double GetHole() {return Hole;}
 private:  
  double Gain;
  double Pitch;
  double AvalangeSigma;
  TRandom *rng;
  int Pac;
  double Hole;
  double XStep,YStep;
  double TraslX,TraslY;
};

#endif
