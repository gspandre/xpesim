#ifndef TVARPOLFLUX_HH
#define TVARPOLFLUX_HH
#include "MC.h"

class TVarPolFlux

{ 
 public:
  TVarPolFlux(double G,double I);
  ~TVarPolFlux(){;}
  void                       DrawPol();
  double                     DegreetoRad(double AngD);
 
 private:
 
  TF1   *S1;
  TF1   *S2;
  TGraph *S1t;
  TGraph *S2t;
  TF1   *Angle;
  TGraph *S;
  double Gamma,Ir,SS;
  
};
#endif
  
