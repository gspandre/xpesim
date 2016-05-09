#ifndef TSOURCE_HH
#define TSOURCE_HH

#include <TF1.h>
#include <pair.h>
#include "TElement.h"

/// \brief Class describing a source.
class TSource{
 public:
  TSource(Double_t ENERGY, Double_t POLARIZATION_ANGLE, Double_t POLARIZATION_DEGREE,
	  TRandom *RND, Bool_t VERBOSE=1);
  ~TSource();
  
  inline TF1* GetThetaDistribution()     {return ThetaDistribution;}
  inline TF1* GetPhiDistribution()       {return PhiDistribution;}
  inline void SetEnergy(Double_t ENERGY) {Energy = ENERGY;}
  void        SetPolarizationAngle(Double_t POLARIZATION_ANGLE);
  void        SetPolarizationDegree(Double_t POLARIZATION_DEGREE); 
  Double_t    GetEnergy();
  void        SetPattern(TString Set);
  std::pair<double,double>   GetConversionXY();

 private:
  Bool_t   Verbose;
  Double_t Energy;
  Double_t PolarizationAngle;
  Double_t PolarizationDegree;
  TF1*     ThetaDistribution;
  TF1*     PhiDistribution;
  TRandom* rnd;
  TString PatternType;
};

#endif
