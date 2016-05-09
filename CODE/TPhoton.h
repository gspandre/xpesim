#ifndef TPHOTON_HH
#define TPHOTON_HH

#include "TSource.h"
#include "TXYZ.h"
#include <TF1.h>
#define ELECTRON_MASS 511.0

/// \brief Class describing a photon.
/// A photon is basically identified by its energy and conversion point.
class TPhoton{
 public:
  TPhoton(TSource *SOURCE, TXYZ CONVERSION_POINT, TRandom *RND, Bool_t VERBOSE=1);
  ~TPhoton();

  inline Double_t GetEnergy()           {return Energy;}

  inline TXYZ GetConversionPoint() {return ConversionPoint;}

  Double_t        GetPhotoelectronTheta();
  Double_t        GetPhotoelectronPhi();
  Double_t        GetAugerElectronTheta();
  Double_t        GetAugerElectronPhi();

 private:
  Bool_t   Verbose;
  TSource* Source;
  Double_t Energy;

  TXYZ     ConversionPoint;

  Double_t PolarizationAngle;
  Double_t PolarizationDegree;
  TF1*     ThetaDistribution;
  TF1*     PhiDistribution;
  TRandom* rnd;
};

#endif
