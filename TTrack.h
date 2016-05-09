#ifndef TTRACK_HH
#define TTRACK_HH

#include "MC.h"
#include "TPhoton.h"
#include "TGasMixture.h"
#include "TXYZ.h"
#include "TDimension.h"


/// \brief Class describing a track.
/// A track is a track of a photoelectron (and, possibly, an Auger electron),
/// ejected by a photon, into a gas mixture. The basic interesting variables
/// are the 3D points in which the primary electrons are created.

class TTrack{
 public:
  TTrack(TPhoton *PHOTON, TGasMixture *MIXTURE, TRandom *RND, Bool_t VERBOSE=1);
  ~TTrack() {;}
  inline Int_t                 GetnPrimaryElectrons()  {return nPrimaryElectrons;}
  inline std::vector<TXYZ>     GetPrimaryIonizationV() {return PrimaryIonizationV;}
  inline std::vector<Int_t>    GetnPairspath() {return nPairspath;}
  inline std::vector<TXYZ>     GetPhotoelectronScatteringV() {return PhotoelectronScatteringV;}
  inline std::vector<TXYZ>     GetAugerScatteringV() {return AugerScatteringV;}
  void                         SetDimension(Double_t,Double_t,Double_t,Double_t);//Dim
  void                         PlotPath();
  void                         PlotPrimaryIonization();
  void                         Drift();
  void                         PropagatePhotoelectron(); 
 
  inline double                GetPhotoelectronPhi() {return PhotoelectronPhi; }
  inline double                GetPhotoelectronTheta() {return PhotoelectronTheta; }
  inline double                GetPhotoelectronPraticalRange() {return PhotoelectronPraticalRange; }
  inline double                GetPhotoelectronTrueRange() {return PhotoelectronTrueRange; }
  inline double                GetAugerPraticalRange() {return AugerPraticalRange; }
  inline double                GetAugerTrueRange() {return AugerTrueRange; }
  inline double                GetAugerEnergy() {return augerEnergy; }
  inline double                GetPhotoelEnergy(){return pheEnergy; }
 

  inline std::vector<std::pair<double,double> >    
    GetDiffElectronPosition() {return  DiffElectronPosition;}///diff
  inline int                   GetAugerCheck() {return  AugerCheck;}
  int nphe;
  int naug;

 private:
  
  Double_t              xdim;//dim
  Double_t              ydim;
  Double_t              zdim;
  Double_t              zmin;
  TXYZ                  ConversionPoint;
  TRandom*              rnd;
  Bool_t                Verbose;
  TPhoton*              Photon;
  TGasMixture*          Mixture;
  TElement*             ConvertingElement;
  Int_t                 nPrimaryElectrons;
  Int_t                 nPhotoelectronPrimaryElectrons;
  Int_t                 nPhotoelectronElasticScatterings;
  Int_t                 nAugerPrimaryElectrons;
  Int_t                 nAugerElasticScatterings;
  double                PhotoelectronPhi;
  double                PhotoelectronTheta;
  std::vector<TXYZ> PhotoelectronScatteringV;
  std::vector<TXYZ> AugerScatteringV;
  std::vector<std::pair<double,double> >  DiffElectronPosition;//diff
  std::vector<Double_t> PhotoelectronEnergy;
  std::vector<Double_t> AugerEnergy;
  std::vector<Int_t> nPairspath;
  std::vector<TXYZ> PrimaryIonizationV;
  Double_t              pheEnergy;
  Double_t              augerEnergy;
  
  Double_t              DistanceFromStart;
  Double_t              MaxDistanceFromStart;
  Double_t              PhotoelectronTrueRange;
  Double_t              PhotoelectronPraticalRange;
  Double_t              PhotoelectronStartToEndRange;
  Double_t              AugerTrueRange;
  Double_t              AugerPraticalRange;
  Double_t              AugerStartToEndRange;
  Double_t              ResidualEnergy;
  Double_t              CX;
  Double_t              CY;
  Double_t              CZ;
  void                  EvaluateNextStep(TString MODE="MOTT");
  Int_t                 GetnPairs(Double_t ENERGY_LOSS);
  int                    AugerCheck;/////
};

#endif
