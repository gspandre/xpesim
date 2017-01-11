#ifndef TEXPERIMENT_HH
#define TEXPERIMENT_HH

#include "MC.h"
#include "TGasMixture.h"
#include "TDetector.h"
#include "TReadout.h"
#include "TDimension.h"
#include "TGem.h"
#include "THexagon3d.h"
#include "TCluster.h"
#include "TTreeAnalysis.h"

class TExperiment
{
 public:
  TExperiment(TRandom *RND,TString suffix="");
  ~TExperiment(){;}
  void     Delete();
  void     SetMixID(Int_t MIXID);
  void     SetSource(Double_t POLARIZATIONANGLE,Double_t POLARIZATIONDEGREE,Double_t  PHOTONENERGY);
  void     SetSource(Double_t POLARIZATIONANGLE,Double_t POLARIZATIONDEGREE,TH1D *obsSpectrum);
  void     SetThickness(double THICKNESS);
  void     SetPressure(double PRESSURE);
  void     Generalsetup(double GP = 50., double PitchSet = 50.);
  TString  GetName();
  TString  GetNameforFile();
  TString  GetMixType();
  void     EventsTree(int Number=1000, TString File= "Eventstree");
  void     G4MCEventsLoop(int FirstEvt, int Number, TString File);
  void     AnalyzeTree();
  void     ReadResult(double &mucut,double &Eff,double &Effcut,double &mu,double &FitProbability, double  &FitProbabilityCut,double &mufirst);
  double   GetEfficiencyMixture(double myene=0);
  TGraph  *GetEfficiencyGraph();

  Double_t GetElectronRange(double Energy);
  Double_t GetPhotoelectricCrossSection(Double_t Energy);
  Double_t GetAbsorptionLenght(Double_t  Energy);
  std::vector<TGraph*> GetAbsProbGraph(Double_t Energy);

  TH1D     *GetClusterH(double PEnergy);
  TH1D     *GetRecordDimH(double PEnergy);
  TH1D     *pheRange, *pheTrueRange;
  TH1D     *augRange, *augTrueRange;
  //
  inline void  SetDGem_Pitch(double GP){ Gem_Pitch =  GP;} 
  inline void  SetDPitch(double P){ Pitch =  P;}
  inline void  SetMaxConvertedEvents(long int Nmax){MaxGeneratedEvents =  Nmax;}
  //
  int nauger, nphel;

 private:
  bool SetAstroSource;
  double Gem_Pitch ;
  double   Pitch ;
  Int_t MixID;
  long int MaxGeneratedEvents;
  Double_t PolarizationAngle;
  Double_t PolarizationDegree;
  double PhotonEnergy;
  double Thickness;
  double Pressure;
  TRandom* rnd;
  TDimension* Dimension;
  TGasMixture* Mixture;
  TSource* Source;
  TH1D *ObsSpectrum;
  TString m_suffix;
  std::vector<TGraph*> AbsorptionProbabilityGr;
  //std::vector<TXYZ> G4PrimaryIonizationV;
};
#endif
