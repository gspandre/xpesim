#ifndef TEXPERIMENT_HH
#define TEXPERIMENT_HH
#include <iostream>
#include <fstream>
#include <pair.h>
#include "TRandom3.h"
#include "TGasMixture.h"
#include "TDetector.h"
#include "TTree.h"
#include "TFile.h"
#include <vector>
#include "TReadout.h"
#include "TH2D.h"
#include "TDimension.h"
#include "TGem.h"
#include "TView.h"
#include "THexagon3d.h"
#include "TCluster.h"
#include "TPaveText.h"
#include "TText.h"
//#include "TMDP.h"
//#include "TVarPolFlux.h"
#include "TSystem.h"
#include "TTreeAnalysis.h"
class TExperiment
{
 public:
  TExperiment(TRandom3 *RND,TString suffix="");
  ~TExperiment(){;}
  void Delete();
  void     SetMixID(Int_t MIXID);
  void     SetSource(Double_t POLARIZATIONANGLE,Double_t POLARIZATIONDEGREE,Double_t  PHOTONENERGY);
  void     SetSource(Double_t POLARIZATIONANGLE,Double_t POLARIZATIONDEGREE,TH1D *obsSpectrum);
  void     SetThickness(double THICKNESS);
  void     SetPressure(double PRESSURE);
  void     Generalsetup(double GP = 50., double PitchSet = 50.);
  TString  GetName();
  void     EventsTree( int Number=1000, TString File= "Eventstree" );
  void     AnalyzeTree();
  void     ReadResult(double &mucut,double &Eff,double &Effcut,double &mu,double &FitProbability, double  &FitProbabilityCut,double &mufirst);
  double   GetEfficiencyMixture(double myene=0);
  TGraph  *GetEfficiencyGraph();
  
  TH1D     *GetClusterH(double PEnergy);
  TH1D     *GetRecordDimH(double PEnergy);
  //
  inline void  SetDGem_Pitch(double GP){ Gem_Pitch =  GP;} 
  inline void  SetDPitch(double P){ Pitch =  P;}
  inline void  SetMaxConvertedEvents(long int Nmax){MaxGeneratedEvents =  Nmax;}
  //
  
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
  TRandom3* rnd;
  TDimension* Dimension;
  TGasMixture* Mixture;
  TSource* Source;
  TH1D *ObsSpectrum;
  TString m_suffix;
  
};
#endif
