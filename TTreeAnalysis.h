#ifndef TTREEANALYSIS_HH
#define TTREEANALYSIS_HH

#include "MC.h"
#include "TGasMixture.h"
#include "TDetector.h"
#include "TReadout.h"
#include "TDimension.h"
#include "TGem.h"
#include "THexagon3d.h"
#include "TCluster.h"

class TTreeAnalysis
{
 public:
  TTreeAnalysis (TString FileName,TString DataTreeName);
  ~TTreeAnalysis()
    {
      delete rnd ;
      delete Dimension;
      delete Mixture;
      delete myDetector;
      // delete t;
    }

  std::vector<ADC>           ReadCluster(int NCluster);
  void                       ClusterAnalysis(int NCluster);
  void                       TreeAnalysis(Bool_t ParSelectFix = false );
  void                       SetClusterPar(int Clusterdim);
  void                       SetFixedClusterPar(double ain,double bin, double cin);
  inline int                 GetClr(){return Clr;}
 
  void                       DeleteAll();
  void                       Save();
  void                       FillH(double ModulationDegree=1.0);
  void                       ClusterViewComp(int NCluster);
  
  void GetResults_SIMPLE(double &Mu, double &MuErr, double &Phi, double &PhiError, double &CutEfficiency, double &FitProbabilityCut);
  void GetResults_CUT(double &Mu, double &MuErr, double &Phi, double &PhiError, double &CutEfficiency, double &FitProbabilityCut);

 private:
  TTree*                     t;
  TFile *UpdateFile;

  std::vector<double>        RealTheta;
  int Clusterdim;
  int Channel[1500];
  int Charge[1500];
  int DigiCharge[1500];
  double Phi;
  int                        Nentries,MixID,Clr ;
  double                     ClusterPulse;
  double                     Pitch;
  
  TDimension*                Dimension;
  TRandom*                   rnd;
  TGasMixture*               Mixture;
  TGasMixture*               MixtureRif;
  TDetector*                 myDetector;
  TReadout                   Readout;             
  double                     a,b,c,ElADC,Avalance,WIonization;//MANCA C!
  TString                    keyName;
  int                        NStatistic;
  double                     XShift;
  double                     YShift;
  double                     XIn;
  double                     YIn;
  
  //  TFile*                      F;
  TTree*                     treef;
  TH1D*                      ThirdMomHi;
  TH1D*                      SecondRatioHi;
  TH1D*                      ModulationCutH;
  TH1D*                      ModulationSimpleH;
  TF1*                       fufi1;
  TH1D*                      ClusterDim;
  
  TH1D*                      RawClusterDimH;
  TH1D*                      RecordDimH;
  TH1D*                      ModulationFirstStep;
  TH1D*                      RecordClusterDimH;
  double                     Pressure;
};
#endif
