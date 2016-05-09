#ifndef FASTSIMULATOR_HH
#define FASTSIMULATOR_HH

#include "MC.h"
#include "TExperiment.h"

class FastSimulator
{
 public:
  FastSimulator(bool draw=1);
  ~FastSimulator()
    {
      bigfile->Close();
      delete rnd;
      filetree->Close();
      delete Experiment;
    }
  
  void SetOption(int MixID,double Thickness,double Pressure,int PixmapS = 50,TString BigName="BigFile.root");  
  void SetEnergy(double Energy,  TString TreeName = "Eventstree");
  
  inline void SetPolarizationFraction(double frac) {PolarizationFraction=frac;}
  inline void SetPolarizationAngle(double deg) {PolarizationAngle=deg;}

  //  void SetSpectrum(TH1D *spectrum){Spectrum=spectrum;}
  TH1F *AngularResp();
  TF1 *PhiDistribution();
  void ProcessEvents(int N);
  void GetResults(double &Mu,double &MuE, double &AngleErr);
 private:
  double PI;
  bool Draw;
  TCanvas *a;
  TString AnalysisName,DataName;
  
  TFile *filetree,*bigfile;
  TTree *treemix,*bigtree;
  TExperiment *Experiment;
  TRandom3 *rnd;
  //  TH1D *Spectrum;
  int MIXID,PIXMAP;
  double  THICKNESS,PRESSURE;
  TString BigName;
  double PolarizationFraction;
  double PolarizationAngle;
  double ModulationFactor,ModulationFactorError,AngleError;
};
#endif
