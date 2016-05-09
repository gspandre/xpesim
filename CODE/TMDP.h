#ifndef TMDP_HH
#define TMDP_HH
#include "TMarker.h"
#include "TPaveText.h"
#include "TMarker.h"
#include "TF1.h"
#include <TGraph.h>
#include <iostream>
#include <fstream>
#include <TCanvas.h>
#include "TGraph.h"
#include "pair.h"
#include "TFile.h"
#include "TH1D.h"
#include <vector>
#include "TGaxis.h"
#include "TText.h"
#include "TPolyLine.h"
#include "TPaveLabel.h"
#include "TAxis.h"
class TMDP
{
 public:
  TMDP();
  ~TMDP() {;} 
  double                                MDP(TGraph *aeff, TF1 *flux, TGraph *Effgas,double InEne,double FinEne,double T);
  TF1                                   *HerculesSpectrum();
  TF1                                   *Spectrum(double SpectralIndex, double FluxinCrab);
  TF1                                   *Spectrum2(double SpectralIndex,  double Constant);
  TF1                                   *CrabSpectrum();
  
  TH1D  *Spectrum(TF1 *sourceSpectrum, double emin=1.0,double emax=30.0); 
  TH1D  *Spectrum(TF1 *sourceSpectrum, TGraph *aeff, double emin=1.0,double emax=30.0); 
  TH1D  *Spectrum(TF1 *sourceSpectrum, TGraph *aeff, TGraph *window, double emin=1.0,double emax=30.0); 
  TH1D  *Spectrum(TF1 *sourceSpectrum, TGraph *aeff, TGraph *window, TGraph *mixEff, double emin=1.0,double emax=30.0); 

  
  void                                  DrawFlux(TString Source);
  void                                  PlotDiffMdp();
  void                                  ReadData(TString FileName);
  void                                  ReadDatafromGraph(TGraph *MuvsEn, TGraph *EffvsEn,int  MIXID,double THICKNESS,double PRESSURE);
  inline TGraph                         *GetEfficiencyGraph(){return  EfficiencyGraph;}
  
  //TGraph                                *MirrorXEUS(double  WindowEfficiency);
  //TGraph                                *MirrorSIMBOLIXPLUS(double  WindowEfficiency);
  TGraph                                *MirrorGraph(double  WindowEfficiency);
  
  inline int                            GetN_bins(){return N_bins;}
  void                                  MDPCanvas();
  void                                  FillEfficiencyGraph();
  void                                  FillMuGraph();
  TGraph                                *FillMDPGraph(TF1 *Source,double Time = 3600.*24.);
  inline TString                        GetMixtureName(){return MixtureName; }
  void                                  CompactPlot();
  double                                EvalFluxParameter(double SpectralIndex,double FluxCrab,double MinEnergy=2.,double MaxEnergy = 10.);
  double                  EvalMDP2_10kev(double SpectralIndex,double FluxinCrab,double Time);
  void                    FluxMDPGraph(TString Name);
  void                   AddSource(double Flux,double Alpha,TString Name,double Time=24);
  void                   AddVariableSource(double Fluxmin,double Fluxmax,double Alphamin,double Alphamax,TString Name,int Marker=0,int color=2,double Time=24);
  void                   FindMW(TGraph *Graph);
  inline double          GetGraphMax(){return GraphMax;}
  inline double          GetGraphWsx(){return GraphWsx;}
  inline double          GetGraphWdx(){return GraphWdx;}
  inline double          GetGraphMaxX(){return GraphMaxX; }
  TGraph                 *MuEpGraph();
  double                 Rate(double InEne,double FinEne,TString Source="CRAB");
  double                 Rate(TF1 *flux, double InEne,double FinEne);
  double                 Rate(TH1D* flux, double InEne,double FinEne);

  inline TGraph          *GetWindowEfficiency(){return Windowgr; }
  inline double          GetMaximumMuEp(){return MaximumMuEp;}
  void                   SetMirror(TString MirrorN);
 private:
  TString MirrorName ;
  TGraph*  mu;
  TGraph*  Mug;
  TGraph*  Windowgr;
  // TGraph*  MirrorX;
  //TGraph*  MirrorS;
  TString MixtureName;

   

  TGraph*  EfficiencyGraph;
  TF1*     Her;
  TF1*     SpectrumI;
  double Pressure;
  double MinimumEnergy;
  double MaximumEnergy;
  double Thickness;
  int    MixID;
  int  N_bins;
  std::vector<double> EnergyNVector;
  std::vector<double> MuNVector;
  std::vector<double> TotalEfficiencyVector;
  /* double* EnergyVector;
  double* TotalEfficiency;
  double* MuVector;*/
  double GraphMax;
  double GraphWsx;
  double GraphWdx;
  double GraphMaxX;
  double MaximumMuEp;
};

#endif
