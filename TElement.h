#ifndef TELEMENT_HH
#define TELEMENT_HH

#include "MC.h"
#define MOLAR_VOLUME          22414.0
#define AVOGADRO_NUMBER       6.02E23
#define N_ENERGY_BINS         75
#define ELEMENTS_FILE_NAME    "./ELEMENTS.DAT"

#define MIN_PHOTOELECTRIC_DATA_ENERGY         1.0
#define MAX_PHOTOELECTRIC_DATA_ENERGY         1000.0
#define MIN_ANALYTIC_CROSS_SECTIONS_ENERGY    0.05

/// \brief Class describing a chemical element.
class TElement {
 public:
  TElement(TString ELEMENT_NAME, TRandom *RND, Bool_t VERBOSE=1);
  ~TElement();
  inline TString  GetChemicalSymbol()         {return ChemicalSymbol;}
  inline Int_t    GetAtomicNumber()           {return AtomicNumber;}
  inline Double_t GetAtomicWeight()           {return AtomicWeight;}
  inline void     SetDensity(Double_t DENSITY){Density = DENSITY;}
  inline Double_t GetDensity()                {return Density;}
  inline Double_t GetkEdge()                  {return kEdge;}
  inline Double_t GetFluorescenceYield()      {return FluorescenceYield;}
  inline Double_t GetMeanIonizationPotential(){return MeanIonizationPotential;}
  void            PlotPhotoelectricCrossSection();
  Double_t        GetPhotoelectricCrossSection(Double_t ENERGY);
  Double_t        GetRutherfordScreeningFactor(Double_t ENERGY);
  Double_t        GetRutherfordTotalCrossSection(Double_t ENERGY);
  Double_t        GetMottTotalCrossSection(Double_t ENERGY);
  Double_t        GetElasticMeanFreePath(Double_t ENERGY, Double_t DENSITY, TString MODE="MOTT");
  Double_t        GetScatteringAngle(Double_t ENERGY, TString MODE="MOTT");
  Double_t        GetStoppingPower(Double_t ENERGY);
  void            PlotStoppingPower(); 
  void            PlotMottcrossSection();
  void            PlotScatteringA();
 private:
  Bool_t          Verbose;
  TRandom*        rnd;
  TString         ChemicalSymbol;
  Int_t           AtomicNumber;
  Double_t        AtomicWeight;
  Double_t        Density;
  Double_t        kEdge;
  Double_t        FluorescenceYield;
  Double_t        MeanIonizationPotential;
  TCanvas*        PhotoelectricCrossSectionCanvas;
  TGraph*         PhotoelectricCrossSectionGraph;
  Int_t           IdentifyElement(TString ELEMENT_NAME);
  Int_t           EvaluatePhotoelectricCrossSection();
  void            EvaluateMeanIonizationPotential();
};

#endif
