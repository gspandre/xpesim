#ifndef TGASMIXTURE_HH
#define TGASMIXTURE_HH

#include "MC.h"
#include "TElement.h"
#include "TCompound.h"
#include "TDimension.h"

#define MIXTURES_FILE_NAME "./MIXTURES.DAT"

/// \brief Class decribing a gas mixture.
class TGasMixture {
 public:
  TGasMixture(Int_t MixtureID, TRandom *RND,TDimension *Dimension,Bool_t VERBOSE=1);
  ~TGasMixture()
    {
      for(int i = 0; i<(int)Compounds.size();i++)
	delete Compounds[i];
      for(int i = 0; i<(int)Elements.size();i++)
	delete Elements[i];
      delete ElectronRangeFunction;
    }

  inline Double_t                GetPressure()       {return Pressure;}
  inline Double_t                GetDensity()        {return Density;}
  inline Double_t                GetDiffusionSigma() {return DiffusionSigma;}
  inline Double_t                GetWIonization()    {return WIonization;}
  inline Int_t                   GetnCompounds()     {return nCompounds;}
  inline std::vector<TCompound*> GetCompounds()      {return Compounds;}
  inline Int_t                   GetnElements()      {return nElements;}
  inline std::vector<TElement*>  GetElements()       {return Elements;}
  Double_t                       GetPhotoelectricCrossSection(Double_t ENERGY);
  Double_t                       GetAbsorptionLenght(Double_t ENERGY);
  Double_t                       GetEfficiency(Double_t ENERGY); //Double_t GAP_THICKNESS);
  Double_t                       GetElectronRange(Double_t ENERGY);
  void                           PlotElectronRange();
  TElement*                      GetConvertingElement(Double_t ENERGY);
  Double_t                       GetElasticMeanFreePath(Double_t ENERGY, TString MODE="MOTT");
  TElement*                      GetScatteringElement(Double_t ENERGY, TString MODE="MOTT");
  Double_t                       GetStoppingPower(Double_t ENERGY);
  void                           PlotPhotoelectricCrossSection();
  inline TString                 GetName(){return  MixtureName;}
  std::vector<TGraph*>           GetAbsProbGraph(Double_t ENERGY);

 private:
  TString                        MixtureName;
  Bool_t                         Verbose;
  TRandom*                       rnd;
  Double_t                       Pressure;
  Double_t                       DiffusionSigma;
  Double_t                       WIonization;
  Double_t                       Density;
  Int_t                          nCompounds;
  std::vector<TCompound*>        Compounds;
  Int_t                          nElements;
  std::vector<TElement*>         Elements;
  TF1*                           ElectronRangeFunction;
  std::vector<TGraph*>           AbsorptionProbabilityGr;
  
  void CheckFractions();

  Double_t GAP_THICKNESS;
};

#endif
