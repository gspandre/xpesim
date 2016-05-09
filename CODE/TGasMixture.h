#ifndef TGASMIXTURE_HH
#define TGASMIXTURE_HH

#include "TElement.h"
#include "TCompound.h"
#include <TF1.h>
#include "TDimension.h"
#include "TText.h"
#define MIXTURES_FILE_NAME "./MIXTURES.DAT"

/// \brief Class decribing a gas mixture.
class TGasMixture {
 public:
  TGasMixture(Int_t MixtureID, TRandom *RND,TDimension *Dimension,Bool_t VERBOSE=1);
  ~TGasMixture()
    {
      for(int i = 0; i<Compounds.size();i++)
	delete Compounds[i];
      for(int i = 0; i<Elements.size();i++)
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
  void                           PlotEfficiency(Double_t ENERGY_MIN=MIN_PHOTOELECTRIC_DATA_ENERGY,
						Double_t ENERGY_MAX=MAX_PHOTOELECTRIC_DATA_ENERGY);
						//Int_t N_POINTS=75, Double_t GAP_THICKNESS=1.0);
  Double_t                       GetElectronRange(Double_t ENERGY);
  void                           PlotElectronRange();
  TElement*                      GetConvertingElement(Double_t ENERGY);
  Double_t                       GetElasticMeanFreePath(Double_t ENERGY, TString MODE="MOTT");
  TElement*                      GetScatteringElement(Double_t ENERGY, TString MODE="MOTT");
  Double_t                       GetStoppingPower(Double_t ENERGY);
  void                           PlotPhotoelectricCrossSection();
  inline TString                 GetName(){return  MixtureName;}
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
  
  
  void CheckFractions();

  Double_t GAP_THICKNESS;
};

#endif
