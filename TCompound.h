#ifndef TCOMPOUND_HH
#define TCOMPOUND_HH

#define COMPOUNDS_FILE "./COMPOUNDS.DAT"
#include "MC.h"
/// \brief Class describing a chemical compound.
class TCompound {
 public:
  TCompound(TString COMPOUND_NAME, Bool_t VERBOSE=1);
  ~TCompound();
  inline TString              GetCompoundName()              {return CompoundName;}
  inline Double_t             GetFraction()                  {return Fraction;}
  inline void                 SetFraction(Double_t FRACTION) {Fraction = FRACTION;}
  inline Int_t                GetnElementsInCompound()       {return nElementsInCompound;}
  inline std::vector<TString> GetElementsInCompound()        {return ElementsInCompound;} 
  inline std::vector<Int_t>   GetnAtomsInCompound()          {return nAtomsInCompound;}
   //modifica calcolo WI
  inline double               GetCompoundIonsNumber()       {return  NIonsT;}
  inline double               GetCompoundStoppingPower()     {return  StoppingP;}
  
 private:
  Bool_t               Verbose;
  TString              CompoundName;
  Double_t             Fraction;
  Int_t                nElementsInCompound;
  std::vector<TString> ElementsInCompound;
  std::vector<Int_t>   nAtomsInCompound;
  //modifica calcolo WI
  double               NIonsT;
  double               StoppingP;
};

#endif
