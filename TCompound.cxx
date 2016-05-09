#include "TCompound.h"


//*******************************************************************************
/// \brief Basic constructor.                                                   *
//*******************************************************************************
TCompound::TCompound(TString COMPOUND_NAME, Bool_t VERBOSE)
{
  Verbose = VERBOSE;
  CompoundName = "NONE";
  // Open the file containing the compound informations and skip the file header..
  ifstream CompoundsFile;
  CompoundsFile.open(COMPOUNDS_FILE);
  TString Dummy = "";
  while (Dummy != "BEGIN:") CompoundsFile >> Dummy;
  // Read the compound properties.
  TString ElementTemp;
  Int_t nAtomsTemp;
  while (!CompoundsFile.eof()){
    CompoundsFile >> CompoundName;
    CompoundsFile >> nElementsInCompound;
    //modifica Calcolo WI
    CompoundsFile >>  NIonsT;  //total ionization/cm
    CompoundsFile >>  StoppingP;   // de/dx in ev/cm
    //
    for (Int_t i=0; i<nElementsInCompound; i++){
      CompoundsFile >> ElementTemp;
      ElementsInCompound.push_back(ElementTemp);
      CompoundsFile >> nAtomsTemp;
      nAtomsInCompound.push_back(nAtomsTemp);
    }
    if (CompoundName == COMPOUND_NAME){
      CompoundsFile.close();
      if (Verbose){
	std::cout << "Compound name:  " << CompoundName << std::endl;
	std::cout << "Number of elements: " << nElementsInCompound << std::endl;
	for (Int_t j=0; j<nElementsInCompound; j++){
	  std::cout << "Element: " << ElementsInCompound[j] << ", ";
	  std::cout << "number of atoms: " << nAtomsInCompound[j] << std::endl;
	}
      }
      return;
    }
    else{
      ElementsInCompound.clear();
      nAtomsInCompound.clear();
    }
  }
  // If the compound is not available...
  std::cout << COMPOUND_NAME << " is not available. Please edit by hand the file " <<
    COMPOUNDS_FILE << std::endl;
  CompoundsFile.close();
  CompoundName = "NONE";
  return;
}


//*******************************************************************************
/// \brief Basic destructor.                                                    *
//*******************************************************************************
TCompound::~TCompound()
{

}
