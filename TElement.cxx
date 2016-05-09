#include "TElement.h"

//*******************************************************************************
/// \brief Basic constructor.                                                   *
//*******************************************************************************
TElement::TElement(TString ELEMENT_NAME, TRandom *RND, Bool_t VERBOSE)
{
  Verbose = VERBOSE;
  rnd = RND;
  if (IdentifyElement(ELEMENT_NAME)){
    std::cout << "Exiting..." << std::endl;
    return;
  }
  if (EvaluatePhotoelectricCrossSection()){
    std::cout << "Exiting..." << std::endl;
    return;
  }
  EvaluateMeanIonizationPotential();
  Density = AtomicWeight/MOLAR_VOLUME;  
  // The Density is re-calculated in TGasMixture (line 92) as partial density and Set here (SetDensity)
}


//*******************************************************************************
/// \brief Basic destructor.                                                    *
//*******************************************************************************
TElement::~TElement()
{

}


//*******************************************************************************
/// \brief Retrieves basic informations on the element.                         *
/// Reads the file containing element informations and retrieves the chemical   *
/// symbol, atomic number, and so on...                                         *
//*******************************************************************************
Int_t TElement::IdentifyElement(TString ELEMENT_NAME)
{
  // Open the file containing the element informations and skip the file header.
  ifstream ElementsFile;
  ElementsFile.open(ELEMENTS_FILE_NAME);
  TString Dummy = "";
  while (Dummy != "BEGIN:") ElementsFile >> Dummy;
  // Read the element properties.
  ChemicalSymbol = "NONE";
  while (!ElementsFile.eof()){
    ElementsFile >> ChemicalSymbol;
    ElementsFile >> AtomicNumber;
    ElementsFile >> AtomicWeight;
    ElementsFile >> kEdge;
    ElementsFile >> FluorescenceYield;
    if (ChemicalSymbol == ELEMENT_NAME){
      ElementsFile.close();
      if (Verbose){
	std::cout << std::endl;
	std::cout << "Chemical symbol:    " << ChemicalSymbol << std::endl;
	std::cout << "Atomic number:      " << AtomicNumber << std::endl;
	std::cout << "Atomic weight:      " << AtomicWeight << std::endl;
	std::cout << "kEdge (keV):        " << kEdge << std::endl;
	std::cout << "Fluorescence yield: " << FluorescenceYield << std::endl;
      }
      return 0;
    }
  }
  // If the element is not available...
  std::cout << ELEMENT_NAME << " is not available. Please edit by hand the file " <<
    ELEMENTS_FILE_NAME << "." << std::endl;
  ElementsFile.close();
  return -1;
}


//*******************************************************************************
/// \brief Evaluates the photoelectric cross section.                           *
/// Reads the files containing the photoelectric cross section for a given      *
/// energy grid (one per element, taken from the NIST database) and fits the    *
/// data with a modified power law. Correct treating of lines to be             *
/// implemented, yet.                                                           *
/// Energy to be provided in keV, cross section returned in cm^2/g.             *
//*******************************************************************************
Int_t TElement::EvaluatePhotoelectricCrossSection()
{
  Double_t Energy[N_ENERGY_BINS];
  Double_t PhotoelectricCrossSection[N_ENERGY_BINS];
  ifstream DataFile;
  TString DataFileName = "./DATA/" + ChemicalSymbol + ".txt";
  DataFile.open(DataFileName);
  TString Dummy = "";
  // Skip the stuff before the energy bins.
  while (Dummy != "BEGIN:") DataFile >> Dummy;
  for (Int_t i=0; i<N_ENERGY_BINS; i++){
    DataFile >> Energy[i] >> Dummy >> Dummy >> PhotoelectricCrossSection[i]
	     >> Dummy >> Dummy >> Dummy >> Dummy;
    Energy[i] *= 1000.0;
  }
  if (Verbose){
    std::cout << std::endl << "Photoelectric cross section (cm^2/g)..." << std::endl;
    for (Int_t i=0; i<N_ENERGY_BINS; i++){
      std::cout << Energy[i] << "   " << PhotoelectricCrossSection[i] << std::endl;
    }
  }
  PhotoelectricCrossSectionGraph =
    new TGraph(N_ENERGY_BINS, Energy, PhotoelectricCrossSection);
  return 0;
}


//*******************************************************************************
/// \brief Plots the photoelectric cross section.                               *
//*******************************************************************************
void TElement::PlotPhotoelectricCrossSection()
{
  TCanvas *PhotoelectricCrossSectionCanvas =
    new TCanvas("PhotoelectricCrossSection", "PhotoelectricCrossSection");
  PhotoelectricCrossSectionCanvas->SetFillColor(10);
  gPad->SetLogx();
  gPad->SetLogy();


  PhotoelectricCrossSectionGraph->SetMarkerStyle(26);
  PhotoelectricCrossSectionGraph->SetMarkerSize(0.7);
  PhotoelectricCrossSectionGraph->SetMarkerColor(1);
  PhotoelectricCrossSectionGraph->GetXaxis()->SetTitle("Photon energy (keV)");
  PhotoelectricCrossSectionGraph->GetYaxis()->SetTitle("Photoelectric cross section (cm^2/g)");
  PhotoelectricCrossSectionGraph->Draw("ap");
}


//*******************************************************************************
/// \brief Returns the photoelectric cross section for a given energy.          *
/// The number is taken from the fit to data with a modified power law.         *
/// Energy to be provided in keV, cross section returned in cm^2/g.             *
//*******************************************************************************
Double_t TElement::GetPhotoelectricCrossSection(Double_t ENERGY)
{
  if ((ENERGY < MIN_PHOTOELECTRIC_DATA_ENERGY) || (ENERGY > MAX_PHOTOELECTRIC_DATA_ENERGY)){
    std::cout << "Warning: TElement::GetPhotoelectricCrossSection." << std::endl;
    std::cout << "Photoelectric cross section data only available in the " <<
      MIN_PHOTOELECTRIC_DATA_ENERGY << "-" << MAX_PHOTOELECTRIC_DATA_ENERGY << " keV range."
	      << std::endl << "Returning -1." << std::endl;
    return -1;

  }
  // It seems that TGraph::Eval(x) fails while trying to retrieve the value of the last point.
  if (ENERGY == MAX_PHOTOELECTRIC_DATA_ENERGY) ENERGY *= 0.999999999999999;
  return PhotoelectricCrossSectionGraph->Eval(ENERGY);
}


//*******************************************************************************
/// \brief Returns the mean ionization potential.                               *
/// This parametrization is from Berger and Seltzer (1964); see Joy's book      *
/// for the reference. The ionization potential is returned in keV.             *
//*******************************************************************************
void TElement::EvaluateMeanIonizationPotential()
{
  MeanIonizationPotential = 0.001*(9.76*AtomicNumber + 58.5/pow((double)AtomicNumber, 0.19));
}


//*******************************************************************************
/// \brief Returns the Rutherford screening factor at a given anergy.           *
/// This parametrization is from Bishop (1976); see Joy's book for the          *
/// reference.                                                                  *
/// Energy to be provided in keV (accurate down to 50 eV).                      *
//*******************************************************************************
Double_t TElement::GetRutherfordScreeningFactor(Double_t ENERGY)
{
  if (ENERGY < MIN_ANALYTIC_CROSS_SECTIONS_ENERGY){
    std::cout << "Warning: TElement::GetRutherfordScreeningFactor." << std::endl;
    std::cout << "No data available for evaluating screening factor below " <<
      MIN_ANALYTIC_CROSS_SECTIONS_ENERGY << " keV."
	      << std::endl << "Returning -1." << std::endl;
    return -1;
  }
  return (3.4E-3)*pow((double)AtomicNumber, 0.67)/ENERGY;
}


//*******************************************************************************
/// \brief Returns the Rutherford integral cross ection at a given energy.      *
/// This parametrization is given to Newbury and Myklebust (1981); see          *
/// Joy's book for the reference.                                               *
/// Energy to be provided in keV, cross section returned in cm^2/atom           *
/// (accurate down to 50 eV).                                                   *
//*******************************************************************************
Double_t TElement::GetRutherfordTotalCrossSection(Double_t ENERGY)
{
  if (ENERGY < MIN_ANALYTIC_CROSS_SECTIONS_ENERGY){
    std::cout << "Warning: TElement::GetRutherfordTotalCrossSection." << std::endl;
    std::cout << "No data available for evaluating Rutherford cross section below " <<
      MIN_ANALYTIC_CROSS_SECTIONS_ENERGY << " keV."
	      << std::endl << "Returning -1." << std::endl;
    return -1;
  }
  Double_t alpha = GetRutherfordScreeningFactor(ENERGY);
  return (5.21E-21)*(pow((double)AtomicNumber, 2.0)/pow(ENERGY, 2.0))*
    (4*kPI/(alpha*(1 + alpha)))*pow(((ENERGY + 511)/(ENERGY + 1024)), 2.0);
}




//*******************************************************************************
/// \brief Returns the Mott integral cross section at a given energy.           *
/// This parametrization is given to Browning (1992); see Joy's book            *
/// for the reference.                                                          *
/// Energy to be provided in keV, cross section returned in cm^2/atom.          *
/// (accurate down to 50 eV).                                                   *
//*******************************************************************************
Double_t TElement::GetMottTotalCrossSection(Double_t ENERGY)
{
  if (ENERGY < MIN_ANALYTIC_CROSS_SECTIONS_ENERGY){
    std::cout << "Warning: TElement::GetMottTotalCrossSection." << std::endl;
    std::cout << "No data available for evaluating Mott cross section below " <<
      MIN_ANALYTIC_CROSS_SECTIONS_ENERGY << " keV."
	      << std::endl << "Returning -1." << std::endl;
    return -1;
  }
  Double_t u = log10(8*ENERGY*pow((double)AtomicNumber, -1.33));
  return (4.7E-18)*(pow((double)AtomicNumber, 1.33) + 0.032*pow((double)AtomicNumber, 2.0))/
    ((ENERGY + 0.0155*pow((double)AtomicNumber, 1.33)*pow((double)ENERGY, 0.5))*
     (1 - 0.02*pow((double)AtomicNumber, 0.5)*exp(-u*u)));
}


//*******************************************************************************
/// \brief Returns the mean free path for elastic collisions at a given energy. *
/// Energy to be provided in keV, density in g/cm^3; mode can be either         *
/// RUTHERFORD or MOTT. The mean free path is returned in cm.                   *
//*******************************************************************************
Double_t TElement::GetElasticMeanFreePath(Double_t ENERGY, Double_t DENSITY, TString MODE)
{
  if (ENERGY < MIN_ANALYTIC_CROSS_SECTIONS_ENERGY){
    std::cout << "Warning: TElement::GetElasticMeanFreePath." << std::endl;
    std::cout << "No data available for evaluating mean free path below " <<
      MIN_ANALYTIC_CROSS_SECTIONS_ENERGY << " keV."
	      << std::endl << "Returning -1." << std::endl;
    return -1;
  }
  Double_t MeanFreePath = AtomicWeight/(AVOGADRO_NUMBER*DENSITY);
  if (MODE == "RUTHERFORD"){
    MeanFreePath /= GetRutherfordTotalCrossSection(ENERGY);
    return MeanFreePath;
  }
  else if (MODE == "MOTT"){
    MeanFreePath /= GetMottTotalCrossSection(ENERGY);
    return MeanFreePath;
  }
  else{
    std::cout << "Warning: TElement::GetElasticMeanFreePath." << std::endl;
    std::cout << "Invalid mode in evaluating the elastic mean free path."
	      << std::endl << "Returning -1." << std::endl;
    return -1;
  }
}


//*******************************************************************************
/// \brief Returns a random scattering angle at a given energy.                 *
/// Energy to be provided in keV; mode can be either RUTHERFORD or MOTT.        *
/// The cosine of the angle is returned.                                        *
//*******************************************************************************
Double_t TElement::GetScatteringAngle(Double_t ENERGY, TString MODE)
{
  if (ENERGY < MIN_ANALYTIC_CROSS_SECTIONS_ENERGY){
    std::cout << "Warning: TElement::GetScatteringAngle." << std::endl;
    std::cout << "No data available for evaluating scattering angles below " <<
      MIN_ANALYTIC_CROSS_SECTIONS_ENERGY << " keV."
	      << std::endl << "Returning -1." << std::endl;
    return -1;
  }//acos
  if (MODE == "RUTHERFORD"){
    return acos(1. - (2.*GetRutherfordScreeningFactor(ENERGY)*rnd->Uniform())/
	    (1. + GetRutherfordScreeningFactor(ENERGY) - rnd->Uniform()));
  }//acos
  else if (MODE == "MOTT"){
    Double_t CorrectionFactor = 2.2 - ((92.0 - AtomicNumber)/92.0);
    return acos(1. - (2.*GetRutherfordScreeningFactor(ENERGY)*pow((rnd->Uniform()), CorrectionFactor))/
	    (1. + GetRutherfordScreeningFactor(ENERGY) - rnd->Uniform())); 
  }
  else{
    std::cout << "Warning: TElement::GetScatteringAngle." << std::endl;
    std::cout << "Invalid mode in evaluating the scattering angle."
	      << std::endl << "Returning -1." << std::endl;
    return -1;
  }
}


//*******************************************************************************
/// \brief Returns the stopping power at a given energy.                        *
/// This parametrization is from Joy and Luo (1989); see Joy's book for the     *
/// reference.                                                                  *
/// Energy to be provided in keV, outcome in keV/cm.                            *
//*******************************************************************************
Double_t TElement::GetStoppingPower(Double_t ENERGY)
{
  if (ENERGY < MIN_ANALYTIC_CROSS_SECTIONS_ENERGY){
    std::cout << "Warning: TElement::GetStoppingPower." << std::endl;
    std::cout << "No data available for evaluating stopping power below " <<
      MIN_ANALYTIC_CROSS_SECTIONS_ENERGY << " keV."
	      << std::endl << "Returning -1." << std::endl;
    return -1;
  }
  return 78500*AtomicNumber*log(1.166*(ENERGY + 0.85*MeanIonizationPotential)/
				MeanIonizationPotential)/(ENERGY*AtomicWeight)*Density; 
}

//validazione Stoppingpower,OK//
void TElement:: PlotStoppingPower()
{
  //TRandom *rnd=new TRandom();
  //TElement *E=new TElement(ELEMENT_NAME,rnd,1);
  SetDensity(10);
  const int N=30; 
  Double_t Stoppingpower[N];
  Double_t Energy[N];

  for (int i=0; i<30; i++)
    {
      Energy[i]=1.0*(i+1);
      Double_t StoppingP=GetStoppingPower(Energy[i]);
      Double_t D=GetDensity();
      Stoppingpower[i]=StoppingP/D;
    }
  TGraph *GrafStoppingpower = new TGraph(N,Energy,Stoppingpower );
  //TCanvas *hCanvas = new TCanvas("stoppingpower"," stoppingpower", 500, 500, 700, 700);
  GrafStoppingpower->SetMinimum(1);
  GrafStoppingpower->Draw("al*");
  std::cout<<Stoppingpower[1]<<std::endl;
  GrafStoppingpower->GetXaxis()->SetTitle("Energy (Kev)");
  GrafStoppingpower->GetYaxis()->SetTitle("Stoppingpower (Kev*cm^2/gr)");
  
}
//validazione cross section
void TElement:: PlotMottcrossSection()
{const int N=100; 
  Double_t MottCS[N];
  Double_t En[N];
for (int i=0; i<100; i++)
  {En[i]=(i+1.0)/10.0;
  MottCS[i]=GetMottTotalCrossSection(En[i]);
  }
TGraph *GrafMottCS = new TGraph(N,En,MottCS );
GrafMottCS->SetMinimum(0.01);
GrafMottCS->Draw("al*");
}
void TElement::PlotScatteringA()//comportamento anomalo a piccoli angoli?
{TH1D *h1=new TH1D("h","distrAng",1000,0.01,240);
 const int N=10000; 
 Double_t Angle[N];
 for (int i=0; i<10000; i++){
   Angle[i]=360*(GetScatteringAngle(1,"MOTT"))/(kPI*2);
   h1->Fill(Angle[i]);}
 h1->Draw();
}
