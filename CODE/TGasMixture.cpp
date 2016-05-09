#include "TGasMixture.h"

//*******************************************************************************
/// \brief Basic constructor.                                                   *
//*******************************************************************************

TGasMixture::TGasMixture(Int_t MixtureID, TRandom *RND,TDimension *Dimension,Bool_t VERBOSE)
{
 
  GAP_THICKNESS = (Dimension->GetZ_Drift()) - (Dimension->GetZ_Gem());
  // Initialize some variables.
  Verbose = VERBOSE;
  ifstream DataFile(MIXTURES_FILE_NAME);
  rnd = RND;
  TString  CompoundName = "";
  TString  DummyString  = "";
  Double_t Frac;
  // Skip file header.
  while (DummyString != "Sigma@1atm")      DataFile >> DummyString;
  // Go straight to the interesting mixture.
  for (Int_t i=0; i<MixtureID*7; i++) {
    DataFile >> DummyString;
    // std::cout<<DummyString<<std::endl;
  }
  // Read relevant properties of the compounds in the mixture.
  for (Int_t i=0; i<3; i++)
    {
      DataFile >> CompoundName;
      DataFile >> Frac;
      if (CompoundName == "-") continue;
      TCompound *myCompound = new TCompound(CompoundName, 0);
      myCompound->SetFraction(Frac);
      Compounds.push_back(myCompound);  
    }
  // DataFile >> Pressure;
  Pressure = Dimension->GetPressure();
  DataFile >> DiffusionSigma;
  // Scale the diffusion sigma according to the pressure.
  DiffusionSigma /= sqrt(Pressure); 
  DataFile.close();
  nCompounds = Compounds.size();
  CheckFractions();
  // Print out informations.
 
  std::cout << std::endl << "Retrieving mixture informations..." << std::endl;
  
  std::cout << "Number of compounds in the mixture: " << nCompounds << std::endl;
  double TotalStopPower = 0;
  double TotalSecondary = 0;
  
  for (Int_t i=0; i<nCompounds; i++)
    {
      std::cout << "Compound name: " << Compounds[i]->GetCompoundName()
		<< ", fraction: " << Compounds[i]->GetFraction() << std::endl;
      MixtureName+=Compounds[i]->GetCompoundName();
      
      MixtureName+=Compounds[i]->GetFraction();
      MixtureName+="-";
      //modifica calcolo WI
      TotalStopPower +=  Compounds[i]->GetFraction() * Compounds[i]->GetCompoundStoppingPower();
      TotalSecondary +=  Compounds[i]->GetFraction() * Compounds[i]->GetCompoundIonsNumber();
      //
    }
  //modifica calcolo WI
  WIonization =   TotalStopPower/ TotalSecondary;
  //
  std::cout<<"Gap Thickness = "<<GAP_THICKNESS<<std::endl;
  std::cout<<"Mixture pressure: " << Pressure << std::endl;
  std::cout<<"Diffusion sigma (from the value @ 1atm): "   
	   << DiffusionSigma << " um (for 1 cm drift)" << std::endl;
  std::cout <<"W ionization: " << WIonization << " eV" << std::endl;
  MixtureName+=Pressure;
  MixtureName+=" atm, ";
  MixtureName+=GAP_THICKNESS;
  MixtureName+=" cm ";
  
  
  // Fill the list of elements in the mixture and evaluate relevant properties.
  Density   = 0.0;
  nElements = 0;
  for (Int_t i=0; i<nCompounds; i++){
    if (Verbose) std::cout << std::endl << "Retrieving informations about "
			   << Compounds[i]->GetCompoundName() << "..." << std::endl;
    for (Int_t j=0; j<Compounds[i]->GetnElementsInCompound(); j++){
      TElement *myElement = new TElement(Compounds[i]->GetElementsInCompound()[j], rnd, 0);
      Density += myElement->GetAtomicWeight()*Compounds[i]->GetnAtomsInCompound()[j]*
	Compounds[i]->GetFraction()*Pressure/(MOLAR_VOLUME);
      Bool_t AlreadyThere = 0;
      for (Int_t k=0; k<nElements; k++){
	if (myElement->GetChemicalSymbol() == Elements[k]->GetChemicalSymbol()){
	  AlreadyThere = 1;
	  Elements[k]->SetDensity(myElement->GetAtomicWeight()*Compounds[i]->GetnAtomsInCompound()[j]*
				  Compounds[i]->GetFraction()*Pressure/(MOLAR_VOLUME) +
				  Elements[k]->GetDensity());
	}
      }
      if (!AlreadyThere){
	nElements += 1;
	Elements.push_back(myElement);
	Elements[nElements-1]->SetDensity(myElement->GetAtomicWeight()*
					  Compounds[i]->GetnAtomsInCompound()[j]*
					  Compounds[i]->GetFraction()*Pressure/(MOLAR_VOLUME));
      }
    }
  }
  // Print out informations.
  if (Verbose){
    std::cout << std::endl << "Density of the mixture: " << Density << " g/cm^3" << std::endl;
    std::cout << "Number of elements in the mixture: " << nElements << std::endl;
    for (Int_t i=0; i<nElements; i++){
      std::cout << "Element: " << Elements[i]->GetChemicalSymbol()
		<< ", partial density: " << Elements[i]->GetDensity() << " g/cm^3" << std::endl;
    }
  }

  // Evaluate electron range as a function of the energy.
  ElectronRangeFunction = new TF1("RangeLaw","[0]*pow(10, (-5.1 + 1.358*log10(x) + 0.215*(pow(log10(x),2)) -0.043*(pow(log10(x),3))))",MIN_PHOTOELECTRIC_DATA_ENERGY,MAX_PHOTOELECTRIC_DATA_ENERGY);
  ElectronRangeFunction->SetNameTitle("Electron Range","Electron Range");
  ElectronRangeFunction->FixParameter(0, 10.0/Density); 
 
}


//*******************************************************************************
/// \brief Checks that the compound fractions are normalized to unity.          *
/// If that's not the case, the fractions are renormalized.                     *
//*******************************************************************************
void TGasMixture::CheckFractions()
{
  Double_t FractionsSum = 0;
  for (Int_t i=0; i<nCompounds; i++) FractionsSum += Compounds[i]->GetFraction();
  if ((FractionsSum < 0.99) || (FractionsSum > 1.01)){
    std::cout << std::endl << "Warning: sum of components fractions (";
    for (Int_t i=0; i<nCompounds-1; i++) std::cout << Compounds[i]->GetFraction() << ", ";
    std::cout << Compounds[nCompounds-1]->GetFraction() << ") is different from unity." << std::endl;
    for (Int_t i=0; i<nCompounds; i++)
      Compounds[i]->SetFraction(Compounds[i]->GetFraction()/FractionsSum);
    std::cout << "Setting fractions to (";
    for (Int_t i=0; i<nCompounds-1; i++) std::cout << Compounds[i]->GetFraction() << ", ";
    std::cout << Compounds[nCompounds-1]->GetFraction() << ") instead." << std::endl;
  }
}


//*******************************************************************************
/// \brief Returns the Photoelectric cross section at a given energy.           *
/// Energy to be provided in keV, cross section returned in cm^2/g.             *
//*******************************************************************************
Double_t TGasMixture::GetPhotoelectricCrossSection(Double_t ENERGY)
{
  if ((ENERGY < MIN_PHOTOELECTRIC_DATA_ENERGY) || (ENERGY > MAX_PHOTOELECTRIC_DATA_ENERGY)){
    std::cout << "Warning: TGasMixture::GetPhotoelectricCrossSection." << std::endl;
    std::cout << "Photoelectric cross section data only available in the " <<
      MIN_PHOTOELECTRIC_DATA_ENERGY << "-" << MAX_PHOTOELECTRIC_DATA_ENERGY << " keV range."
	      << std::endl << "Returning -1." << std::endl;
    return -1;

  }
  Double_t CrossSection = 0;
  for (Int_t i=0; i<nElements; i++){
    CrossSection += Elements[i]->GetPhotoelectricCrossSection(ENERGY)*
      Elements[i]->GetDensity()/Density;
  }
  return CrossSection;
}


//*******************************************************************************
/// \brief Returns the photoelectric absorption lenght at a given energy.       *
/// Energy to be provided in keV, absorption lenght returned in cm.             *
//*******************************************************************************
Double_t TGasMixture::GetAbsorptionLenght(Double_t ENERGY)
{
  if ((ENERGY < MIN_PHOTOELECTRIC_DATA_ENERGY) || (ENERGY > MAX_PHOTOELECTRIC_DATA_ENERGY)){
    std::cout << "Warning: TGasMixture::GetAbsorptionLenght." << std::endl;
    std::cout << "Absorption lenght data only available in the " <<
      MIN_PHOTOELECTRIC_DATA_ENERGY << "-" << MAX_PHOTOELECTRIC_DATA_ENERGY << " keV range."
	      << std::endl << "Returning -1." << std::endl;
    return -1;
  }
  return 1.0/(Density*GetPhotoelectricCrossSection(ENERGY));
}


//*******************************************************************************
/// \brief Returns the detection efficiency at given energy and gap thickness.  *
/// Energy to be provided in keV.                                               *
//*******************************************************************************
Double_t TGasMixture::GetEfficiency(Double_t ENERGY)
{
  if ((ENERGY < MIN_PHOTOELECTRIC_DATA_ENERGY) || (ENERGY > MAX_PHOTOELECTRIC_DATA_ENERGY)){
    std::cout << "Warning: TGasMixture::GetEfficiency." << std::endl;
    std::cout << "Efficiency data only available in the " <<
      MIN_PHOTOELECTRIC_DATA_ENERGY << "-" << MAX_PHOTOELECTRIC_DATA_ENERGY << " keV range."
	      << std::endl << "Returning -1." << std::endl;
    return -1;
  }
  return (1 - exp(-GAP_THICKNESS/GetAbsorptionLenght(ENERGY)));
}


//*******************************************************************************
/// \brief Plots the detection efficiency and related quantities.               *
//*******************************************************************************
void TGasMixture::PlotEfficiency(Double_t ENERGY_MIN, Double_t ENERGY_MAX)// Int_t N_POINTS,Double_t GAP_THICKNESS)
{
  TCanvas *EfficiencyCanvas =
    new TCanvas("Efficiency", "Efficiency", 10, 10, 900, 800);
  EfficiencyCanvas->SetFillColor(10);
  EfficiencyCanvas->Divide(2, 2);
  const Int_t nPoints = N_ENERGY_BINS;//N_POINTS;
  Double_t Energy[nPoints];
  Double_t PhotoelectricCrossSection[nPoints];
  Double_t AbsorptionProbability[(const Int_t)nElements][nPoints];
  Double_t AbsorptionLenght[nPoints];
  Double_t Efficiency[nPoints];
  if (Verbose){
    std::cout << std::endl << "Efficiency data:" << std::endl;
    std::cout << "Energy(keV)  PhCrossSection(cm^2/g) " 
	      << " AbsLenght(cm)  Eff(" << GAP_THICKNESS << " cm gap)" << std::endl;
  }
  for (Int_t i=0; i<nPoints; i++){
    //    Energy[i] = ENERGY_MIN + i*(ENERGY_MAX - ENERGY_MIN)/(nPoints-1);
    Energy[i] = ENERGY_MIN *pow(ENERGY_MAX/ENERGY_MIN,1.0*i/(nPoints-1));
    PhotoelectricCrossSection[i] = GetPhotoelectricCrossSection(Energy[i]);
    for (Int_t j=0; j<nElements; j++){
      AbsorptionProbability[j][i] =
	(Elements[j]->GetPhotoelectricCrossSection(Energy[i])*Elements[j]->GetDensity())/
	(GetPhotoelectricCrossSection(Energy[i])*Density);
    }
    AbsorptionLenght[i] = GetAbsorptionLenght(Energy[i]);
    Efficiency[i] = GetEfficiency(Energy[i]);
    if (Verbose) std::cout << Energy[i] << "  " << PhotoelectricCrossSection[i] << "  "
			   << AbsorptionLenght[i] << "  " << Efficiency[i] << std::endl;
  }
  EfficiencyCanvas->cd(1);
  gPad->SetLogx();
  gPad->SetLogy();
  TGraph *PhotoelectricCrossSectionGraph = new TGraph(nPoints, Energy, PhotoelectricCrossSection);
  PhotoelectricCrossSectionGraph->SetTitle("PhotoelectricCrossSection");
  PhotoelectricCrossSectionGraph->SetMarkerStyle(26);
  PhotoelectricCrossSectionGraph->SetMarkerSize(0.7);
  PhotoelectricCrossSectionGraph->SetMarkerColor(1);
  PhotoelectricCrossSectionGraph->GetXaxis()->SetTitle("Photon energy (keV)");
  PhotoelectricCrossSectionGraph->GetYaxis()->SetTitle("Photoelectric cross section (cm^2/g)");
  PhotoelectricCrossSectionGraph->GetYaxis()->SetTitleOffset(1.1);
  //PhotoelectricCrossSectionGraph->SetRange(1.,50.);
  PhotoelectricCrossSectionGraph->Draw("ap");
  EfficiencyCanvas->cd(2);
  gPad->SetLogx();
  gPad->SetLogy();
  TGraph *AbsorptionProbabilityGraphs[(const Int_t)nElements];
  for (Int_t i=0; i<nElements; i++){
    AbsorptionProbabilityGraphs[i] = new TGraph(nPoints, Energy, AbsorptionProbability[i]);
  }
  AbsorptionProbabilityGraphs[0]->SetTitle("AbsorptionProbabilities");
  AbsorptionProbabilityGraphs[0]->SetMarkerStyle(26);
  AbsorptionProbabilityGraphs[0]->SetMarkerSize(0.7);
  AbsorptionProbabilityGraphs[0]->SetMarkerColor(1);
  AbsorptionProbabilityGraphs[0]->GetXaxis()->SetTitle("Photon energy (keV)");
  AbsorptionProbabilityGraphs[0]->GetYaxis()->SetTitle("AbsorptionProbabilities");
  AbsorptionProbabilityGraphs[0]->GetYaxis()->SetTitleOffset(1.1);
  AbsorptionProbabilityGraphs[0]->SetMaximum(1.0);
  AbsorptionProbabilityGraphs[0]->SetMinimum(0.001);
  AbsorptionProbabilityGraphs[0]->SetName(Elements[0]->GetChemicalSymbol());
  AbsorptionProbabilityGraphs[0]->Draw("ap");
  for (Int_t i=1; i<nElements; i++){
    AbsorptionProbabilityGraphs[i]->SetMarkerStyle(26);
    AbsorptionProbabilityGraphs[i]->SetMarkerSize(0.7);
    AbsorptionProbabilityGraphs[i]->SetMarkerColor(i+1);
    AbsorptionProbabilityGraphs[i]->SetName(Elements[i]->GetChemicalSymbol());
    AbsorptionProbabilityGraphs[i]->Draw("p");
  }
  EfficiencyCanvas->cd(3);
  gPad->SetLogx();
  gPad->SetLogy();
  TGraph *AbsorptionLenghtGraph = new TGraph(nPoints, Energy, AbsorptionLenght);
  AbsorptionLenghtGraph->SetTitle("AbsorptionLenght");
  AbsorptionLenghtGraph->SetMarkerStyle(26);
  AbsorptionLenghtGraph->SetMarkerSize(0.7);
  AbsorptionLenghtGraph->SetMarkerColor(1);
  AbsorptionLenghtGraph->GetXaxis()->SetTitle("Photon energy (keV)");
  AbsorptionLenghtGraph->GetYaxis()->SetTitle("Absorption lenght (cm)");
  AbsorptionLenghtGraph->GetYaxis()->SetTitleOffset(1.1);
  AbsorptionLenghtGraph->Draw("ap");
  EfficiencyCanvas->cd(4);
  gPad->SetLogx();
  gPad->SetLogy();
  TGraph *EfficiencyGraph = new TGraph(nPoints, Energy, Efficiency);
  EfficiencyGraph->SetTitle("Efficiency");
  EfficiencyGraph->SetMarkerStyle(26);
  EfficiencyGraph->SetMarkerSize(0.7);
  EfficiencyGraph->SetMarkerColor(1);
  EfficiencyGraph->GetXaxis()->SetTitle("Photon energy (keV)");
  EfficiencyGraph->GetYaxis()->SetTitle("Efficiency");
  EfficiencyGraph->GetYaxis()->SetTitleOffset(1.1);
  EfficiencyGraph->Draw("ap");
}


//*******************************************************************************
/// \brief Returns the electron range as a function of the energy.              *
/// Energy to be provided in keV.                                               *
//*******************************************************************************
Double_t TGasMixture::GetElectronRange(Double_t ENERGY)
{
  return ElectronRangeFunction->Eval(ENERGY);
}


//*******************************************************************************
/// \brief Plots the electron range as a function of the energy.                *
//*******************************************************************************
void TGasMixture::PlotElectronRange()
{
  //
  const Int_t np = N_ENERGY_BINS;
  Double_t elenergy[np];
  Double_t ElasticFreePathMott[np];
  for (int i=0;i<N_ENERGY_BINS;i++)
    {
      elenergy[i] = 1+i;
      //elenergy[i] =MIN_PHOTOELECTRIC_DATA_ENERGY *pow(MAX_PHOTOELECTRIC_DATA_ENERGY/MIN_PHOTOELECTRIC_DATA_ENERGY,1.0*i/(np-1));
      ElasticFreePathMott[i] = 10000*GetElasticMeanFreePath(elenergy[i],"MOTT");//da cm a micron
    }
  //
  TCanvas *ElectronRangeCanvas = new TCanvas("ElectronRange", "ElectronRange",1200,800);
  ElectronRangeCanvas->SetFillColor(10);
  ElectronRangeCanvas->Divide(2,2);//
  ElectronRangeCanvas->cd(1);//
  ElectronRangeFunction->GetYaxis()->SetLimits(0.001,1000);
  // gPad->SetLogx();
  //gPad->SetLogy();
  ElectronRangeFunction->SetRange(1.,50.);
  ElectronRangeFunction->GetXaxis()->SetTitle("Electron energy (keV)");
  ElectronRangeFunction->GetYaxis()->SetTitle("Pratical range (mm)");
  ElectronRangeFunction->SetLineWidth(0);
 
  ElectronRangeFunction->Draw("ap");
  
  ElectronRangeCanvas->cd(2);//
  TGraph *ElasticFreePath = new TGraph(np, elenergy, ElasticFreePathMott);
  //gPad->SetLogx();//
  //gPad->SetLogy();//
  ElasticFreePath->GetXaxis()->SetTitle("Electron energy (keV)");
  ElasticFreePath->GetYaxis()->SetTitle("Elastic meanfreepath(#mum)");
  ElasticFreePath->GetXaxis()->SetRangeUser(1.,50.);
  ElasticFreePath->GetYaxis()->SetRangeUser(1.,100.);
  ElasticFreePath->SetTitle("Mott scattering");
  ElasticFreePath->Draw("ap");
  std::cout<<" CHECK"<<std::endl;
  ElectronRangeCanvas->cd(3);
  const Int_t nPoints = N_ENERGY_BINS;//N_POINTS;
  Double_t Energy[nPoints];
  Double_t PhotoelectricCrossSection[nPoints];
  Double_t Efficiency[nPoints];
  if (Verbose){
    std::cout << std::endl << "Efficiency data:" << std::endl;
    std::cout << "Energy(keV)  PhCrossSection(cm^2/g) " 
	      << " AbsLenght(cm)  Eff(" << GAP_THICKNESS << " cm gap)" << std::endl;
  }
  for (Int_t i=0; i<nPoints; i++){
    double ENERGY_MIN=1.;
    double ENERGY_MAX=50.;
    Energy[i] = ENERGY_MIN *pow(ENERGY_MAX/ENERGY_MIN,1.0*i/(nPoints-1));
    PhotoelectricCrossSection[i] = GetPhotoelectricCrossSection(Energy[i]);
    Efficiency[i] = GetEfficiency(Energy[i]);
  }
  TGraph *PhotoelectricCrossSectionGraph = new TGraph(nPoints, Energy, PhotoelectricCrossSection);
  PhotoelectricCrossSectionGraph->SetTitle("PhotoelectricCrossSection");
  PhotoelectricCrossSectionGraph->SetMarkerStyle(2);
  PhotoelectricCrossSectionGraph->SetMarkerSize(0.7);
  PhotoelectricCrossSectionGraph->SetMarkerColor(2);
  PhotoelectricCrossSectionGraph->GetXaxis()->SetTitle("Photon energy (keV)");
  PhotoelectricCrossSectionGraph->GetYaxis()->SetTitle("Photoelectric cross section (cm^2/g)");
  PhotoelectricCrossSectionGraph->GetYaxis()->SetTitleOffset(1.2);
  PhotoelectricCrossSectionGraph->Draw("ap");
  
  ElectronRangeCanvas->cd(4);
  
  TGraph *EfficiencyGraph = new TGraph(nPoints, Energy, Efficiency);
  EfficiencyGraph->SetTitle("Efficiency");
  EfficiencyGraph->SetMarkerStyle(2);
  EfficiencyGraph->SetMarkerSize(0.7);
  EfficiencyGraph->SetMarkerColor(2);
  EfficiencyGraph->GetXaxis()->SetTitle("Photon energy (keV)");
  EfficiencyGraph->GetYaxis()->SetTitle("Efficiency");
  EfficiencyGraph->GetYaxis()->SetRangeUser(0,1);
  EfficiencyGraph->GetYaxis()->SetTitleOffset(1.1);
  TString Name = GetName();
  TText *T = new TText(5,0.5,Name);
 
  EfficiencyGraph->Draw("ap");
  T->SetTextFont(12);
  T->Draw();

}


//*******************************************************************************
/// \brief Extract the element converting an incoming photon.                   *
/// The element in the mixture is extracted according to the relative values of *
/// the photoelectric cross sections at a given energy.                         *
/// Energy to be provided in keV.                                               *
//*******************************************************************************
TElement* TGasMixture::GetConvertingElement(Double_t ENERGY)
{
  Double_t Random = rnd->Uniform();
  Double_t Dummy  = 0.0;
  for (Int_t i=0; i<nElements; i++){
    Double_t AbsorptionProbability =
      (Elements[i]->GetPhotoelectricCrossSection(ENERGY)*Elements[i]->GetDensity())/
      (GetPhotoelectricCrossSection(ENERGY)*Density);
    if (Random < Dummy + AbsorptionProbability) return Elements[i];
    Dummy += AbsorptionProbability;
  }
  std::cout << "Warning: TGasMixture::GetConvertingElement." << std::endl;
  std::cout << "Exiting the element loop and returning last element of the mixture." << std::endl;
  std::cout << "You may consider revising the routine..." << std::endl;
  return Elements[nElements-1];
}


//*******************************************************************************
/// \brief Returns the the total elestic mean free path for the mixture.        *
/// The mean free path is evaluated by summing the inverse of the mean free     *
/// paths for all the elements.                                                 *
/// Energy to be provided in keV, mode cam be either RUTHERFORD or MOTT.        *
//*******************************************************************************
Double_t TGasMixture::GetElasticMeanFreePath(Double_t ENERGY, TString MODE)
{
  if (ENERGY < MIN_ANALYTIC_CROSS_SECTIONS_ENERGY){
    std::cout << "Warning: TGasMixture::GetElesticMeanFreePath." << std::endl;
    std::cout << "No data available for evaluating mean free path below "
	      << MIN_ANALYTIC_CROSS_SECTIONS_ENERGY << " eV."
	      << std::endl << "Returning -1." << std::endl;
    return -1;
  }
  Double_t MeanFreePath = 0.0;
  if ((MODE == "RUTHERFORD") || (MODE == "MOTT")){
    for (Int_t i=0; i<nElements; i++){
      MeanFreePath +=
	1.0/(Elements[i]->GetElasticMeanFreePath(ENERGY, Elements[i]->GetDensity(), MODE));
    }
    return 1.0/MeanFreePath;
  }
  else{
    std::cout << "Warning: TGasMixture::GetElasticMeanFreePath." << std::endl;
    std::cout << "Invalid mode in evaluating the elastic mean free path."
	      << std::endl << "Returning -1." << std::endl;
    return -1;
  }
}


//*******************************************************************************
/// \brief Extract the element scattering a travelling photoelectron.           *
/// The element in the mixture is extracted according to the relative values of *
/// the elastic cross sections at a given energy.                               *
/// Energy to be provided in keV, mode can be either RUTHERFORD OR MOTT.        *
//*******************************************************************************
TElement* TGasMixture::GetScatteringElement(Double_t ENERGY, TString MODE)
{
  Double_t Random = rnd->Uniform();
  Double_t Dummy  = 0.0;
  for (Int_t i=0; i<nElements; i++){
    Double_t ScatteringProbability = GetElasticMeanFreePath(ENERGY, MODE)/
      Elements[i]->GetElasticMeanFreePath(ENERGY, Elements[i]->GetDensity(), MODE);
    if (Random < Dummy + ScatteringProbability) return Elements[i];
    Dummy += ScatteringProbability;
  }
  std::cout << "Warning: TGasMixture::GetScatteringElement." << std::endl;
  std::cout << "Exiting the element loop and returning last element of the mixture." << std::endl;
  std::cout << "You may consider revising the routine..." << std::endl;
  return Elements[nElements-1];
}


//*******************************************************************************
/// \brief Returns the stopping power at a given energy.                        *
/// Energy to be provided in keV, outcome in keV/cm.                            *
//*******************************************************************************
Double_t TGasMixture::GetStoppingPower(Double_t ENERGY)
{
  if (ENERGY < MIN_ANALYTIC_CROSS_SECTIONS_ENERGY){
    std::cout << "Warning: TGasMixture::GetStoppingPower." << std::endl;
    std::cout << "No data available for evaluating stopping power below "
	      << MIN_ANALYTIC_CROSS_SECTIONS_ENERGY << " eV."
	      << std::endl << "Returning -1." << std::endl;
    return -1;
  }
  Double_t StoppingPower = 0.0;
  for (Int_t i=0; i<nElements; i++){
    StoppingPower += Elements[i]->GetStoppingPower(ENERGY);
  }
  return StoppingPower;

}
//validazione Photoelectriccrosssection per le mixtures
void  TGasMixture::PlotPhotoelectricCrossSection()
{
  const int N=30; 
  Double_t MPCSection[N];
  Double_t PE[N];
  for (int i=0; i<30; i++){
    PE[i]=1.0*(i+1);
    MPCSection[i]=GetPhotoelectricCrossSection(PE[i]);}
  TGraph *GrafMPCSection = new TGraph(N,PE,MPCSection);
  
  GrafMPCSection->Draw("al*");
  std::cout<<MPCSection[1]<<std::endl;
}
//validazione stoppingpower per le mixtures

