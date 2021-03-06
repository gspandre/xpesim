#include "TTrack.h"


//*******************************************************************************
/// \brief Basic constructor.                                                   *
//*******************************************************************************
TTrack::TTrack(TPhoton *PHOTON, TGasMixture *MIXTURE, TRandom *RND, Bool_t VERBOSE)
{
  //Dimension         = new TDimension();
  rnd               = RND;
  Photon            = PHOTON;
  Mixture           = MIXTURE;
  Verbose           = VERBOSE;
  nPrimaryElectrons = 0;
  ConvertingElement = Mixture->GetConvertingElement(Photon->GetEnergy());
  /*SetDimension (Dimension->GetGem_Radius()/2 ,
		Dimension->GetGem_Radius()/2,
		Dimension->GetZ_Drift(),
		Dimension->GetZ_Gem() );*/
  DistanceFromStart    = 0.0;
  MaxDistanceFromStart = 0.0;
  PhotoelectronPhi = 0.0; //cambio scope 
  PhotoelectronTheta = 0.0; //cambio scope 

  //dimensioni in cm
  /* xdim=1;
     ydim=1;
     zdim=1;
     zmin=0;*/
  //delete Dimension;//???????
}

void TTrack::SetDimension(Double_t a,Double_t b,Double_t c,Double_t d)
{
  xdim=a;
  ydim=b;
  zdim=c;
  zmin=d;
}
//propagate the photoelectron  
void TTrack::PropagatePhotoelectron()
{ 
  nPhotoelectronElasticScatterings = 0;
  PhotoelectronTrueRange           = 0.0;
  MaxDistanceFromStart             = 0.0;
  nPrimaryElectrons =0;
  ResidualEnergy  =  Photon->GetEnergy() - ConvertingElement->GetkEdge();
  //
  
  TString Elementino = ConvertingElement->GetChemicalSymbol();
  if( Elementino == "H")
  {std::cout<<" eccoci "<<Elementino<<std::endl;
  }
  //
  if(ResidualEnergy<=0) 
    std::cout<<"Stai caando fori dal vaso!! ResidualEnergy = "<<ResidualEnergy<<std::endl;
  ConversionPoint =  Photon->GetConversionPoint();
  nPairspath.push_back(0);//per far coincidere i vettori pairs e Photoelectronscattering
  //PrimaryIonizationV.push_back(TXYZ(0,0,0)); ELIMINARE OUT
  PhotoelectronScatteringV.push_back(ConversionPoint);
  PhotoelectronEnergy.push_back(ResidualEnergy);
  PhotoelectronTheta = Photon->GetPhotoelectronTheta();//cambio scope
  PhotoelectronPhi   = Photon->GetPhotoelectronPhi();//
  CX = sin(PhotoelectronTheta)*cos(PhotoelectronPhi);
  CY = sin(PhotoelectronTheta)*sin(PhotoelectronPhi);
  CZ = cos(PhotoelectronTheta);
  
  if (Verbose){
    std::cout << std::endl << "Photoelectron path (elastic scattrings coordinates):" << std::endl;
    std::cout << "X (cm), Y (cm), Z (cm), ResidualEnergy (keV)" << std::endl;
  }
  
  
  while (ResidualEnergy > MIN_ANALYTIC_CROSS_SECTIONS_ENERGY && 
	 PhotoelectronScatteringV[nPhotoelectronElasticScatterings].X()<xdim && 
	 PhotoelectronScatteringV[nPhotoelectronElasticScatterings].X()>-xdim &&
	 PhotoelectronScatteringV[nPhotoelectronElasticScatterings].Y()<ydim && 
	 PhotoelectronScatteringV[nPhotoelectronElasticScatterings].Y()>-ydim &&
	 PhotoelectronScatteringV[nPhotoelectronElasticScatterings].Z()<zdim &&
	 PhotoelectronScatteringV[nPhotoelectronElasticScatterings].Z()>zmin)
    
    {
      EvaluateNextStep("PHOTOELECTRON");
      TXYZ  ScatteringPoint = PhotoelectronScatteringV[nPhotoelectronElasticScatterings-1];
      
      DistanceFromStart =
	sqrt(pow(ScatteringPoint.X() - ConversionPoint.X(), 2.0) +
	     pow((ScatteringPoint.Y() - ConversionPoint.Y()), 2.0) +
	     pow((ScatteringPoint.Z() - ConversionPoint.Z()), 2.0));
      

      if (DistanceFromStart > MaxDistanceFromStart){
	MaxDistanceFromStart = DistanceFromStart;
	PhotoelectronPraticalRange = MaxDistanceFromStart;
      }
    }
  // Last electrons are created in the coordinates of the last collision.
  if(ResidualEnergy <= MIN_ANALYTIC_CROSS_SECTIONS_ENERGY && PrimaryIonizationV.size()>0)
    {
      Int_t ResidualPhotoelectronPairs = GetnPairs(ResidualEnergy);
      // std::cout<<" SIZE =  "<< PrimaryIonizationV.size()<<" "<< ResidualPhotoelectronPairs<<" "<<nPrimaryElectrons<<std::endl;
      for (Int_t i=0; i<ResidualPhotoelectronPairs; i++)
	{
	  PrimaryIonizationV.push_back(PrimaryIonizationV[nPrimaryElectrons-1]);
	}
      nPrimaryElectrons += ResidualPhotoelectronPairs;
    }
  PhotoelectronStartToEndRange = DistanceFromStart;
  
  if (Verbose){
    std::cout << "Number of collisions: " << nPhotoelectronElasticScatterings << std::endl;
    std::cout << "True range:           " << PhotoelectronTrueRange << " cm" << std::endl;
    std::cout << "Pratical range:       " << PhotoelectronPraticalRange << " cm" << std::endl;
    std::cout << "Start-to-end range:   " << PhotoelectronStartToEndRange << " cm" << std::endl; 
  }
  
  // Propagate the Auger electron (if any).
  nAugerElasticScatterings  = 0;
  AugerTrueRange            = 0.0;
  MaxDistanceFromStart      = 0.0;
  AugerCheck = 0;////
  if (rnd->Uniform() > ConvertingElement->GetFluorescenceYield()){
    //if (0.0> ConvertingElement->GetFluorescenceYield()){
    ResidualEnergy = ConvertingElement->GetkEdge();
    // std::cout<<""<<ResidualEnergy<<std::endl;
    //ResidualEnergy = 0.5;
    AugerCheck = 1;/////////115
    AugerScatteringV.push_back(ConversionPoint);
    AugerEnergy.push_back(ResidualEnergy);
    
    Double_t AugerTheta = Photon->GetAugerElectronTheta();
    Double_t AugerPhi   = Photon->GetAugerElectronPhi();
    CX = sin(AugerTheta)*cos(AugerPhi);
    CY = sin(AugerTheta)*sin(AugerPhi);
    CZ = cos(AugerTheta);
    if (Verbose){
      std::cout << std::endl << "Auger path (elastic scattrings coordinates):" << std::endl;
      std::cout << "X (cm), Y (cm), Z (cm), ResidualEnergy (keV)" << std::endl;
    }
    while (
	   ResidualEnergy > MIN_ANALYTIC_CROSS_SECTIONS_ENERGY &&
	   PhotoelectronScatteringV[nPhotoelectronElasticScatterings].X()<xdim && 
	   PhotoelectronScatteringV[nPhotoelectronElasticScatterings].X()>-xdim &&
	   PhotoelectronScatteringV[nPhotoelectronElasticScatterings].Y()<ydim && 
	   PhotoelectronScatteringV[nPhotoelectronElasticScatterings].Y()>-ydim &&
	   PhotoelectronScatteringV[nPhotoelectronElasticScatterings].Z()<zdim &&
	   PhotoelectronScatteringV[nPhotoelectronElasticScatterings].Z()>zmin)
      {
	EvaluateNextStep("AUGER");
	TXYZ  ScatteringPoint = AugerScatteringV[nAugerElasticScatterings-1];
	DistanceFromStart =
	  sqrt(pow(ScatteringPoint.X() - ConversionPoint.X(), 2.0) +
	       pow((ScatteringPoint.Y() - ConversionPoint.Y()), 2.0) +
	       pow((ScatteringPoint.Z() - ConversionPoint.Z()), 2.0));
	if (DistanceFromStart > MaxDistanceFromStart){
	  MaxDistanceFromStart = DistanceFromStart;
	  AugerPraticalRange = MaxDistanceFromStart;
	}
      }
    // Last electrons are created in the coordinates of the last collision.
    if(ResidualEnergy <= MIN_ANALYTIC_CROSS_SECTIONS_ENERGY && PrimaryIonizationV.size()>0)   
      { Int_t ResidualAugerPairs = GetnPairs(ResidualEnergy);
      for (Int_t i=0; i<ResidualAugerPairs; i++)
	{
	  PrimaryIonizationV.push_back(PrimaryIonizationV[nPrimaryElectrons-1]);   
	}
      nPrimaryElectrons += ResidualAugerPairs;
      }
    AugerStartToEndRange = DistanceFromStart;
    if (Verbose){
      std::cout << "Number of collisions: " << nAugerElasticScatterings << std::endl;
      std::cout << "True range:           " << AugerTrueRange << " cm" << std::endl;
      std::cout << "Pratical range:       " << AugerPraticalRange << " cm" << std::endl;
      std::cout << "Start-to-end range:   " << AugerStartToEndRange << " cm" << std::endl; 
    }
  }
  if (Verbose){
    std::cout << std::endl << "Number of primary electrons: " << nPrimaryElectrons
	      << " (" << 1000*Photon->GetEnergy()/nPrimaryElectrons
	      << " eV per pair)" << std::endl;
  }
}


//*******************************************************************************
/// \brief Evaluate the coordinates of the next step in photoelectron path.     *
/// The formula are taken from Joy's book - note that there is an error in the  *
/// book: V1=AM*sin(Phi) -> V1=AN*sin(Phi).                                     *
/// The routine also takes care of generating the pimary ionization.            * 
//*******************************************************************************
void TTrack::EvaluateNextStep(TString MODE)
{
  TXYZ ScatteringPoint;
  if (Verbose){
    if (MODE == "PHOTOELECTRON"){
      ScatteringPoint = PhotoelectronScatteringV[nPhotoelectronElasticScatterings];
      std::cout << ScatteringPoint.X() << "  "
		<< ScatteringPoint.Y() << "  "
		<< ScatteringPoint.Z() << "  "
		<< ResidualEnergy << std::endl;
    }
    else if (MODE == "AUGER"){
      ScatteringPoint = AugerScatteringV[nAugerElasticScatterings];
      
      std::cout << ScatteringPoint.X() << "  "
		<< ScatteringPoint.Y() << "  "
		<< ScatteringPoint.Z() << "  "
		<< ResidualEnergy << std::endl;
    }
  }
  Double_t Path =rnd->Exp(Mixture->GetElasticMeanFreePath(ResidualEnergy, "MOTT"));
  
  
  Double_t Phi  =
    Mixture->GetScatteringElement(ResidualEnergy, "MOTT")->GetScatteringAngle(ResidualEnergy, "MOTT");
  
  
  Double_t Psi  = rnd->Uniform(0, 2*kPI);
  Double_t V1,V2;
  if(CZ == 0)
    {
      V1   = 0;
      V2   = sin(Phi);
    }
  else
    {
      Double_t AN   = -CX/CZ;
      Double_t AM   = 1.0/sqrt(1 + AN*AN);
      V1   = AM*sin(Phi);
      V2   = AN*AM*sin(Phi);
    }
  Double_t V3   = cos(Psi);
  Double_t V4   = sin(Psi);
  Double_t CA   = CX*cos(Phi) + V1*V3 + CY*V2*V4;
  Double_t CB   = CY*cos(Phi) + V4*(CZ*V1 - CX*V2);
  Double_t CC   = CZ*cos(Phi) + V2*V3 - CY*V1*V4;
  CX = CA;
  CY = CB;
  CZ = CC;
  Double_t EnergyLoss = Path*Mixture->GetStoppingPower(ResidualEnergy); 
  ResidualEnergy -= EnergyLoss;
  //  if (ResidualEnergy< MIN_ANALYTIC_CROSS_SECTIONS_ENERGY){std::cout<<MODE<<" "<<ResidualEnergy<<"aleoo "<<PhotoelectronScatteringV.size()<<std::endl;
  
  if (MODE == "PHOTOELECTRON")
    {
      ScatteringPoint = PhotoelectronScatteringV[nPhotoelectronElasticScatterings];
      double x = ScatteringPoint.X();
      double y = ScatteringPoint.Y();
      double z = ScatteringPoint.Z();
      
      ScatteringPoint.SetXYZ(x + Path*CA, y + Path*CB, z + Path*CC);
      
      PhotoelectronScatteringV.push_back(ScatteringPoint);
      PhotoelectronEnergy.push_back(ResidualEnergy);
      
      Int_t nPairs = GetnPairs(EnergyLoss);
      nPrimaryElectrons += nPairs;
      nPairspath.push_back(nPairs);//coppie in un path
      for (Int_t i=0; i<nPairs; i++){
	Double_t Position = rnd->Uniform();
	TXYZ IonizationPoint = TXYZ(x + Path*CA*Position, y + Path*CB*Position, z + Path*CC*Position);
	
	PrimaryIonizationV.push_back(IonizationPoint);
      }
      nPhotoelectronElasticScatterings ++;
      PhotoelectronTrueRange += Path;
    }
  else if (MODE == "AUGER"){
    ScatteringPoint = AugerScatteringV[nAugerElasticScatterings];
    double x = ScatteringPoint.X();
    double y = ScatteringPoint.Y();
    double z = ScatteringPoint.Z();
    
    TXYZ ScatteringPoint = TXYZ(x + Path*CA, y + Path*CB, z + Path*CC);
    AugerScatteringV.push_back(ScatteringPoint);
    
    AugerEnergy.push_back(ResidualEnergy);
    Int_t nPairs = GetnPairs(EnergyLoss);
    nPrimaryElectrons += nPairs;
    for (Int_t i=0; i<nPairs; i++){
      Double_t Position = rnd->Uniform();
      TXYZ IonizationPoint = TXYZ(x + Path*CA*Position, y + Path*CB*Position, z + Path*CC*Position);
      PrimaryIonizationV.push_back(IonizationPoint);
    }
    nAugerElasticScatterings ++;
    AugerTrueRange += Path;
  }
}

//*******************************************************************************
/// \brief Returns the number of e-ion pairs generated between two collisions.  *
//*******************************************************************************
Int_t TTrack::GetnPairs(Double_t ENERGY_LOSS)
{
  //return (Int_t)(0.54 + rnd->Exp(1000*ENERGY_LOSS/Mixture->GetWIonization()));
  // new version Compound Poisson distribution
  Double_t MeanNumberSecondary =  1000.*ENERGY_LOSS/(Mixture->GetWIonization());
  Int_t NumberPrimary = rnd->Poisson(MeanNumberSecondary/3.);
  Int_t NumberSecondary = 0;
  for (int i=0; i<NumberPrimary; i++)
    {
      NumberSecondary +=rnd->Poisson(3.);
    }
  return NumberSecondary;
}


//*******************************************************************************
/// \brief Plots the path of the photoelectron and of the Augerelectron...      *
/// ...according to the positions of the elestic scatterings.                   *
//*******************************************************************************
void TTrack::PlotPath()
{ //if (PhotoelectronScatteringV.size()<2){PropagatePhotoelectron();}
  const Int_t NP = nPhotoelectronElasticScatterings;
  Double_t XP[NP];
  Double_t YP[NP];
  Double_t ZP[NP];
  for (Int_t i=0; i<NP; i++){
    TXYZ PhotoelectronPoint =  PhotoelectronScatteringV[i];
    XP[i] = PhotoelectronPoint.X()*10000.0;
    YP[i] = PhotoelectronPoint.Y()*10000.0;
    ZP[i] = PhotoelectronPoint.Z()*10000.0;

  }
  const Int_t NA = nAugerElasticScatterings;
  Double_t XA[NA];
  Double_t YA[NA];
  Double_t ZA[NA];
  for (Int_t i=0; i<NA; i++){
    TXYZ AugerPoint =  AugerScatteringV[i];
    XA[i] = AugerPoint.X()*10000.0;
    YA[i] = AugerPoint.Y()*10000.0;
    ZA[i] = AugerPoint.Z()*10000.0;
    
  }
  Double_t Dimension = 0.01*(Int_t)(PhotoelectronTrueRange*30*10000);
  TCanvas *TrackPathCanvas = new TCanvas("TrackPath", "TrackPath", 50, 50, 900, 800);
  TrackPathCanvas->SetFillColor(10);
  TrackPathCanvas->Divide(2, 2);
  TrackPathCanvas->cd(1);
  TGraph *XYPGraph = new TGraph(NP, XP, YP);
  XYPGraph->SetTitle("XYProjection");
  XYPGraph->GetHistogram()->SetBins(100, ConversionPoint.X()-Dimension, ConversionPoint.X()+Dimension);
  XYPGraph->SetMinimum(ConversionPoint.Y()-Dimension);
  XYPGraph->SetMaximum(ConversionPoint.Y()+Dimension);
  // XYPGraph->GetHistogram()->SetBins(100, -5000,5000 );
  //XYPGraph->SetMaximum(5000);
  // XYPGraph->SetMinimum(-5000);
  XYPGraph->GetXaxis()->SetTitle("X coordinate (um)");
  XYPGraph->GetYaxis()->SetTitle("Y coordinate (um)");
  XYPGraph->GetYaxis()->SetTitleOffset(1.3);
  XYPGraph->SetLineColor(4);
  XYPGraph->SetLineWidth(2);
  XYPGraph->Draw("alp");
  TGraph *XYAGraph = new TGraph(NA, XA, YA);
  XYAGraph->SetLineColor(2);
  XYAGraph->SetLineWidth(2);
  XYAGraph->Draw("l");
  TrackPathCanvas->cd(2);
  TGraph *XZPGraph = new TGraph(NP, XP, ZP);
  XZPGraph->SetTitle("XZProjection");
  XZPGraph->GetHistogram()->SetBins(100, ConversionPoint.X()-Dimension, ConversionPoint.X()+Dimension);
  //XZPGraph->SetMinimum(ConversionPoint.Z()-Dimension);
  //XZPGraph->SetMaximum(ConversionPoint.Z()+Dimension);
  XZPGraph->GetXaxis()->SetTitle("X coordinate (um)");
  XZPGraph->GetYaxis()->SetTitle("Z coordinate (um)");
  XZPGraph->GetYaxis()->SetTitleOffset(1.3);
  XZPGraph->SetLineColor(4);
  XZPGraph->SetLineWidth(2);
  XZPGraph->Draw("alp");
  TGraph *XZAGraph = new TGraph(NA, XA, ZA);
  XZAGraph->SetLineColor(2);
  XZAGraph->SetLineWidth(2);
  XZAGraph->Draw("l");
  TrackPathCanvas->cd(3);
  TGraph *YZPGraph = new TGraph(NP, YP, ZP);
  YZPGraph->SetTitle("YZProjection");
  //  YZPGraph->GetHistogram()->SetBins(100, ConversionPoint.Y()-Dimension,ConversionPoint.Y()+ Dimension);
  //  YZPGraph->SetMinimum(ConversionPoint.Z()-Dimension);
  //  YZPGraph->SetMaximum(ConversionPoint.Z()+Dimension);
  YZPGraph->GetXaxis()->SetTitle("Y coordinate (um)");
  YZPGraph->GetYaxis()->SetTitle("Z coordinate (um)");
  YZPGraph->GetYaxis()->SetTitleOffset(1.3);
  YZPGraph->SetLineColor(4);
  YZPGraph->SetLineWidth(2);
  YZPGraph->Draw("alp");
  TGraph *YZAGraph = new TGraph(NA, YA, ZA);
  YZAGraph->SetLineColor(2);
  YZAGraph->SetLineWidth(2);
  YZAGraph->Draw("l");
  TrackPathCanvas->cd(4);
}


//*******************************************************************************
/// \brief Plots the distribution of primary ionization.                        *
//*******************************************************************************
void TTrack::PlotPrimaryIonization()
{ if (PhotoelectronScatteringV.size()<2){PropagatePhotoelectron();}
  const Int_t N = nPrimaryElectrons;
  Double_t X[N];
  Double_t Y[N];
  Double_t Z[N];
  for (Int_t i=0; i<N; i++){
    TXYZ PrimaryIonizationPoint =  PrimaryIonizationV[i];
    X[i] = PrimaryIonizationPoint.X()*10000.0; //cm-> um
    Y[i] = PrimaryIonizationPoint.Y()*10000.0;
    Z[i] = PrimaryIonizationPoint.Z()*10000.0;
    /*
      X[i] = PrimaryIonizationX[i]*10000.0;
      Y[i] = PrimaryIonizationY[i]*10000.0;
      Z[i] = PrimaryIonizationZ[i]*10000.0; 
    */
  }
  Double_t Dime = 0.01*(Int_t)(PhotoelectronTrueRange*30*10000);
  TCanvas *PrimaryIonizationCanvas =
    new TCanvas("PrimaryIonization", "PrimaryIonization", 60, 60, 910, 910);
  PrimaryIonizationCanvas->SetFillColor(10);
  PrimaryIonizationCanvas->Divide(2, 2);
  PrimaryIonizationCanvas->cd(1);
  TGraph *XYGraph = new TGraph(N, X, Y);
  XYGraph->SetTitle("XYProjection");
  XYGraph->GetHistogram()->SetBins(100, ConversionPoint.X()-Dime,ConversionPoint.X()+ Dime);
  XYGraph->SetMinimum(ConversionPoint.Y()-Dime);
  XYGraph->SetMaximum(ConversionPoint.Y()+Dime);
  XYGraph->GetXaxis()->SetTitle("X coordinate (um)");
  XYGraph->GetYaxis()->SetTitle("Y coordinate (um)");
  XYGraph->GetYaxis()->SetTitleOffset(1.3);
  XYGraph->SetLineColor(4);
  XYGraph->SetLineWidth(2);
  XYGraph->SetMarkerStyle(26);
  XYGraph->SetMarkerSize(0.15);
  XYGraph->Draw("ap");
  PrimaryIonizationCanvas->cd(2);
  TGraph *XZGraph = new TGraph(N, X, Z);
  XZGraph->SetTitle("XZProjection");
  XZGraph->GetHistogram()->SetBins(100, ConversionPoint.X()-Dime, ConversionPoint.X()+Dime);
  XZGraph->SetMinimum(ConversionPoint.Z()*10000-Dime);
  XZGraph->SetMaximum(ConversionPoint.Z()*10000+Dime);
  XZGraph->GetXaxis()->SetTitle("X coordinate (um)");
  XZGraph->GetYaxis()->SetTitle("Z coordinate (um)");
  XZGraph->GetYaxis()->SetTitleOffset(1.3);
  XZGraph->SetLineColor(4);
  XZGraph->SetLineWidth(2);
  XZGraph->SetMarkerStyle(26);
  XZGraph->SetMarkerSize(0.15);
  XZGraph->Draw("ap");
  PrimaryIonizationCanvas->cd(3);
  TGraph *YZGraph = new TGraph(N, Y, Z);
  YZGraph->SetTitle("YZProjection");
  YZGraph->GetHistogram()->SetBins(100, ConversionPoint.Y()-Dime,ConversionPoint.Y() + Dime);
  YZGraph->SetMinimum(ConversionPoint.Z()*10000-Dime);
  YZGraph->SetMaximum(ConversionPoint.Z()*10000+Dime);
  YZGraph->GetXaxis()->SetTitle("Y coordinate (um)");
  YZGraph->GetYaxis()->SetTitle("Z coordinate (um)");
  YZGraph->GetYaxis()->SetTitleOffset(1.3);
  YZGraph->SetLineColor(4);
  YZGraph->SetLineWidth(2);
  YZGraph->SetMarkerStyle(26);
  YZGraph->SetMarkerSize(0.15);
  YZGraph->Draw("ap");
  PrimaryIonizationCanvas->cd(4);
}

/*void TTrack::Drift(double Z_MIN)
  {
  double DiffusionSigma = Mixture->GetDiffusionSigma(); //mu/sqrt(cm) 
  //  std::cout<<"DiffusionSigma = "<<DiffusionSigma<<std::endl;
  for(int i = 0; i < nPrimaryElectrons;i++)
  {
  TXYZ ElectronPosition = PrimaryIonizationV[i];
  
  if (ElectronPosition.Z() > Z_MIN)
	{
	double sigma = DiffusionSigma*sqrt(ElectronPosition.Z()-Z_MIN)/1e4; // cm
	//	  std::cout<<ElectronPosition.X()<<" "<<ElectronPosition.Y()<<" "<<ElectronPosition.Z()<<" "<<sigma<<std::endl;
	
	ElectronPosition.SetXYZ(ElectronPosition.X()+rnd->Gaus(0,sigma),
	ElectronPosition.Y()+rnd->Gaus(0,sigma),
	Z_MIN);
	
	PrimaryIonizationV[i] = ElectronPosition;
	}
	}
	}  */

void TTrack::Drift()
{
  double DiffusionSigma = Mixture->GetDiffusionSigma();
  double Xelectron,Yelectron,Zelectron,sigma;
  std::vector<TXYZ>::iterator pos;
  for(pos =  PrimaryIonizationV.begin(); pos !=  PrimaryIonizationV.end();++pos)
    {
      Xelectron=(*pos).X();
      Yelectron=(*pos).Y();
      Zelectron=(*pos).Z();
      if(Zelectron>zmin) 
      	{
	  
	  sigma = DiffusionSigma*sqrt(Zelectron-zmin)/10000;
	  Xelectron+=rnd->Gaus(0,sigma);
	  Yelectron+=rnd->Gaus(0,sigma);
	  DiffElectronPosition.push_back(std::make_pair( Xelectron, Yelectron));
	  
	}
      /* else
	{std::cout<<"Under Gem"<<std::endl;
	//std::cout<<"Diff electron size"<< DiffElectronPosition.size()<<std::endl;
	}*/
    }
}




