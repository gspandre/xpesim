#include <iostream>
#include <fstream>
#include "TRandomGenerator.h"
#include "TRandom.h"
#include "TGasMixture.h"
#include "TDetector.h"
#include "TTree.h"
#include "TFile.h"
#include <vector>
#include "TReadout.h"
#include "TH2D.h"
#include "TDimension.h"
#include "TGem.h"
#include "TView.h"
#include "THexagon3d.h"
#include "TCluster.h"
#include "TPaveText.h"
#include "TText.h"
#include "TMDP.h"
#include "TTreeAnalysis.h"
#include "TVarPolFlux.h"
#include "TExperiment.h"
/////
// Default macro parameters.
Int_t    MixID                = 1;//
Double_t PhotonEnergy         = 6;//!!!!!!!5.9
Double_t PolarizationAngle    = 90.0;//GRADI!0.0;
Double_t PolarizationDegree   = 0.0 ;//1.0

TH1D *Spectrum;
// Init relevant variables.
TRandomGenerator *RandomGenerator = new TRandomGenerator();
//RandomGenerator->SetRandomSeed();
TRandom *rnd = RandomGenerator->GetRandomGenerator();
TDimension *Dimension=new TDimension();

TGasMixture *Mixture = new TGasMixture(MixID, rnd,Dimension,1);
TSource *Source = new TSource(PhotonEnergy, PolarizationAngle, PolarizationDegree, rnd, 0);
//aggiungo conversione dipendente dalla sezione d'urto fotoelettrica
Double_t Zconversion=1.5;
////////inizio funzioni Her

//TMDP FHER;
//TF1 *HerD = FHER.HerculesSpectrum();

///////fine Funzioni Her
TXYZ ConversionPoint(0.0,0.0,Zconversion);

void SetSpectrum(TH1D *h)
{
  Spectrum = h;
}

void SetDimensionMix(double Pressure, double Thickness, double Gem_Pitch,double Pitch)
{
  Dimension->SetPressure(Pressure);
  Dimension->SetZ_Drift(Thickness + Dimension->GetZ_Gem());
  Dimension-> SetGem_Pitch(Gem_Pitch);
  Dimension->SetPitch(Pitch);
  //
  if(Pitch==80)
    Dimension->SetPixel(160,138);
  else Dimension->SetPixel(352,300);//PER GLORIA
  Mixture = new TGasMixture(MixID, rnd,Dimension,1);
}

void  ComputeConversionPoint(double Energy = PhotonEnergy)
{
  double Lambda=1./(Mixture->GetPhotoelectricCrossSection(Energy)*(Mixture->GetDensity()));
  double Zconversion=Dimension->GetZ_Drift()-rnd->Exp(Lambda);
  ConversionPoint.SetXYZ(0.0,0.0,Zconversion);
}

//Calcola la probabilita' che il rame sul GEM riemetta un fotone
double ProbabilityCuFluorescence(double Energy = PhotonEnergy)
{
  double Density = 8.920;//gr/cm^3
  const int N = 5;
  double energy[N] = {8.980,10,15,20,30};
  double sigmaphoto[N] = {277, 214, 73.1, 33.1, 10.5};
  TGraph *GCu = new TGraph(N,energy,sigmaphoto);
  double Sigma = GCu->Eval(Energy);
  double Lambda = 1./(Sigma*Density);//cm
  double AbsorptionProb = 1 - exp(- 0.0004/Lambda );
  AbsorptionProb *= 0.5;//buchi nel GEM diminuiscono S efficace;
  //posizione emissione fotone
  AbsorptionProb *= 0.407;//Fluorescence Cu;
  double count=0;
  for (int i=0; i<100000; i++)
    {
      double X = rnd->Uniform(-0.5,0.5);
      double Y = rnd->Uniform(-0.5,0.5);
      //direzione emissione
      double theta = acos(rnd->Uniform(-1,1));
      double phi = rnd->Uniform(-TMath::Pi(),TMath::Pi());
      //distanza conversione
      double LambdaFluo = 1./(Mixture->GetPhotoelectricCrossSection(8)*(Mixture->GetDensity()));
      double ConversionDistance=rnd->Exp(LambdaFluo);
      X += ConversionDistance*sin(theta)*cos(phi);
      Y += ConversionDistance*sin(theta)*sin(phi);
      double Z = 0.6 +  ConversionDistance * cos(theta);//Dimension->GetZ_Drift()
      if (X>-0.5 && X<0.5 && Y<0.5 && Y>-0.5 && Z<Dimension->GetZ_Drift() && Z>0.6)
	{
	  //ConversionPoint.SetXYZ(X,Y,Z);
	  count++;
	}
    }
  std::cout<<"Riemission prob (%)= "<<AbsorptionProb*100 <<std::endl;
  std::cout<<"% absorbed  = "<<count/1000.<<std::endl;
  AbsorptionProb*=count/100000.;
  return AbsorptionProb; 
} 

void ConversionPattern(double Energy) // NOT USED NdG 13/12/06
{
  double X = 0.;
  double Y = 0.;
  double nHole = rnd->Uniform(0.,1.);//(0.,9.);
 
  if (nHole<0.49999)
   { 
     X = -0.2;
     Y = 0.;
   }
  else 
    {
      X = 0.;
      Y = 0.;
    }
  double H = 0.4;//0.4
  double Radius = 0.04;//0.3
  double PosTheta = 0.;//Theta direzione fotone
  double PosRad = 0.;//raggio
  double PosPhi = 0.;//Phi direzione fotone
  double PosAngle = 0.;//coordinata angolare del punto di emissione  
  double Conv = 1.;
  double ConvMax = 0.1;
  double Zconversion = 1.;
  while ( pow(PosRad*cos(PosAngle)+ H*sin(PosTheta)*cos(PosPhi),2) + pow(PosRad*sin(PosAngle)+ H*sin(PosTheta)*sin(PosPhi),2)>  (Radius* Radius)    ||   
	  Conv > ConvMax || 
	  Zconversion < 0.6)   
    {// esce dal loop quando il punto di conversione e' nello spessore di assorb. (>0.6), in un raggio di 40 micron, e conv<Convmax
      Energy = rnd->Gaus(5.89,0.1);//Gaus(5.9,0.2);
      double PEnergy=0;
      double RndS = rnd->Uniform();
         //Fe
      if (RndS<0.96)//0.96 Ne
	{
	  PEnergy =rnd->Gaus(5.89,0.1);//Gaus(5.9,0.2);
	}
      else if (RndS<0.98) PEnergy=rnd->Gaus(4.1,0.3);//995 Ne
      else PEnergy = rnd->Gaus(5,0.3);
      
      // 5.4
      //PEnergy =rnd->Gaus(5.4,0.1);
      //Energy = PEnergy;

      PhotonEnergy = Energy;
      PosTheta = rnd->Uniform(0.,2.);
      PosPhi = rnd->Uniform(0.,2.* TMath::Pi());
      PosRad = rnd->Uniform(0., Radius);
      PosAngle = rnd->Uniform(0.,2.*TMath::Pi());
      ComputeConversionPoint(Energy);
      Zconversion = ConversionPoint.Z();
      Conv = -Zconversion + Dimension->GetZ_Drift();
      ConvMax = (Dimension->GetZ_Drift()-0.6)/cos(PosTheta);
    }
  
  //Zconversion e' calcolato in S di Rif con 0 sul readout
  double Ztot = Dimension->GetZ_Drift() - Zconversion + H;
  X+=PosRad*cos(PosAngle)+( Ztot )*sin(PosTheta)*cos(PosPhi);
  Y+=PosRad*sin(PosAngle)+( Ztot )*sin(PosTheta)*sin(PosPhi);
  // std::cout<< "  X = "<<X<<"  Y = "<<Y<<"Pos angle"<<PosAngle<<std::endl;
  //
  /*double NRad = Ztot * sin(atan(0.3/4.)) + 0.03;
    double PosRad = rnd->Uniform(0.,NRad); 
    double PosPhi = rnd->Uniform(0.,2.* TMath::Pi());
    X += PosRad*cos(PosPhi);
    Y += PosRad*sin(PosPhi);*/ 
  ConversionPoint.SetXYZ(X ,Y , Zconversion); 
  
}

void SetMixtureID(Int_t ID)
{
  MixID = ID;
  Mixture = new TGasMixture(MixID, rnd,Dimension,1);
}

void TestMixture()
{
  //Mixture->PlotEfficiency();
  Mixture->PlotElectronRange();
}

void SetPhotonEnergy(Double_t ENERGY)
{
  Source->SetEnergy(ENERGY);
}


void SetConversionPoint(Double_t x, Double_t y, Double_t z)
{
  ConversionPoint.SetXYZ(x,y,z);
}


void SetPolarizationAngle(Double_t ANGLE)
{
  Source->SetPolarizationAngle(ANGLE);
}


void SetPolarizationDegree(Double_t DEGREE)
{
  Source->SetPolarizationDegree(DEGREE);
}

void TestAngularDistributions(Int_t nTracks=3000)
{
  TCanvas *TestAnglesCanvas = new TCanvas("TestAngularDistributions",
					  "TestAngularDistributions", 50, 50, 900, 900);
  
  const int N = nTracks;
  double XAuger[N];
  double YAuger[N];
  double ZAuger[N];
  TestAnglesCanvas->SetFillColor(10);
  TH1D *PhotoelectronThetaHistogram = new TH1D("PhotoelectronTheta", "PhotoelectronTheta", 60, 0, kPI);
  TH1D *PhotoelectronPhiHistogram = new TH1D("PhotoelectronPhi", "PhotoelectronPhi", 60, 0, 2*kPI);
  TH1D *AugerThetaHistogram = new TH1D("AugerTheta", "AugerTheta", 60, -1, 1);
  TH1D *AugerPhiHistogram = new TH1D("AugerPhi", "AugerPhi", 60, 0, 2*kPI);
  for (Int_t i=0; i<nTracks; i++){
    if (!((i+1)%100) && i != 0) std::cout << i+1 << " photons generated." << std::endl;
    TPhoton *Photon = new TPhoton(Source, ConversionPoint, rnd, 0);
    PhotoelectronThetaHistogram->Fill(Photon->GetPhotoelectronTheta());
    PhotoelectronPhiHistogram->Fill(Photon->GetPhotoelectronPhi());
    AugerThetaHistogram->Fill(cos(Photon->GetAugerElectronTheta()));/////
    AugerPhiHistogram->Fill(Photon->GetAugerElectronPhi());
    double Phi = Photon->GetAugerElectronPhi(); 
    double Theta = Photon->GetAugerElectronTheta();
    XAuger[i] = sin(Theta)*cos(Phi);
    YAuger[i] = sin(Theta)*sin(Phi);
    ZAuger[i] = cos(Theta);
    delete Photon;
  }
  TestAnglesCanvas->Divide(2, 2);
  TestAnglesCanvas->cd(1);
  PhotoelectronThetaHistogram->Draw();
  // PhotoelectronThetaHistogram->Fit(Source->GetThetaDistribution(), "", "", 0.0, kPI);
  TestAnglesCanvas->cd(2);
  PhotoelectronPhiHistogram->Draw();
  //PhotoelectronPhiHistogram->Fit(Source->GetPhiDistribution(), "", "", 0.0, 2*kPI);
  TestAnglesCanvas->cd(3);
  AugerThetaHistogram->Draw();
  //AugerThetaHistogram->Fit("pol0", "", "", 0.0, kPI);
  TestAnglesCanvas->cd(4);
  AugerPhiHistogram->Draw();
  //AugerPhiHistogram->Fit("pol0", "", "", 0.0, 2*kPI);
  TCanvas *c2 = new TCanvas("c2","c2",1200,800);
  c2->Divide(2,2);
  TGraph *g1 = new TGraph(N, XAuger, YAuger);
  c2->cd(1);
  g1->Draw("ap");
  c2->cd(2);
  TGraph *g2 = new TGraph(N, ZAuger, XAuger);
  g2->Draw("ap");
  TGraph *g3 = new TGraph(N, ZAuger, YAuger);
  c2->cd(3);
  g3->Draw("ap");
  c2->cd(4);

}


void TestPrimaryIonization(Int_t nTracks=1000)
{
  TCanvas *TestIonizationCanvas = new TCanvas("TestIonization", "TestIonization", 50, 50, 600, 500);
  TestIonizationCanvas->SetFillColor(10);
  Double_t nPairsAverage = PhotonEnergy*1000.0/Mixture->GetWIonization();
  Int_t nPairsMin = (Int_t)(nPairsAverage*0.5);
  Int_t nPairsMax = (Int_t)(nPairsAverage*1.5);
  Int_t nBins     = (nPairsMax - nPairsMin)/3;
  TH1D *nPairsHistogram = new TH1D("PrimaryElectrons", "PrimaryElectrons",
				   nBins, nPairsMin, nPairsMax);
  
  for (Int_t i=0; i<nTracks; i++){
    if (!((i+1)%100) && i != 0) std::cout << i+1 << " tracks generated." << std::endl;
    TPhoton *Photon = new TPhoton(Source, ConversionPoint, rnd, 0);
    TTrack *Track = new TTrack(Photon, Mixture, rnd, 0);
    Track->PropagatePhotoelectron(); 
    nPairsHistogram->Fill(Track->GetnPrimaryElectrons());
  }
  nPairsHistogram->Draw();
}

void PlotTrack()
{
  std::cout<<ConversionPoint.Z()<<std::endl;
  TPhoton *Photon = new TPhoton(Source, ConversionPoint, rnd, 0);
  TH1D *TRangeH = new TH1D("TrueRange","TrueRange",1000,0,0.5);
  for (int i = 0; i<1000; i++)
    {
      TTrack Track(Photon, Mixture, rnd, 0);
      // TTrack *Track   = new TTrack(Photon, Mixture, rnd, 0);
      Track.SetDimension (Dimension->GetGem_Radius()/2 ,
			  Dimension->GetGem_Radius()/2,
			  Dimension->GetZ_Drift(),
			  Dimension->GetZ_Gem() );
      Track.PropagatePhotoelectron(); 
      //Track->PlotPath();
      //Track->PlotPrimaryIonization();
      double TrueRange =  Track.GetPhotoelectronTrueRange();
      //std::cout<<"TrueRange ="<<TrueRange<<std::endl;
      TRangeH->Fill(TrueRange);
    }
  TRangeH->Draw();
}


void DrawPath()
{
 TDetector myDetector(Mixture,rnd,0,Dimension);
 TPhoton Photon(Source, ConversionPoint, rnd, 0);
 TTrack Track(&Photon, Mixture, rnd, 0);
}

//crea un tree di dati montecarlo con tutte le informazioni sulle tracce e i primari fino alla GEM
void GetTree(int Nevents=100, TString FileName = "treetracks.root")
{
  std::cout<<"Photon enrgy  = "<< PhotonEnergy<<std::endl;
  TTree *tree=new TTree("Tree","datisim");
  Double_t xp[1000];
  Double_t yp[1000];
  Double_t zp[1000];
  Double_t xAuger[500];
  Double_t yAuger[500];
  Double_t zAuger[500];
  Int_t nAuger;
  Int_t nPa[1000];
  Double_t X[2000];
  Double_t Y[2000];
  Double_t Z[2000];
  Int_t NP;
  Int_t N;
  Int_t ntr;
  Int_t k;
  Int_t Nd;
  Double_t xe[2000];//dif
  Double_t ye[2000];//dif
  double Phi;
  double Theta;
  double ZC;
  double XI=0;
  double YI=0;
  tree->Branch("NP",&NP,"NP/I");
  tree->Branch("xPhotoel",xp,"x[NP]/D");
  tree->Branch("yPhotoel",yp,"y[NP]/D");
  tree->Branch("zPhotoel",zp,"z[NP]/D");
  tree->Branch("InitialPhi",&Phi,"Phi/D");
  tree->Branch("InitialTheta",&Theta,"Theta/D");
  tree->Branch("nA",&nAuger,"nAuger/I");
  tree->Branch("xAuger",xAuger,"xAuger[nAuger]/D");
  tree->Branch("yAuger",yAuger,"yAuger[nAuger]/D");
  tree->Branch("zAuger",zAuger,"zAuger[nAuger]/D");
  tree->Branch("N",&N,"N/I");
  tree->Branch("nPrimaryI",nPa,"nPa[NP]/I");
  tree->Branch("XPrimaryI",X,"X[N]/D");
  tree->Branch("YPrimaryI",Y,"Y[N]/D");
  tree->Branch("ZPrimaryI",Z,"Z[N]/D");
  tree->Branch("ntr",&ntr,"ntr/I");
  tree->Branch("nUndetectedPhoton",&k,"k/I");
  tree->Branch("Ndrift",&Nd,"Nd/I");
  tree->Branch("xPrimaryElectronPosition",xe,"xe[N]/D");
  tree->Branch("yPrimaryElectronPosition",ye,"ye[N]/D");
  tree->Branch("ConversionZ",&ZC,"ZC/D");
  tree->Branch("XInitial",&XI ,"XI/D");
  tree->Branch("YInitial",&YI ,"YI/D");
  for ( ntr=0; ntr<Nevents; ntr++ )
    { if (ntr %100 == 0)std::cout<<"--------------"<<ntr<<std::endl;
    ComputeConversionPoint();
    k=0;
    while (ConversionPoint.Z()<0.6)
      { 
	ComputeConversionPoint();
	k++;
      }
    XI = ConversionPoint.X();
    YI =ConversionPoint.Y();
    ZC = ConversionPoint.Z();
    TPhoton Photon(Source,ConversionPoint, rnd, 0);
    TTrack Track(&Photon, Mixture, rnd, 0);
    Track.SetDimension (Dimension->GetGem_Radius()/2 ,
		Dimension->GetGem_Radius()/2,
		Dimension->GetZ_Drift(),
		Dimension->GetZ_Gem() );
      Track.PropagatePhotoelectron(); 
      Phi = Track.GetPhotoelectronPhi();
      Theta = Track.GetPhotoelectronTheta();
      Track.Drift();
      
      NP = Track.GetPhotoelectronScatteringV().size();

      for (int i=0; i<NP; i++)
	{
	  TXYZ PhotoEPoint = Track.GetPhotoelectronScatteringV()[i];
	  xp[i] = PhotoEPoint.X()*10000;
	  yp[i] = PhotoEPoint.Y()*10000;
	  zp[i] = PhotoEPoint.Z()*10000;
	  nPa[i]=Track.GetnPairspath()[i]; 
	}
      
      nAuger=Track.GetAugerScatteringV().size();
      for  (int i=0; i<nAuger; i++)
	{
	  TXYZ AugerScatteringPoint=Track.GetAugerScatteringV()[i];
	  xAuger[i]=AugerScatteringPoint.X()*10000;
	  yAuger[i]=AugerScatteringPoint.Y()*10000;
	  zAuger[i]=AugerScatteringPoint.Z()*10000;
	}   
      N = Track.GetPrimaryIonizationV().size();
      for (int i=0; i<N; i++)
	{
	  TXYZ PrimaryIPoint = Track.GetPrimaryIonizationV()[i];
	  X[i] = PrimaryIPoint.X()*10000.0; 
	  Y[i] = PrimaryIPoint.Y()*10000.0;
	  Z[i] = PrimaryIPoint.Z()*10000.0;
	  nPa[i]=Track.GetnPairspath()[i];
	}
      Nd = Track.GetDiffElectronPosition().size();
      for (int i=0; i<Nd; i++)
	{
	  std::pair<double,double> position=Track.GetDiffElectronPosition()[i];
	  xe[i]=position.first*10000;
	  ye[i]=position.second*10000;
	}
            
      tree->Fill();  
    }
  
  TFile myFile(FileName,"RECREATE");
  tree->Write();
  myFile.Close();
  tree->Scan(); 
  
}

void ValidationStoppingPMix()
{
  const int N=60; 
  Double_t PhotonE[N];
  Double_t StoppingP[N];
  for (int i=0; i<60; i++)
    {
      PhotonE[i]=1.0+0.5*(i);
      StoppingP[i]=(Mixture->GetStoppingPower(PhotonE[i]));///Mixture->GetDensity();
      std::cout<<" energy = "<< PhotonE[i] <<" Stopping P =  "<< StoppingP[i]<<std::endl;
    }
  TGraph *GrafStoppingP = new TGraph(N,PhotonE,StoppingP);
  GrafStoppingP->GetXaxis()->SetTitle("Photon energy (keV)");
  GrafStoppingP->GetYaxis()->SetTitle("dE/dx (Kev/cm) ");//Kev cm^2/gr");
  GrafStoppingP->Draw("ap");
 
}

void TestReadout()
{ 
  TPhoton Photon(Source,ConversionPoint, rnd, 0);
  TTrack Track(&Photon, Mixture, rnd, 0);
  Track.PropagatePhotoelectron(); 
  Track.Drift();
  TReadout Gino(Dimension);
  
  TH2D *hxy = new TH2D("hxy","hxy",200,0,10000,200,0,10000);
  TH2D *hij = new TH2D("hij","h2",200,0,10000,200,0,10000);
  int N = Track.GetPrimaryIonizationV().size();
  double xe[800];
  double ye[800];
  std::pair<int,int> Pixel;
  std::pair<int,int> PixelCenterXY;//
   for (int i=0; i<N; i++)
     {   
      std::pair<double,double> position=Track.GetDiffElectronPosition()[i];
      xe[i]=position.first;//per far coincidere le origini dei sistemi di riferimento
      ye[i]=position.second;
      Pixel=Gino.PositiontoMatrix(xe[i],ye[i]);
      Gino.ChannelfromIJ(Pixel.first,Pixel.second); 
      int channeltest = Gino.GetChannel(); 
      std::pair<int,int> ijtest = Gino.IJfromChannel(channeltest);
      Gino.CenterHexagontoXY(Pixel);
      /* std::pair<int,int>*/ PixelCenterXY=Gino.GetCenterHexagontoXY();
      Gino.IJtoClusterChannel(Pixel.first,Pixel.second);
      std::pair<int,int> ClusterChannel =Gino.GetClusterChannelfromIJ();
      Gino.ClusterChanneltoIJ( ClusterChannel.first,ClusterChannel.second );
      std::pair<int,int> IJ =Gino.GetIJfromClusterChannel();
      hxy->Fill(xe[i]*10000+sqrt(3.)*80*161/4,ye[i]*10000.+80.*139./2.);
      hij->Fill(PixelCenterXY.first,PixelCenterXY.second);
      std::cout<<channeltest  <<" <-channel   "<< ijtest.first<< ijtest.second<<
	  "<-i, j " <<"Pixeli-> "<<Pixel.first<<" "<<Pixel.second<< "I="<<IJ.first<<" J="<<IJ.second<<std::endl;      
       }
  std::cout<< xe[5] <<"  "<< ye[5]<<" "<<std::endl;
  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  hij->Draw();
  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  hxy->Draw();
}    

void ROut()
{
  TPhoton Photon(Source,ConversionPoint, rnd, 0);
  TTrack Track(&Photon, Mixture, rnd, 0);
  Track.PropagatePhotoelectron(); 
  Track.Drift();
  TReadout Gino(Dimension);
  TGem Gem(rnd,Dimension);
  TH2D *hxy = new TH2D("hxy","hxy",200,0,10000,200,0,10000);
  TH2D *hij = new TH2D("hij","h2",200,0,10000,200,0,10000);
  double xe[800];
  double ye[800];
  std::pair<int,int> Pixel;
  std::pair<int,int> PixelCenterXY;
  std::vector<std::pair<double,double> > position=Track.GetDiffElectronPosition();
  for (int i=0; i<position.size();i++)
    {
      xe[i]=position[i].first;
      ye[i]=position[i].second;
      hxy->Fill(xe[i]*10000+sqrt(3.)*80*161/4,ye[i]*10000+80*139/2);
    }
 
  std::vector<std::pair<double,double> > SecE=Gem.DiffusionofSecondaryElectrons(position);//gemdiffusion
  std::cout<<"Secondarysize "<<SecE.size() <<std::endl;
  for (int i=0; i<SecE.size(); i++)
    {
      std::pair<double,double> PSecE=SecE[i];
      Pixel=Gino.PositiontoMatrix(PSecE.first,PSecE.second);
      Gino.CenterHexagontoXY(Pixel);
      PixelCenterXY=Gino.GetCenterHexagontoXY();
      hij->Fill(PixelCenterXY.first,PixelCenterXY.second);
      // std::cout<<"ss "<<PixelCenterXY.first<<" "<<PixelCenterXY.second<<std::endl;	
    }
   
  std::cout<< xe[5] <<"  "<< ye[5]<<" "<<std::endl;
  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  hij->Draw();
  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  hxy->Draw();
}     

void  Exa()
{
  ComputeConversionPoint();
  while (ConversionPoint.Z()<0.6)
    { 
      ComputeConversionPoint();
      }
  TPhoton Photon(Source,ConversionPoint, rnd, 0);
  TTrack Track(&Photon, Mixture, rnd, 0);
  Track.PropagatePhotoelectron();
  Track.Drift();
  TGem Gem(rnd,Dimension);
  std::vector<std::pair<double,double> >  DifEl=Track.GetDiffElectronPosition();
  std::vector<std::pair<double,double> >  SecEl=Gem.DiffusionofSecondaryElectrons(DifEl);
  TReadout Readout(Dimension);
  TDetector myDetector(Mixture,rnd,0,Dimension);
  std::vector<ADC> Digi=myDetector.mySampling(SecEl ,Readout);
  //std::cout<<Track.GetPhotoelectronTrueRange()<<std::endl;
  myDetector.myDrawSignal(Digi,Readout);
  
}

void Test(double a,double b)
{
  TReadout fufi(Dimension);
  std::pair<int,int> ij=  fufi.PositiontoMatrix(a,b);//(a+50./sqrt(3))/10000.,b/10000.);
  std::cout<<"ij di ab="<<ij.first<<"  "<<ij.second<<std::endl;
  fufi.CenterHexagontoXY(ij);
  std::pair<int,int> centerex=fufi.GetCenterHexagontoXY();
  fufi.ChannelfromIJ(ij.first,ij.second);
  int c = fufi.GetChannel();
  std::cout<<"center di ab "<<centerex.first<<"  "<<centerex.second<<std::endl;
  std::cout<<"channel"<< c <<std::endl;
  std::pair<int,int> ijb =   fufi.IJfromChannel(c);
  fufi.CenterHexagontoXY(ijb);
  std::pair<int,int> cHex=fufi.GetCenterHexagontoXY();
  std::cout<<"center di ab da channel "<<cHex.first<<"  "<<cHex.second<<std::endl;
}

void Center(int i, int j)
{
  TReadout fufi(Dimension);
  std::pair<int,int> ij=make_pair(i,j);
  fufi.CenterHexagontoXY(ij);
  std::pair<double,double> centerex=fufi.GetCenterHexagontoXY();
  std::cout<<"center "<<centerex.first<<"  "<<centerex.second<<std::endl;
}

void TestX()
{
  std::vector<std::pair<double,double> > Vect; 
  for (int i=0; i<1000; i++)
    { 
      double Xv=i/10000.+100./10000.;
      double Yv=i/10000.+100./10000.;
      std::pair<double,double> xx = make_pair(Xv,Yv);
      Vect.push_back(xx);
      }
  
  for (int i=0; i<1000; i++)
    { 
      double Xv=i/10000.+100./10000;
      double Yv=-(i/10000.+100./10000.)+1200./10000.;
      std::pair<double,double> xx = make_pair(Xv,Yv);
      Vect.push_back(xx);
    }
  TReadout Readout(Dimension);
  TGem Gem(rnd,Dimension);
  std::vector<std::pair<double,double> >  SecEl=Gem.DiffusionofSecondaryElectrons(Vect);
  TDetector myDetector(Mixture,rnd,0,Dimension);
  std::vector<ADC> Digi=myDetector.mySampling(SecEl,Readout);
  
  myDetector.myDrawSignal(Digi,Readout);
}

void TestChannel(int Channel)
{
   TReadout fufi(Dimension);
   std::pair<int,int> ij = fufi.IJfromChannel(Channel);
   std::cout<<"ij= "<<ij.first<<"  "<<ij.second<<std::endl;
   fufi.ChannelfromIJ(ij.first, ij.second);
   std::cout<<fufi.GetChannel()<<std::endl;
}

void TestCircle(double Z)
{
  std::vector<std::pair<double,double> > Vect; 
  double Sigma = Mixture->GetDiffusionSigma() * sqrt(Z)/10000;
  for (int i=0; i<1000; i++)
    { 
      double Xv=500./10000. * sin(i/(2.*kPI));
      double Yv=500./10000. * cos(i/(2.*kPI));
      Xv+=rnd->Gaus(0,Sigma);
      Yv+=rnd->Gaus(0,Sigma);
      std::pair<double,double> xx = make_pair(Xv,Yv);
      Vect.push_back(xx);
      }
  
  for (int i=0; i<1000; i++)
    { 
      double Xv=900./10000. * sin(i/(2.*kPI));
      double Yv=900./10000. * cos(i/(2.*kPI));
     
      
      Xv+=rnd->Gaus(0,Sigma);
      Yv+=rnd->Gaus(0,Sigma);
      std::pair<double,double> xx = make_pair(Xv,Yv);
      Vect.push_back(xx);
    }  
  TReadout Readout(Dimension);
  TGem Gem(rnd,Dimension);
  std::vector<std::pair<double,double> >  SecEl=Gem.DiffusionofSecondaryElectrons(Vect);
  TDetector myDetector(Mixture,rnd,0,Dimension);
  std::vector<ADC> Digi=myDetector.mySampling(SecEl,Readout);
  
  myDetector.myDrawSignal(Digi,Readout);
}

//Energia da sorgente con spettro a potenza inversa
double EnergyfromSpectrum(double SpectralIndex = 2.01)
{
  TF1 *Spectrum = new TF1("Spectrum","[0]*x^(-[1])",0.1,50);
  Spectrum->SetParameters(1, SpectralIndex);
  return Spectrum->GetRandom(2.,10.);
}

//Crea un tree di dati, nella stessa forma dei dati veri, contiene alcune informazioni montecarlo, gli angoli di emissione, l'energia del fotone e la miscela usata
void EventsTree( int Number=1000,double PEnergy = PhotonEnergy, TString File= "Eventstree",Double_t PolarDeg = PolarizationDegree)
{
  TTree *tree=new TTree("EventsTree","dat");
  //  std::cout<<"GEM PITCH = "<<Dimension->GetGem_Pitch()<<" PITCH = "<<Dimension->GetPitch()<<std::endl;
  int  Channel[1000];
  int  Charge[1000];
  int  DigiCharge[1000];
  double Phi;
  double Theta;
  int MID = MixID;
  TReadout Readout(Dimension);
  int Nevent;
  ADC Signal;
  TGem Gem(rnd,Dimension);
  std::cout<<" Gem Pitch =  "<<Gem.GetPitch()<<" ReadOut Pitch= " <<Readout.GetPitchR()<<" GEM Hole = "<<Gem.GetHole()<<std::endl;
  TDetector myDetector(Mixture,rnd,0,Dimension);
  int Clusterdim;
  int k;
  int ACheck;/////
  //Variabili IMAGING
  double XI=0;
  double YI=0;
  double Zconversion = 0.0;
  tree->Branch("Auger",&ACheck,"ACheck/I");///////
  tree->Branch("Nevent",&Nevent,"Nevent/I");
  tree->Branch("Clusterdim",&Clusterdim,"Clusterdim/I");
  tree->Branch("Channel",Channel,"Channel[Clusterdim]/I");
  tree->Branch("Charge",Charge,"Charge[Clusterdim]/I");
  tree->Branch("DigiCharge",DigiCharge,"DigiCharge[Clusterdim]/I");
  tree->Branch("InitialPhi",&Phi,"Phi/D");
  tree->Branch("InitialTheta",&Theta,"Theta/D");
  tree->Branch("MixID",&MID ,"MixID/I");
  tree->Branch("Energy",&PEnergy ,"PEnergy/D");
  //tree->Branch("Signal",&Signal.Channel,"Channel/I:Charge:DigiCharge");
  tree->Branch("Undetected",&k ,"k/I");
  //Branch per IMAGING
  tree->Branch("XInitial",&XI ,"XI/D");
  tree->Branch("YInitial",&YI ,"YI/D");
  tree->Branch("ZInitial",&Zconversion ,"Zconversion/D");
  //
 
  for (Int_t N=0; N<Number; N++ )
    { 
       if(N % 100==0)
	{ 
	  std::cout<<" ********** Nevent = "<<N<<std::endl;
	}
           
       k=0;

      ComputeConversionPoint(PEnergy);
      k = 1;
      while (ConversionPoint.Z()<0.6)
	{
	  ComputeConversionPoint(PEnergy);
	  k++;
	}

      TSource So(PEnergy, PolarizationAngle, PolarDeg, rnd, 0);
      So.SetPattern("UNIFORM");  //======>>>>  Uniforme(tra -0.007 e 0.007) vuol dire tra +-0.7mm 
      std::pair<double,double> XYConversion =  So.GetConversionXY();
      XI = XYConversion.first;
      YI = XYConversion.second;     
      Zconversion = ConversionPoint.Z();     
      //
      TPhoton Photon(&So,ConversionPoint, rnd, 0);
      TTrack Track(&Photon, Mixture, rnd, 0);
      Track.SetDimension (Dimension->GetGem_Radius()/2 ,
			  Dimension->GetGem_Radius()/2,
			  Dimension->GetZ_Drift(),
			  Dimension->GetZ_Gem() );
      Track.PropagatePhotoelectron();
      Phi = Track.GetPhotoelectronPhi();
      ACheck = Track.GetAugerCheck();////////////
      Theta = Track.GetPhotoelectronTheta();
      Track.Drift();
      std::vector<std::pair<double,double> >  DifEl=Track.GetDiffElectronPosition();  // X,Y elettrone della traccia dopo diffusione (calc in TTrack.Drift)
      std::vector<ADC> Digi= myDetector.mySampling(Gem.DiffusionofSecondaryElectrons(DifEl), Readout);
      if(Digi.size()>10)
	{
	  Clusterdim=Digi.size();
	  std::vector<ADC>::iterator pos;
	  int i = 0;
	  for (pos=Digi.begin(); pos!=Digi.end(); ++pos)
	    {
	      Signal=(*pos);
	      Channel[i]=Signal.Channel;
	      Charge[i]=Signal.Charge;
	      DigiCharge[i++]=Signal.DigiCharge;
	    }
	  tree->Fill();
	}
    }
 
  TString F = File;
  F += ".root";
  TFile FileT(F,"RECREATE");
  tree->Write();
  FileT.Close();
  TString FD = File;
  FD += ".txt";
  Dimension->WriteDimensionFile(FD);
}

void EventsTreeSp( int Number=1000,double PEnergy = PhotonEnergy, TString File= "Eventstree",Double_t PolarDeg = PolarizationDegree)
{
  TTree *tree=new TTree("EventsTree","dat");
  //  std::cout<<"GEM PITCH = "<<Dimension->GetGem_Pitch()<<" PITCH = "<<Dimension->GetPitch()<<std::endl;
  int  Channel[1000];
  int  Charge[1000];
  int  DigiCharge[1000];
  double Phi;
  double Theta;
  int MID = MixID;
  TReadout Readout(Dimension);
  int Nevent;
  ADC Signal;
  TGem Gem(rnd,Dimension);
  std::cout<<" Gem Pitch =  "<<Gem.GetPitch()<<" ReadOut Pitch= " <<Readout.GetPitchR()<<" Hole = "<<Gem.GetHole()<<std::endl;
  TDetector myDetector(Mixture,rnd,0,Dimension);
  int Clusterdim;
  int k;
  int ACheck;/////
  //Variabili IMAGING
  double XI=0;
  double YI=0;
  tree->Branch("Auger",&ACheck,"ACheck/I");///////
  tree->Branch("Nevent",&Nevent,"Nevent/I");
  tree->Branch("Clusterdim",&Clusterdim,"Clusterdim/I");
  tree->Branch("Channel",Channel,"Channel[Clusterdim]/I");
  tree->Branch("Charge",Charge,"Charge[Clusterdim]/I");
  tree->Branch("DigiCharge",DigiCharge,"DigiCharge[Clusterdim]/I");
  tree->Branch("InitialPhi",&Phi,"Phi/D");
  tree->Branch("InitialTheta",&Theta,"Theta/D");
  tree->Branch("MixID",&MID ,"MixID/I");
  tree->Branch("Energy",&PEnergy ,"PEnergy/D");
  //tree->Branch("Signal",&Signal.Channel,"Channel/I:Charge:DigiCharge");
  tree->Branch("Undetected",&k ,"k/I");
  //Branch per IMAGING
  tree->Branch("XInitial",&XI ,"XI/D");
  tree->Branch("YInitial",&YI ,"YI/D");
  //
  const int Ngraph =6;;
  double EnV[Ngraph] = {8.979,10,15,20,30,40};
  double CuAbsV[Ngraph] = {277,214,73.1,33.1,10.5,4.52};
  TGraph* CuAbsorption = new TGraph(Ngraph,EnV,CuAbsV);
  

  double Zdrift = Dimension->GetZ_Drift();
  for (Int_t N=0; N<Number; N++ )
    { 
      /////////////calcola energia
      // PEnergy = HerD->GetRandom(35,47);
      ////////////
      
      if(N % 100==0)
	{ 
	  std::cout<<" ********** Nevent = "<<N<<std::endl;
	}
      
      
      double Zconversion= 0.0;
      k=0;
      PEnergy=0;
      while (Zconversion<0.6 || Zconversion>0.6+Zdrift)
	{ 
	  
	  //PEnergy =  PhotonEnergy;
	  PEnergy =  Spectrum->GetRandom();
	  while (PEnergy<2) PEnergy=Spectrum->GetRandom();
	  
	  double Lambda=1./(Mixture->GetPhotoelectricCrossSection(PEnergy)*(Mixture->GetDensity()));
	  Zconversion = Dimension->GetZ_Drift() - rnd->Exp(Lambda);
	  //double GemEmis = rnd->Uniform(0,1);
	  
	  if (PEnergy>8.98 && Zconversion<0.6)
	    {
	      double lambdaCu = 1./(8.9*CuAbsorption->Eval(PEnergy));//in cm
	      if(rnd->Exp(lambdaCu) < 4e-4 && rnd->Uniform(0,1)<0.25)//Assorbimento+Fluorescenza(40%)*Prob.emissione nel gas(~50%)
		{
		  PEnergy=7.999;
		  Lambda=1./(Mixture->GetPhotoelectricCrossSection(PEnergy)*(Mixture->GetDensity()));
		  Zconversion = 0.6 + rnd->Exp(Lambda);//Dimension->GetZ_Drift() - rnd->Exp(Lambda);
		}
	      else Zconversion =0.0;
	    }
	  k++;
	}
      if (PEnergy==7.999) 
	{
	  std::cout<<" BINGO "<<std::endl;
	  PolarDeg=0;
	}
      else if ( PEnergy<5 || PEnergy>6.4 ) PolarDeg=0;
      else PolarDeg=0.99;
      PolarizationAngle = rnd->Gaus(90.,1);
      TSource So(PEnergy, PolarizationAngle, PolarDeg, rnd, 0);
      So.SetPattern("UNIFORM");
      std::pair<double,double> XYConversion =  So.GetConversionXY();
      XI = XYConversion.first;
      YI = XYConversion.second;
      TXYZ ConversionPoint(XI,YI,Zconversion);
	  
      TPhoton Photon(&So,ConversionPoint, rnd, 0);
      TTrack Track(&Photon, Mixture, rnd, 0);
      Track.SetDimension (Dimension->GetGem_Radius()/2 ,
			  Dimension->GetGem_Radius()/2,
			  Dimension->GetZ_Drift(),
			  Dimension->GetZ_Gem() );
      Track.PropagatePhotoelectron();
      Phi = Track.GetPhotoelectronPhi();
      ACheck = Track.GetAugerCheck();////////////
      Theta = Track. GetPhotoelectronTheta();
      Track.Drift();
      std::vector<std::pair<double,double> >  DifEl=Track.GetDiffElectronPosition();
      std::vector<ADC> Digi= myDetector.mySampling(Gem.DiffusionofSecondaryElectrons(DifEl), Readout);
      if(Digi.size()>10)
	{
	  Clusterdim=Digi.size();
	  std::vector<ADC>::iterator pos;
	  int i =0;
	  for (pos=Digi.begin(); pos!=Digi.end(); ++pos)
	    {
	      Signal=(*pos);//Digi[i];
	      Channel[i]=Signal.Channel;
	      Charge[i]=Signal.Charge;
	      DigiCharge[i++]=Signal.DigiCharge;
	    }
	  tree->Fill();
	}
    }
 
  TString F = File;
  F += ".root";
  TFile FileT(F,"RECREATE");
  tree->Write();
  FileT.Close();
  TString FD = File;
  FD += ".txt";
  Dimension->WriteDimensionFile(FD);
}

void TestCluster(int a,int b) 
{
  std::vector<ADC> Test;
  ADC uno;
  ADC due;
  uno.Charge=a;
  due.Charge=b;
  uno.Channel=10971;
  due.Channel=6001;
  Test.push_back(uno);
  Test.push_back(due);
  TReadout Readout(Dimension);
  TCluster Cluster;
  std::vector<Pixel> Testf = Cluster.ADCtoPixel(Test,Readout);
  Cluster.Analysis(Testf);
}
    
void ElasticPath()
{
  TH1D *Pat = new TH1D("Path","Path",1000,0,50);
  for (int i=0 ; i<200 ; i++)
    {
      double Pa = Mixture->GetElasticMeanFreePath(i*1./10.,"MOTT");
      Pa *=10000;
      Pat->Fill(Pa);
      std::cout<<Pa<<std::endl;
    }  
  Pat->Draw();
}

//Riapre il tree di GetTree e prosegue la simulazione fino alla lettura del Readout 
void EventsfromTreeTracks(TString FileName,TString FileInput ="treetracks"  )
{
  FileInput +=".root";
  TFile *fft= new TFile(FileInput);
  TTree *tft=(TTree*)fft->Get("Tree");
  /* tft->SetBranchStatus("*",0);
  tft->SetBranchStatus("xPrimaryElectronPosition",1);
  tft->SetBranchStatus("yPrimaryElectronPosition",1);
  tft->SetBranchStatus("Ndrift",1);
  tft->SetBranchStatus("InitialPhi",1);*/
  
  
  int nCluster=(int)tft->GetEntries();
  int  Channel[1000];
  int  Charge[1000];
  int  DigiCharge[1000];
  double Phi;
  double Theta;
  double XI=0;
  double YI=0;
  TReadout Readout(Dimension);
  Int_t Nevent;
  ADC Signal;
  TGem Gem(rnd,Dimension);
  TDetector myDetector(Mixture,rnd,0,Dimension);
  std::cout<<"GEM Pitch= "<<Dimension->GetGem_Pitch()<<" Pix Pitch = " <<Dimension->GetPitch()<<std::endl;
  int Clusterdim;
  Double_t xElDrifted[2000];
  Double_t yElDrifted[2000];
  int Sizediff;
  double InitialPhi;
  int MID = MixID;
  double PEnergy = PhotonEnergy;
 
  tft->SetBranchAddress("xPrimaryElectronPosition",xElDrifted);
  tft->SetBranchAddress("yPrimaryElectronPosition",yElDrifted);
  tft->SetBranchAddress("Ndrift",&Sizediff);
  tft->SetBranchAddress("InitialPhi",&InitialPhi);
  tft->SetBranchAddress("XInitial",&XI);
  tft->SetBranchAddress("YInitial",&YI);

  TTree *treeft = new TTree("EventsTree","dat");//"EventsTree","dat");
  treeft->Branch("Nevent",&Nevent,"Nevent/I");
  treeft->Branch("Clusterdim",&Clusterdim,"Clusterdim/I");
  treeft->Branch("Channel",Channel,"Channel[Clusterdim]/I");
  treeft->Branch("Charge",Charge,"Charge[Clusterdim]/I");
  treeft->Branch("DigiCharge",DigiCharge,"DigiCharge[Clusterdim]/I");
  treeft->Branch("InitialPhi",&Phi,"Phi/D");
  treeft->Branch("MixID",&MID ,"MixID/I");
  treeft->Branch("PhotonEnergy ",&PEnergy ,"PhotonE/D");
  treeft->Branch("XInitial",&XI ,"XI/D");
  treeft->Branch("YInitial",&YI ,"YI/D");
  for ( Nevent=0; Nevent<nCluster; Nevent++ )
   { 
     tft->GetEntry(Nevent);
     if (Nevent %100 == 0)std::cout<<"--------------"<<Nevent<<std::endl;
     Phi = InitialPhi;
     //XI = 0.0;
     //YI = 0.0;
     std::vector<std::pair<double,double> >  DifEl;
     for (int i = 0; i<Sizediff ; i++)
       {
	 //converte in cm le coordinate degli elettroni diffusi nel gas
	 std::pair<double,double> xyElDrifted=make_pair(xElDrifted[i]/10000,yElDrifted[i]/10000);
	 DifEl.push_back(xyElDrifted);
       }
     std::vector<ADC> Digi= myDetector.mySampling(Gem.DiffusionofSecondaryElectrons(DifEl), Readout);
     Clusterdim=Digi.size();
     std::vector<ADC>::iterator pos;
     int i =0;
     for (pos=Digi.begin(); pos!=Digi.end(); ++pos)
       {
	 Signal=(*pos);//Digi[i];
	 Channel[i]=Signal.Channel;
	 Charge[i]=Signal.Charge;
	 DigiCharge[i++]=Signal.DigiCharge;
       }
     treeft->Fill();
     
   }
  // treeft->Scan(); 
  TString FileNameD = FileName;
  FileName +=".root";
  TFile File(FileName,"RECREATE");//"Digital.root","RECREATE");
  treeft->Write();
  File.Close();
    FileNameD +=".txt";
  Dimension->WriteDimensionFile(FileNameD);
  delete tft;
  delete fft;
}

void MuvsEnergy(int Number = 1000)
{  
  double mu[13];
  double Energy[13];
  //TF1 *fufi = new TF1("fufi","[0]+[1]*(cos(x+[2]))^2",0,2*kPI);
  double Theta2;
  TReadout Readout(Dimension);
  ADC Signal;
  TGem Gem(rnd,Dimension);
  TDetector myDetector(Mixture,rnd,0,Dimension);
  TF1 *fufi = new TF1("fufi","[0]+[1]*(cos(x+[2]))^2",0,2*kPI);
  fufi->SetParLimits(2,0,kPI);
  //fufi->SetParameters(1.,1.,0.0);
  for (int PhotonE = 3 ; PhotonE <16 ; PhotonE++)
    {
      Source->SetEnergy(PhotonE);      
      //      TSource *Sour = new TSource(PhotonE, PolarizationAngle, PolarizationDegree, rnd, 0);
     
      TH1D *Theta = new TH1D("Theta","Theta",500,-1000,1000);
      for ( int Nevent=0; Nevent<Number; Nevent++ )
	{ 
	  while (ConversionPoint.Z()<0.6)
	    { 
	      ComputeConversionPoint();
	      //  k++;
	    }
	  TPhoton Photon(Source,ConversionPoint, rnd, 0);
	  TTrack Track(&Photon, Mixture, rnd, 0);
	  Track.PropagatePhotoelectron();
	  Track.Drift();
	  std::vector<std::pair<double,double> >  DifEl=Track.GetDiffElectronPosition();
	  std::vector<ADC> Digi= myDetector.mySampling(Gem.DiffusionofSecondaryElectrons(DifEl), Readout);
	  if(Digi.size()>10)
	    { 
	       TCluster Cluster(1.2 , 2.5 , 200);
	      std::vector<Pixel> D = Cluster.ADCtoPixel (Digi,Readout);
	      Cluster.Analysis(D);
	      Theta2 = Cluster.GetTheta2();
	      Theta->Fill(Theta2);
	    }
	 
	}
      Theta->Fit("fufi");
      double A    = fufi->GetParameter(0);
      double B    = fufi->GetParameter(1);
      double phi0 = fufi->GetParameter(2);
      mu[PhotonE-3]  = B / (2.*A+B);
      std::cout<< " mu = "<<mu[PhotonE-3]<<std::endl;
      Energy[PhotonE-3] = PhotonE*1.;
      std::cout<<" Energy =  "<<Energy[PhotonE-3]<<std::endl;
      delete Theta;
      
    }
   TGraph *Mu = new TGraph(13,Energy , mu);
   Mu->GetXaxis()->SetTitle("PhotonEnergy Kev");
   Mu->GetYaxis()->SetTitle("Mu");
   Mu->Draw("aL*");
}
	  
void Set(double a,double b)
{
  TGem Gem(rnd,Dimension);
  std::pair<double,double> aa = Gem.GemSampling(a,b);
  std::cout<<" X = "<<aa.first<<" Y = "<<aa.second<<std::endl;
}

void EnergyScan(int ClusterNumber = 1000,double EnergyMin=10,double EnergyMax=20)
{
  char Name[100];
  TString File ="KrDME6_4_drift1003cm3atm_Kev";
  double Energy = EnergyMin;
  while  (Energy <= EnergyMax)
    {
      std::cout<<"----------------Energy----------------"<<Energy<<std::endl;
      sprintf(Name,"DME_1cm0.5atm_Kev%.1f",Energy);
      EventsTree(ClusterNumber,Energy,Name);//EventsTree(ClusterNumber,Energy,FileName);
      std::cout<<Name<<""<<std::endl;
      Energy +=0.5;
    }
}

 void MixScan(int ClusterNumber ,int BEGIN,int END,double Energy)
{ 
  char Name[100];
  TString File ="KrDME6_4_drift1003cm3atm_Kev";
  MixID = BEGIN;
  while  (MixID <= END)
    {
      SetMixtureID(MixID);
     
      std::cout<<"----------------Energy----------------"<<Energy<<std::endl;
      sprintf(Name,"DME_1cm0.5atm_KevMIXID%.1i",MixID);
      EventsTree(ClusterNumber,Energy,Name);//EventsTree(ClusterNumber,Energy,FileName);
      std::cout<<Name<<""<<std::endl;
      MixID +=1;
    }
}     


/* sprintf(Name,"EventsTree_Kev.1f.root",Energy);*/
void  RuthScattering(TString ElementName ="C",double Energy = 1 )
{
  TCanvas *ca = new TCanvas ("ca","ca",1200,800);
  ca->Divide(2,2);
  TElement *Ele = new TElement(ElementName,rnd ,1);
  TH1D *h1=new TH1D("h1","distrAng",100,0.01,240);
  TH1D *h2=new TH1D("h2","distrAngR",100,0.01,240);
  const int N=10000; 
  Double_t Angle[N];
  Double_t AngleR[N];
  for (int i=0; i<10000; i++){
    Angle[i]= 360*(Ele->GetScatteringAngle(Energy,"MOTT"))/(kPI*2);
    AngleR[i]=360*(Ele->GetScatteringAngle(Energy,"RUTHERFORD"))/(kPI*2);
    h1->Fill(Angle[i]);
    h2->Fill(AngleR[i]);
  }
  h1->SetLineColor(2);
  ca->cd(1);
  h2->Draw();
  ca->cd(2);
  h1->Draw();
  ca->cd(3);
  h2->Draw();
  h1->Draw("same");
  // TF1 *FF = new TF1("FF","[0]*TMath::Gaus(x,0,[1])",0,200);
  //  h2->Fit("FF","E","",2,200);
}

void ParC(int Mixid,double Factor=45000)
{
  Mixture = new TGasMixture(MixID, rnd,Dimension,1);
  const int N = 200;
  double Pulse[N];
  double Parc[N];
  double Energy[N];
  double WI = Mixture->GetWIonization();
  for (int i=1; i<201 ; i++)
    {
      Pulse[i] = i*5100./201.;
      Energy[i] = Pulse[i] * (50./380.) * WI/1000.;
      Parc[i]=Factor * Mixture->GetElasticMeanFreePath(Energy[i], "MOTT");
    
    }
  TGraph *G = new TGraph(N,Pulse,Parc);
  G->Draw("ap");
}

void GetTreeBasic(double En = PhotonEnergy,int Nevents=10000, TString FileName = "Basictree.root")
{
  std::cout<<"Photon energy  = "<< En<<std::endl;
  TTree *tree=new TTree("Tree","datisim");
  Int_t nAuger; 
  Int_t NP;
  Int_t N;
  Int_t ntr;
  Int_t k;
  Int_t Nd;
  Double_t xe[2000];//dif
  Double_t ye[2000];//dif
  double Phi;
  double Theta;
  double ZC;
  double XI=0;
  double YI=0;
  int Ntr;
  double Energy;
  int Mixid;
  tree->Branch("NP",&NP,"NP/I");
  tree->Branch("InitialPhi",&Phi,"Phi/D");
  tree->Branch("Ntr",&Ntr,"Ntr/I");
  tree->Branch("InitialTheta",&Theta,"Theta/D");
  tree->Branch("N",&N,"N/I");
  tree->Branch("Ndrift",&Nd,"Nd/I");
  tree->Branch("xPrimaryElectronPosition",xe,"xe[Nd]/D");
  tree->Branch("yPrimaryElectronPosition",ye,"ye[Nd]/D");
  tree->Branch("ConversionZ",&ZC,"ZC/D");
  tree->Branch("XInitial",&XI ,"XI/D");
  tree->Branch("YInitial",&YI ,"YI/D");
  tree->Branch("Energy",&Energy,"Energy/D");
  tree->Branch("Mix",&Mixid,"Mixid/D");
  Mixid = MixID;
  for ( ntr=0; ntr<Nevents; ntr++ )
    { if (ntr %100 == 0)std::cout<<"--------------"<<ntr<<std::endl;
    
    ComputeConversionPoint(En);
    k=0;
    while (ConversionPoint.Z()<0.6)
      { 
	ComputeConversionPoint();
	k++;
      }
    ZC = ConversionPoint.Z();
    Energy=En;
    TSource So(Energy, PolarizationAngle, PolarizationDegree, rnd, 0);
    
    XI=ConversionPoint.X();
    YI=ConversionPoint.Y();
    /*So.SetPattern("UNIFORM");
    std::pair<double,double> XYConversion =  So.GetConversionXY();
    Ntr = ntr;
    XI = XYConversion.first;
    YI = XYConversion.second;
    TXYZ ConversionPointII(XI,YI,ZC);
    TPhoton Photon(&So,ConversionPointII, rnd, 0);
    */
    TPhoton Photon(&So,ConversionPoint, rnd, 0);
    TTrack Track(&Photon, Mixture, rnd, 0);
    Track.SetDimension (Dimension->GetGem_Radius()/2 ,
		Dimension->GetGem_Radius()/2,
		Dimension->GetZ_Drift(),
		Dimension->GetZ_Gem() );
      Track.PropagatePhotoelectron(); 
      Phi = Track.GetPhotoelectronPhi();
      Theta = Track.GetPhotoelectronTheta();
      Track.Drift();
      
      NP = Track.GetPhotoelectronScatteringV().size();
      N = Track.GetPrimaryIonizationV().size();
      Nd = Track.GetDiffElectronPosition().size();
      for (int i=0; i<Nd; i++)
	{
	  std::pair<double,double> position=Track.GetDiffElectronPosition()[i];
	  xe[i]=position.first*10000;
	  ye[i]=position.second*10000;
	}
      tree->Fill();  
    }
  TFile myFile(FileName,"RECREATE");
  tree->Write();
  myFile.Close();
}

void BasicScan()
{
  double Energy = 2;
  while (Energy<12)
    { 
      TString Name;
      char Dummy[100];
      sprintf(Dummy,"BasicTree7_%.1fkeV",Energy);
      Name =Dummy;
      Name+=".root";
      std::cout<<Name<<std::endl;
      GetTreeBasic(Energy,10000,Name);
      Energy+=2;
      }
}

void FromBasicScan()
{
  double Energy = 2;
  while (Energy<12)
    {
      SetDimensionMix(1,1 ,50,50);
      TString Name;
      char Dummy[100];
      sprintf(Dummy,"BasicTree2_%.1fkeV",Energy);
      Name =Dummy;
      //Name+=".root";
      TString FileName;
      sprintf(Dummy,"M2_5_%.1fkeV",Energy);
      FileName = Dummy;
      EventsfromTreeTracks(FileName,Name);
      Energy+=2;
    }
  
}

void RecFlux(TH1D *PulseH)
{
  const int N = PulseH->GetNbinsX();
  double Energy[N];
  double EffInv[N];
  double Eff[N];
  TH1D *SpREC = new TH1D("SpREC","SpREC",300,1.,25.);
  TH1D *EnH = new TH1D("EnH","EnH",300,1.,25.);
  
  double WIonization = Mixture->GetWIonization();
  double ElADC = Dimension->GetElectronsADC();
  double Avalance = Dimension->GetGain();
    for (int i=0; i<N;i++)
      {
	double Pulse = PulseH->GetBinCenter(i);
	Energy[i] = 0.4 + Pulse * (ElADC/Avalance) * WIonization/1000.;
	Eff[i]=Mixture->GetEfficiency(Energy[i]);
	EffInv[i] = 1./Eff[i];
	double W = EffInv[i] * PulseH->GetBinContent(i);
	SpREC->Fill(Energy[i],W);
	EnH->Fill(Energy[i]);
      }
    TGraph *gEff = new TGraph(N,Energy,Eff);
    TGraph *gEffInv = new TGraph(N,Energy,EffInv);
    TCanvas *c1 = new TCanvas("c1","c1",800,500);
    gEff->Draw("ap");
    TCanvas *c2= new TCanvas("c2","c2",800,500);
    gEffInv->Draw("ap");
    TCanvas *c3= new TCanvas("c3","c3",800,500);
    SpREC->Draw();
    TCanvas *c4= new TCanvas("c4","c4",800,500);
    SpREC->Draw();
}
