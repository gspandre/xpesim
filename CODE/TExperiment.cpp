#include "TExperiment.h"

TExperiment::TExperiment(TRandom3 *RND, TString suffix)
{
 rnd = RND;
 MaxGeneratedEvents=10000000;
 SetAstroSource=false;
 m_suffix=suffix;
}

void TExperiment::Delete()
{
  //  delete rnd;
  delete Dimension;
  delete Mixture;
  delete Source;
}
void  TExperiment::SetMixID(Int_t MIXID)
{
  
  MixID = MIXID;

}
void  TExperiment::SetSource(Double_t POLARIZATIONANGLE,Double_t POLARIZATIONDEGREE,Double_t  PHOTONENERGY)
{
  PolarizationAngle=POLARIZATIONANGLE;
  PolarizationDegree=POLARIZATIONDEGREE;
  PhotonEnergy = PHOTONENERGY;
}

void TExperiment::SetSource(Double_t POLARIZATIONANGLE,Double_t POLARIZATIONDEGREE,TH1D *obsSpectrum)
{
  PolarizationAngle=POLARIZATIONANGLE;
  PolarizationDegree=POLARIZATIONDEGREE;
  ObsSpectrum  = obsSpectrum;
  PhotonEnergy = ObsSpectrum->GetRandom();
  SetAstroSource=true;
}

void  TExperiment::SetThickness(double THICKNESS)
{
  Thickness=THICKNESS;
}

void  TExperiment::SetPressure(double PRESSURE)
{
  Pressure=PRESSURE;
}

void TExperiment::Generalsetup(double GP ,double PitchSet)
{
  //RandomGenerator = new TRandomGenerator();
  //rnd = RandomGenerator->GetRandomGenerator();
  Dimension = new TDimension();
  Dimension->SetZ_Drift(Thickness + Dimension->GetZ_Gem());
  Dimension->SetPressure(Pressure);
  //
  Pitch =PitchSet; 
  Dimension->SetGem_Pitch(GP);
  Dimension->SetPitch(Pitch);
  if (Pitch == 80)
    {
      Dimension->SetPixel(160,138);//PER GLORIA
    }
  else
    {
      Dimension->SetPixel(352,300);
    }
  Mixture = new TGasMixture(MixID, rnd,Dimension,1);
  Source = new TSource(PhotonEnergy, PolarizationAngle, PolarizationDegree, rnd, 0);
}



TString  TExperiment::GetName()
{
  TString Name= "Scan";
  char Dummy[100];
  std::vector<TCompound*> CompoundsVector  = Mixture->GetCompounds();
  std::vector<TCompound*>::iterator pos;
  for (pos = CompoundsVector.begin();pos!=CompoundsVector.end();++pos)
    {
      Name+=(*pos)->GetCompoundName();
      if(pos==(CompoundsVector.end()-1)) sprintf(Dummy,"_%.1f",(*pos)->GetFraction());
      else       sprintf(Dummy,"_%.1f-",(*pos)->GetFraction());
      Name+=Dummy;
    }
  sprintf(Dummy,"_%.1f",Thickness);Name+=Dummy;
  sprintf(Dummy,"cm_%.1fatm",Mixture->GetPressure());Name+=Dummy;
  return Name;
}
 
void TExperiment::EventsTree( int Number, TString File)
{
  TTree *tree=new TTree("EventsTree","dat");
  std::cout<<"GEM PITCH = "<<Dimension->GetGem_Pitch()<<" PITCH = "<<Dimension->GetPitch()<<" FirstDim =   "<<Dimension->GetPixel().first  <<std::endl;
  int  Channel[1500];
  int  Charge[1500];
  int  DigiCharge[1500];
  double Phi;
  double Theta;
  int MID = MixID;
  TReadout Readout(Dimension);
  int Nevent;
  ADC Signal;
  TGem Gem(rnd,Dimension);
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
  tree->Branch("Energy",&PhotonEnergy ,"PhotonEnergy/D");
  //tree->Branch("Signal",&Signal.Channel,"Channel/I:Charge:DigiCharge");
  tree->Branch("Undetected",&k ,"k/I");
  //Branch per IMAGING
  tree->Branch("XInitial",&XI ,"XI/D");
  tree->Branch("YInitial",&YI ,"YI/D");
  //
  int perc=0;  
  int perc_stop=0;
  int GeneratedEvents=0;
  if(SetAstroSource)
    perc_stop=TMath::Max(1,int(MaxGeneratedEvents/100));
  else
    perc_stop=TMath::Max(1,int(Number/100));
  int prc=0;
  for (Int_t N=0; N<Number; N++ )
    { 
      if(GeneratedEvents>=MaxGeneratedEvents)
	break;
      GeneratedEvents++;
      if(SetAstroSource)
	PhotonEnergy=ObsSpectrum->GetRandom();
      double PConv = GetEfficiencyMixture();
      double r1 = rnd->Uniform();
      //std::cout<<" energy: "<<PhotonEnergy<<" -> "<<PConv<<" rnd, "<<r1<<" ngen: "<<GeneratedEvents<<" proc "<<N<<std::endl;
      while (r1>=PConv) 
	{
	  if(SetAstroSource)
	    {
	      PhotonEnergy=ObsSpectrum->GetRandom();
	      PConv = GetEfficiencyMixture();
	    }
	  r1=rnd->Uniform();
	  GeneratedEvents++;
	  //	  std::cout<<" energy: "<<PhotonEnergy<<" -> "<<PConv<<" rnd, "<<r1<<" ngen: "<<GeneratedEvents<<" proc "<<N<<std::endl;
	}
      /////////////calcola energia
      // PhotonEnergy = HerD->GetRandom(35,47);
      ////////////
      
      if(SetAstroSource)
	{
	  if(GeneratedEvents * 100/MaxGeneratedEvents > prc)
	    {
	      prc++;
	      perc = GeneratedEvents/MaxGeneratedEvents*100;
	      std::cout<<" *** "<<prc<<" % done, ** N event conv. = "<<N<<", tot. on the mirror area: "<<GeneratedEvents<<std::endl;
	    }
	}
      else
	{
	  if(N * 100 % Number == 0)
	    {
	      perc = N/Number*100;
	      std::cout<<" *** "<<perc<<" % done, ** N event conv. = "<<N<<", tot. on the mirror area: "<<GeneratedEvents<<std::endl;
	    }
	}
      
      Source->SetPattern("UNIFORM");
      std::pair<double,double> XYConversion =  Source->GetConversionXY();
      XI = XYConversion.first;
      YI = XYConversion.second;
      TXYZ ConversionPoint(XI,YI,0.0);
      double Lambda=1./(Mixture->GetPhotoelectricCrossSection(PhotonEnergy)*(Mixture->GetDensity()));
      double Zconversion= 0.0;
      k=0;
      while (ConversionPoint.Z()<0.6)
	{ 
	  Zconversion = Dimension->GetZ_Drift() - rnd->Exp(Lambda);
	  ConversionPoint.SetXYZ(XI,YI,Zconversion);
	  k++;
	}
      Nevent = N;
      
      
      TPhoton Photon(Source,ConversionPoint, rnd, 0);
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
      
      // commento start here:
      //std::vector<TXYZ> Digi0=Track.GetPrimaryIonizationV();
      //	      std::vector<std::pair<double,double> > pippo;
      //      std::vector<TXYZ>::iterator pos0;
      //     for (pos0=Digi0.begin(); pos0!=Digi0.end(); ++pos0)
      //    {
      //    pippo.push_back(std::make_pair((*pos0).X(),(*pos0).Y()));
      //    }
      // commento end here:
      
      std::vector<ADC> Digi= myDetector.mySampling(Gem.DiffusionofSecondaryElectrons(DifEl), Readout);
      if(Digi.size()>10)//DA 10 a 1
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
  //  tree->Scan(); 
  char Dummy[100];
  TString DirName = "../Trees/";
  TString name = GetName();
  DirName += name;
  TString Command = "mkdir ";
  Command+=DirName;
  gSystem->Exec(Command);
  DirName += "/";
  DirName += File;
  TString F = DirName;
  if(!SetAstroSource)
    {
      sprintf(Dummy,"_%.1fkeV",PhotonEnergy); 
      F+=Dummy;
      F+=m_suffix;
      F+=".root";
    }
  else
    {
      F+="AstroScan";
      F+=m_suffix;
      F+=".root";
    }
  
  TFile FileT(F,"RECREATE");
  tree->Write();
  FileT.Close();
  TString FD = DirName;
  if(!SetAstroSource)
    {
      FD+=Dummy;
      FD += ".txt";
    }
  else
    {
      FD += "AstroScan";
      FD+=m_suffix;
      FD+=".txt";

    }
  Dimension->WriteDimensionFile(FD);
  delete tree;
}

void     TExperiment::AnalyzeTree()
{
  TString Name = "../Trees/";
  Name+= GetName();
  
  if(!SetAstroSource)
    {
      char Dummy[100];
      if (Pitch==50)
	sprintf(Dummy,"/Eventstree_%.1fkeV",PhotonEnergy); 
      else if (Pitch==80)
	sprintf(Dummy,"/Pixmap8_%.1fkeV",PhotonEnergy); //modifica per Pixmap 80
      else std::cout<<"!!Pitch unknow!!  "<<std::endl;
      Name+=Dummy;
    }
  else
    {
      Name+="/EventstreeAstroScan";
      Name+=m_suffix;
    }
  TTreeAnalysis A(Name,"EventsTree");
  A.TreeAnalysis();
  A.FillH();
  A.DeleteAll();
}

void TExperiment::ReadResult(double &mucut,double &Eff,double &Effcut, double &mu, double &FitProbability, double  &FitProbabilityCut, double &mufirst)
{
  TString Name = "../Trees/";
  Name+= GetName();
  char Dummy[100];
  if (Pitch==50)
    sprintf(Dummy,"/Eventstree_%.1fkeVAnalysis.root",PhotonEnergy); 
  else if (Pitch==80)
    sprintf(Dummy,"/Pixmap8_%.1fkeVAnalysis.root",PhotonEnergy);  //modifica per Pixmap 80
  else std::cout<<"!!Pitch unknow!!  "<<std::endl;
  Name+=Dummy;
  TFile file(Name);
  
  TH1D *ModulationFirstStepH = (TH1D*) file.Get("ModulationFirstStep");  
  TH1D *ModulationCutH = (TH1D*) file.Get("ModulationCut");
  TH1D *ModulationSimpleH = (TH1D*)file.Get("ModulationSimple");

  TF1 *fuficut = (TF1*) ModulationCutH->GetListOfFunctions()->FindObject("fufi");
  double A2    = fuficut->GetParameter(0);
  double B2    = fuficut->GetParameter(1);
  double phi02 = fuficut->GetParameter(2);
  double TotalEvents = ModulationSimpleH->GetEntries();
  FitProbabilityCut = fuficut->GetProb();
  mucut   = 100.0 * B2 / (2.*A2+B2);  
  std::cout<<"mu2 Second= "<<mucut<<std::endl;
  //delete ModulationCutH;
  // delete ModulationSimpleH;
  // delete fufi;

  TF1 *fufi = (TF1*) ModulationSimpleH->GetListOfFunctions()->FindObject("fufi");
  double A    = fufi->GetParameter(0);
  double B    = fufi->GetParameter(1);
  double phi  = fufi->GetParameter(2);
  Eff = GetEfficiencyMixture()*100.;
  Effcut = Eff* (ModulationCutH->GetEntries())/TotalEvents;
  FitProbability = fufi->GetProb();
  mu   = 100.0 * B / (2.*A+B);  
  //  delete fufi;

  TF1 *fufif = (TF1*) ModulationFirstStepH->GetListOfFunctions()->FindObject("fufi");
  double Af    = fufif->GetParameter(0);
  double Bf    = fufif->GetParameter(1);
  double phif  = fufif->GetParameter(2);
  mufirst   = 100.0 * Bf / (2.*Af+Bf);
  file.Close();

}

double TExperiment::GetEfficiencyMixture(double myene)
{
  if(myene==0.0)
    return Mixture->GetEfficiency(PhotonEnergy);
  return Mixture->GetEfficiency(myene);
  
}

TGraph  *TExperiment::GetEfficiencyGraph()
{
  double emin    = 1.0;
  double emax    = 100.0;
  int N=1000;
  double *e = new double[N];
  double *eFF = new double[N];
  double de = (emax-emin)/(1.0*N);
  for(int i=0; i<N; i++)
    {
      e[i]=emin+de*i;
      eFF[i]=GetEfficiencyMixture(e[i]);
      std::cout<<e[i]<<" "<<eFF[i]<<std::endl;
    }
  return new TGraph(N,e,eFF);
}



TH1D *TExperiment::GetClusterH(double PEnergy)
{
  TString Name = "../Trees/";
  Name+= GetName();
  char Dummy[100];
  sprintf(Dummy,"/Eventstree_%.1fkeVAnalysis.root",PEnergy); 
  Name+=Dummy;
  TFile file(Name);
  TH1D *ClusterH =(TH1D*)file.Get("ClusterDim");
  ClusterH->SetDirectory(0);
  file.Close();
  return ClusterH;
}

TH1D *TExperiment::GetRecordDimH(double PEnergy)
{
  TString Name = "../Trees/";
  Name+= GetName();
  char Dummy[100];
  sprintf(Dummy,"/Eventstree_%.1fkeVAnalysis.root",PEnergy); 
  Name+=Dummy;
  TFile file(Name);
  TH1D *RecordDimH =(TH1D*)file.Get("RecordDimH");
  RecordDimH->SetDirectory(0);
  file.Close();
  return RecordDimH;
}
