#include "TExperiment.h"

TExperiment::TExperiment(TRandom *RND, TString suffix)
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
  Source = new TSource(PhotonEnergy, PolarizationAngle, PolarizationDegree, rnd, 0);
}

void TExperiment::SetSource(Double_t POLARIZATIONANGLE,Double_t POLARIZATIONDEGREE,TH1D *obsSpectrum)
{
  PolarizationAngle=POLARIZATIONANGLE;
  PolarizationDegree=POLARIZATIONDEGREE;
  ObsSpectrum  = obsSpectrum;
  PhotonEnergy = ObsSpectrum->GetRandom();
  //cout << "SetSource Spect: energy 0 " << PhotonEnergy << endl;///<<++++++++
  Source = new TSource(PhotonEnergy, PolarizationAngle, PolarizationDegree, rnd, 0);
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
  Dimension = new TDimension();
  Dimension->SetZ_Drift(Thickness + Dimension->GetZ_Gem());
  Dimension->SetPressure(Pressure);
  Pitch = PitchSet; 
  Dimension->SetGem_Pitch(GP);
  Dimension->SetPitch(Pitch); 
  if (Pitch == 80)
    {
      Dimension->SetPixel(160,138);
    }
  else
    {
      Dimension->SetPixel(352,300);
    }
  std::cout << "General SETUP: " << std::endl;
  std::cout << "GEM Diameter  (cm) = " << Dimension->GetGem_Radius() << std::endl;
  std::cout << "GEM Pitch  (micron) = " << Dimension->GetGem_Pitch() << std::endl;
  std::cout << "ASIC Pitch (micron) = " << Dimension->GetPitch() << std::endl;
  std::cout << "Pixels Matrix = " << Dimension->GetPixel().first << " x " << Dimension->GetPixel().second << std::endl;
  std::cout << "Absorption Gap (cm) = " << Dimension->GetZ_Drift() - Dimension->GetZ_Gem() << std::endl;
  std::cout << "Transfer Gap   (cm) = " << Dimension->GetZ_Gem() << std::endl;

  Mixture = new TGasMixture(MixID,rnd,Dimension,0);
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

TString  TExperiment::GetNameforFile()
{
  TString Name= "Scan";
  char Dummy[100];
  std::vector<TCompound*> CompoundsVector  = Mixture->GetCompounds();
  std::vector<TCompound*>::iterator pos;
  for (pos = CompoundsVector.begin();pos!=CompoundsVector.end();++pos)
    {
      Name+=(*pos)->GetCompoundName();
      if(pos==(CompoundsVector.end()-1)) sprintf(Dummy,"(%d)",(int)((*pos)->GetFraction()*100));
      else       sprintf(Dummy,"(%d)-",(int)((*pos)->GetFraction()*100));
      Name+=Dummy;
    }  
  sprintf(Dummy,"_%.1fatm",Mixture->GetPressure());Name+=Dummy;
  sprintf(Dummy,"_%.1fcm",Thickness);Name+=Dummy;
  return Name;
}
 
TString  TExperiment::GetMixType()
{
  TString Name;
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
  sprintf(Dummy," %.1f",Thickness);Name+=Dummy;
  sprintf(Dummy,"cm %.1fatm",Mixture->GetPressure());Name+=Dummy;
  return Name;
}

//*******************************************************************************
/// \brief Loop over G4 generated hit and extract readout hit                   *
//*******************************************************************************
void TExperiment::G4MCEventsLoop(int Number, TString File)
{/* This fuction emulates the EventsTree function 
    except for the photoelectron propagation.
    A lot of stuff is duplicated and hard coded, need to cleanup!!!
 */

  
  // Input Geant4 root file and get the list of hit and energy
  TFile *G4File = new TFile(File);
  TTree *G4Tree = (TTree*) G4File->Get("Events");
  int G4Entries = G4Tree->GetEntries() ;
  cout << "G4Entries " << G4Entries  << endl;
  Int_t G4EventID, G4DU;
  Double_t G4X_Hit, G4Y_Hit, G4Z_Hit, G4Energy;
  G4Tree->SetBranchAddress("EventID", &G4EventID);
  G4Tree->SetBranchAddress("DetUnitID", &G4DU);
  G4Tree->SetBranchAddress("X_Hit", &G4X_Hit);
  G4Tree->SetBranchAddress("Y_Hit", &G4Y_Hit);
  G4Tree->SetBranchAddress("Z_Hit", &G4Z_Hit);
  G4Tree->SetBranchAddress("EnergyDeposit", &G4Energy);
  
  // Output root file setup
  int MID = MixID;

  TString F = File;
  F.ReplaceAll(".root","_GPDMC.root");
  TFile FileT(F,"RECREATE");
  
  TTree *tree=new TTree("EventsTree","dat");
  int  Channel[1500];
  int  Charge[1500];
  int  DigiCharge[1500];
  double Phi;
  double Theta;
  double PhotoelectronPracticalRange, PhotoelectronEnergy, AugerPracticalRange, AugerEnergy;
  double PhotoelectronTrueRange, AugerTrueRange;
  double XI, YI, ZI;
  int Nevent, Clusterdim, ACheck, k;
  Int_t NPrimaryElectrons;
  Int_t DU, EvtId;
  tree->Branch("Auger",&ACheck,"ACheck/I");///////
  tree->Branch("Nevent",&Nevent,"Nevent/I");
  tree->Branch("Clusterdim",&Clusterdim,"Clusterdim/I");
  tree->Branch("Channel",Channel,"Channel[Clusterdim]/I");
  tree->Branch("Charge",Charge,"Charge[Clusterdim]/I");
  tree->Branch("DigiCharge",DigiCharge,"DigiCharge[Clusterdim]/I");
  tree->Branch("InitialPhi",&Phi,"Phi/D");
  tree->Branch("InitialTheta",&Theta,"Theta/D");
  tree->Branch("PhotoelPracticalRange",&PhotoelectronPracticalRange,"PhotoelectronPracticalRange/D");
  tree->Branch("PhotoelTrueRange",&PhotoelectronTrueRange,"PhotoelectronTrueRange/D");
  tree->Branch("PhotoelectronEnergy",&PhotoelectronEnergy,"PhotoelectronEnergy/D");
  tree->Branch("AugerPracticalRange",&AugerPracticalRange,"AugerPracticalRange/D");
  tree->Branch("AugerTrueRange",&AugerTrueRange,"AugerTrueRange/D");
  tree->Branch("AugerEnergy",&AugerEnergy,"AugerEnergy/D");
  tree->Branch("MixID",&MID ,"MixID/I");
  tree->Branch("Energy",&PhotonEnergy ,"PhotonEnergy/D");
  tree->Branch("Undetected",&k ,"k/I");
  tree->Branch("XInitial",&XI ,"XI/D");
  tree->Branch("YInitial",&YI ,"YI/D");
  tree->Branch("ZInitial",&ZI ,"ZI/D");
  tree->Branch("NPrimaryElectrons",&NPrimaryElectrons ,"NPrimaryElectrons/I");
  tree->Branch("DU",&DU ,"DU/I");
  tree->Branch("EvtId",&EvtId ,"EvtId/I");
  
  // Need a fake source and photon
  Source = new TSource(1, 0, 0, rnd, 0);
  TXYZ ConversionPoint(0.,0.,0.);
  TPhoton Photon(Source, ConversionPoint, rnd, 0);
  
  // Need detector
  TGem Gem(rnd,Dimension);
  TDetector myDetector(Mixture,rnd,0,Dimension);
  TReadout Readout(Dimension);
  ADC Signal;
  
  // Init vars and start loop
  std::vector<TXYZ> G4PrimaryIonizationV;
  Int_t G4CurrentEvtId = -1;
  Int_t G4NumProcessedEvt = 0;
  int small_size = 0;
  for (int i=0; i<G4Entries; i++)
    {
      G4Tree->GetEntry(i);
      
      if (G4EventID != G4CurrentEvtId) // new event: process previous one and reset variables
	{
	  // if not the first ID, process the event
	  if (G4CurrentEvtId !=-1)
	    {
	      if (G4NumProcessedEvt%10 ==0) {
		cout << "G4MC:: Processing " << G4NumProcessedEvt << " evts. Last Id: " << G4CurrentEvtId << endl;}
	      //cout << "G4MC:: Processing event id " << G4CurrentEvtId << " (" << G4NumProcessedEvt << " evt)" << endl;
	      G4NumProcessedEvt ++;
	      // DO THE PROCESSING HERE!!!!!
	      DU = G4DU;
	      TTrack Track(&Photon, Mixture, rnd, 0);
	      Track.SetDimension (Dimension->GetGem_Radius()/2 ,
				  Dimension->GetGem_Radius()/2,
				  Dimension->GetZ_Drift(),
				  Dimension->GetZ_Gem() );
	      Track.SetPrimaryIonizationV(G4PrimaryIonizationV); //
	      Track.Drift();
	      std::vector<std::pair<double,double> >  DifEl=Track.GetDiffElectronPosition();      
	      std::vector<ADC> Digi= myDetector.mySampling(Gem.DiffusionofSecondaryElectrons(DifEl), Readout);
	      //cout << "G4MC:: CS DIGI SIZE: " << Digi.size() << " DU " << DU <<endl;
	      
	      if(Digi.size()>0)
		{
		  EvtId = G4CurrentEvtId;
		  Clusterdim=Digi.size(); 
		  std::vector<ADC>::iterator pos;
		  int i =0;
		  for (pos=Digi.begin(); pos!=Digi.end(); ++pos)
		    {
		      Signal=(*pos);//Digi[i];
		      if (Signal.Channel>0) // valid channel
			{
			  Channel[i]=Signal.Channel;
			  //cout << "CS Signal.Channel: " << Signal.Channel << endl;
			  Charge[i]=Signal.Charge;
			  DigiCharge[i++]=Signal.DigiCharge;
			}
		    }
		  tree->Fill();
		}
	      else small_size++;
	      
	    }

	  // check if reached the maximum number of events:
	  //cout << " processed evt " << G4NumProcessedEvt << " evt" << endl;
	  if (G4NumProcessedEvt==Number) {break;}

	  // move to next one...
	  //cout << "G4MC:: New Event found, id= " << G4EventID << endl;
	  G4CurrentEvtId = G4EventID;
	  
	  // reset ionization vector:
	  G4PrimaryIonizationV.clear();

	  
	}
      else // fill event hit list
	{
	  //cout << "Id " <<i <<  " G4EventID " << G4EventID << " G4 Hit (X,Y,Z,E) ";
	  //cout << G4X_Hit <<" " << G4Y_Hit <<" " << G4Z_Hit << " " << G4Energy << endl;

	  //
	  // conversion from G4 geom to MC
	  //
	  
	  /* //First version:
	  double Z_MC = G4Z_Hit+174.534+0.5+0.06; // center in Gas cell @ -174.534, add half thickness, add transfer gap; all in cm
	  double X_MC = -99.;
	  double Y_MC = -99.;
	  if (G4DU == 0) //, for DU 0:
	    {
	      X_MC = G4X_Hit; // this is the same, in cm
	      Y_MC = G4Y_Hit-28.; // center in Gas cell @ +28; all in cm
	    }
	  if (G4DU == 1)//, for DU 1: rotating -120 deg puts du1 in the position of du 0
	    {
	      X_MC = (-0.5*G4X_Hit + (sqrt(3)/2)*G4Y_Hit); 
	      Y_MC = ((-sqrt(3)/2)*G4X_Hit -0.5*G4Y_Hit)-28.; 
	    }
	  if (G4DU == 2)//, for DU 1: rotating 120 deg puts du2 in the position of du 0
	    {
	      X_MC = (-0.5*G4X_Hit + (-sqrt(3)/2)*G4Y_Hit); 
	      Y_MC = ((sqrt(3)/2)*G4X_Hit -0.5*G4Y_Hit)-28.; 
	    }
	  */

	  // Second version, coordinate il local (Gas Cell) coordinate
	  double Z_MC = G4Z_Hit+0.5+0.06; // add half thickness, add transfer gap; all in cm
	  double X_MC = G4X_Hit; // already centered in 0
	  double Y_MC = G4Y_Hit;
	  
	  double E_MC = G4Energy; // assume Energy is in keV (as needed in the MC)
	  // Get number of electron/ion pair in this step: pure poisson fluctuations
	  int nAvePairs = 1000*E_MC/Mixture->GetWIonization();
	  int nPairs    = rnd->Poisson(nAvePairs);
	  // fill the ionization vector
	  for (Int_t i=0; i<nPairs; i++){
	    if (std::fabs(X_MC)<0.8 && std::fabs(Y_MC)<0.8) // include only hit close to the active region (which is 15x15 mm)
	      {
		TXYZ IonizationPoint = TXYZ(X_MC, Y_MC, Z_MC);
		G4PrimaryIonizationV.push_back(IonizationPoint);
	      }
	  }
	  
	  //cout << " IN MC GEOMETRY: (X,Y,Z,E) ";
	  //cout << X_MC <<" " << Y_MC <<" " << Z_MC << " " << E_MC << " NPairs " << nPairs
	  //     << " V.size " << G4PrimaryIonizationV.size() << endl; 
	}
    } // end of evt loop
  cout << "G4MC::  Discarded events with 0 pixels: " << small_size << endl;

  
  tree->Write();
  FileT.Close();
  //delete tree;

}

void TExperiment::EventsTree(int Number, TString File)
{
  bool MixEffON = false;
  TTree *tree=new TTree("EventsTree","dat");
  int  Channel[1500];
  int  Charge[1500];
  int  DigiCharge[1500];
  double Phi;
  double Theta;
  double PhotoelectronPracticalRange, PhotoelectronEnergy, AugerPracticalRange, AugerEnergy;
  double PhotoelectronTrueRange, AugerTrueRange;
  int MID = MixID;
  TReadout Readout(Dimension);
  int Nevent;
  ADC Signal;
  TGem Gem(rnd,Dimension);
  TDetector myDetector(Mixture,rnd,0,Dimension);
  int Clusterdim;
  int k;
  int ACheck;
  Int_t NPrimaryElectrons;
  //Variabili IMAGING
  double XI=0;
  double YI=0;
  double ZI=0;
  tree->Branch("Auger",&ACheck,"ACheck/I");///////
  tree->Branch("Nevent",&Nevent,"Nevent/I");
  tree->Branch("Clusterdim",&Clusterdim,"Clusterdim/I");
  tree->Branch("Channel",Channel,"Channel[Clusterdim]/I");
  tree->Branch("Charge",Charge,"Charge[Clusterdim]/I");
  tree->Branch("DigiCharge",DigiCharge,"DigiCharge[Clusterdim]/I");
  tree->Branch("InitialPhi",&Phi,"Phi/D");
  tree->Branch("InitialTheta",&Theta,"Theta/D");
  tree->Branch("PhotoelPracticalRange",&PhotoelectronPracticalRange,"PhotoelectronPracticalRange/D");
  tree->Branch("PhotoelTrueRange",&PhotoelectronTrueRange,"PhotoelectronTrueRange/D");
  tree->Branch("PhotoelectronEnergy",&PhotoelectronEnergy,"PhotoelectronEnergy/D");
  tree->Branch("AugerPracticalRange",&AugerPracticalRange,"AugerPracticalRange/D");
  tree->Branch("AugerTrueRange",&AugerTrueRange,"AugerTrueRange/D");
  tree->Branch("AugerEnergy",&AugerEnergy,"AugerEnergy/D");
  tree->Branch("MixID",&MID ,"MixID/I");
  tree->Branch("Energy",&PhotonEnergy ,"PhotonEnergy/D");
  tree->Branch("Undetected",&k ,"k/I");
  tree->Branch("XInitial",&XI ,"XI/D");
  tree->Branch("YInitial",&YI ,"YI/D");
  tree->Branch("ZInitial",&ZI ,"ZI/D");
  tree->Branch("NPrimaryElectrons",&NPrimaryElectrons ,"NPrimaryElectrons/I");
  //
  int nbin = 1000;
  double xmin = 0;
  double xmax = 6000;
  int nbint = 100;
  double xmint = 0;
  double xmaxt = 200;
  pheRange = new TH1D("PhotoelPracticalRange","PhotoelPracticalRange",nbin,xmin,xmax);
  pheRange->GetXaxis()->SetTitle("Practical range (micron)");
  augRange = new TH1D("AugerPracticalRange","AugerPracticalRange",nbint,xmint,xmaxt);
  augRange->GetXaxis()->SetTitle("Practical range (micron)");
  pheTrueRange = new TH1D("PhotoelTrueRange","Photoel. total track lenght",nbin,xmin,xmax);
  pheTrueRange->GetXaxis()->SetTitle("True range (micron)");
  augTrueRange = new TH1D("AugerTrueRange","Auger total track lenght",nbint,xmint,xmaxt);
  augTrueRange->GetXaxis()->SetTitle("True range (micron)");
  
  //
  int perc=0;  
  int GeneratedEvents=0;
  int prc=0;
  int N;
  int small_size = 0;
  for (N=0; N<Number; N++ )
    { 
      gSystem->ProcessEvents();
      if(GeneratedEvents>=MaxGeneratedEvents)
      	break;
      GeneratedEvents++;
      if(SetAstroSource) PhotonEnergy=ObsSpectrum->GetRandom();

      if(MixEffON){
      // commentare per avere stessi risultati di SimPolScan
      double PConv = GetEfficiencyMixture(PhotonEnergy);
      double r1 = rnd->Uniform(); 
      //cout << "EventsTree: Energy 0" << PhotonEnergy << endl;///<<++++++++      
      while (r1>=PConv) 
	{
	  if(SetAstroSource)
	    {
	      PhotonEnergy=ObsSpectrum->GetRandom();
	      PConv = GetEfficiencyMixture(PhotonEnergy);
	    }
	  r1=rnd->Uniform();        // commentare per avere stessi risultati di SimPolScan
	  GeneratedEvents++;
	}
       //cout << "EventsTree: Energy 1 " << PhotonEnergy << endl;///<<++++++++         
      }
      else GeneratedEvents++;

      if(SetAstroSource)
	{
	  if(N%1000 == 0)
	    {
	      perc = (int)((float)GeneratedEvents/MaxGeneratedEvents*100);
	      std::cout<<" *** "<< prc <<" % done, ** N event conv. = "<< 
		N <<", tot. on the mirror area: "<< GeneratedEvents <<std::endl;
	    }
	}
      else {

	if(N%1000 == 0)
	  {
	    perc = (int)((float)N/Number*100);	
	    std::cout<<" *** "<<perc<<" % done, ** N event conv. = "<<N<<", tot. on the mirror area: "<<GeneratedEvents<<std::endl;
	  }
      }
      Source->SetPattern("UNIFORM"); // 2 options: UNIFORM, GAUS
      std::pair<double,double> XYConversion =  Source->GetConversionXY();
      XI = XYConversion.first;
      YI = XYConversion.second;
      TXYZ ConversionPoint(XI,YI,0.0);
      
      double Lambda=1./(Mixture->GetPhotoelectricCrossSection(PhotonEnergy)*(Mixture->GetDensity()));
      double Zconversion= 0.0;
      k=0;
      
      while (ConversionPoint.Z()<0.06)  // 0.6 is the transfer gap (between GEM and ASIC)
	{ 
	  Zconversion = Dimension->GetZ_Drift() - rnd->Exp(Lambda);
	  ConversionPoint.SetXYZ(XI,YI,Zconversion);
	  k++;
	}
      Nevent = N;
      ZI = Zconversion;
      TPhoton Photon(Source,ConversionPoint, rnd, 0); // check which ph energy here!!!

      //cout << "EventsTree: Energy 2 " << PhotonEnergy << endl;///<<++++++++

      TTrack Track(&Photon, Mixture, rnd, 0);
      Track.SetDimension (Dimension->GetGem_Radius()/2 ,
			  Dimension->GetGem_Radius()/2,
			  Dimension->GetZ_Drift(),
			  Dimension->GetZ_Gem() );
      Track.PropagatePhotoelectron();
      Phi = Track.GetPhotoelectronPhi();
      ACheck = Track.GetAugerCheck();
      Theta = Track.GetPhotoelectronTheta();
      PhotoelectronPracticalRange = Track.GetPhotoelectronPraticalRange();
      PhotoelectronEnergy = Track.GetPhotoelEnergy();
      AugerPracticalRange = Track.GetAugerPraticalRange();
      AugerEnergy = Track.GetAugerEnergy();
      PhotoelectronTrueRange = Track.GetPhotoelectronTrueRange();
      AugerTrueRange = Track.GetAugerTrueRange();


      NPrimaryElectrons = Track.GetNPrimaryElectrons();

      if(Track.nphe)nphel++;
      if(Track.naug)nauger++;

      
      if(PhotoelectronEnergy<5.8 && AugerPracticalRange>0){
      if(PhotoelectronPracticalRange>0)pheRange->Fill(PhotoelectronPracticalRange*10000);
      if(AugerPracticalRange>0)augRange->Fill(AugerPracticalRange*10000);
      if(PhotoelectronTrueRange>0)pheTrueRange->Fill(PhotoelectronTrueRange*10000);
      if(AugerTrueRange>0)augTrueRange->Fill(AugerTrueRange*10000);
      }

      Track.Drift();
      std::vector<std::pair<double,double> >  DifEl=Track.GetDiffElectronPosition();      
      std::vector<ADC> Digi= myDetector.mySampling(Gem.DiffusionofSecondaryElectrons(DifEl), Readout);
      //if(Digi.size()>10) // CS why event with less than 10 pixels are discarded?
      //cout << "CS DIGI SIZE: " << Digi.size() << endl;
      if(Digi.size()>1)
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
      else small_size++;
    } //end of the loop

  cout << "Number of photoel: " << nphel << endl;
  cout << "Number of Auger el.: " << nauger << endl;
  
  if(SetAstroSource){
    perc = (int)((float)GeneratedEvents/MaxGeneratedEvents*100);
      std::cout <<"  *** " << prc << " % done, ** N event conv. = " << N << ", tot. on the mirror area: " << GeneratedEvents <<std::endl;
  }
  else	  
    {
      perc = (int)((float)N/Number*100);
      std::cout<<" *** "<<perc<<" % done, ** N event conv. = "<<N<<", tot. on the mirror area: "<<GeneratedEvents<<std::endl;
    }
  cout << " Discarded events with size smaller than 10 pixels: " << small_size << endl;
  char Dummy[100];
  TString slash;
  TString DirName = "Trees";  // create directory Trees
  TString Command = "mkdir ";
  Command+=DirName;
  gSystem->Exec(Command);

#ifdef R__WIN32
  slash = "\\";
#else 
  slash = "/";
#endif
  
  TString name = GetName();   // create subdirectory
  DirName+=slash;
  DirName += name;
  Command = "mkdir ";
  Command+=DirName;
  gSystem->Exec(Command);
  DirName+=slash;
  DirName += File;
  TString F = DirName;
  sprintf(Dummy,"_%.1fkeV",PhotonEnergy); 
  F+=Dummy;
  F+=".root";
  TFile FileT(F,"RECREATE");
  tree->Write();
  FileT.Close();
  TString FD = DirName;
  FD+=Dummy;
  FD += ".txt";
  Dimension->WriteDimensionFile(FD);
  delete tree;
}

void     TExperiment::AnalyzeTree()
{
  TString Name = "../Trees/";
  Name+= GetName();
  char Dummy[100];
  if (Pitch==50)
    sprintf(Dummy,"/Eventstree_%.1fkeV",PhotonEnergy); 
  else if (Pitch==80)
    sprintf(Dummy,"/Pixmap8_%.1fkeV",PhotonEnergy); //modifica per Pixmap 80
  else std::cout<<"!!Pitch unknow!!  "<<std::endl;
  Name+=Dummy;
  std::cout<<"Analize file .... " << Name+".root" << std::endl;
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
  //double phi02 = fuficut->GetParameter(2);
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
  //double phi  = fufi->GetParameter(2);
  Eff = GetEfficiencyMixture()*100.;
  Effcut = Eff* (ModulationCutH->GetEntries())/TotalEvents;
  FitProbability = fufi->GetProb();
  mu   = 100.0 * B / (2.*A+B);  
  //  delete fufi;

  TF1 *fufif = (TF1*) ModulationFirstStepH->GetListOfFunctions()->FindObject("fufi");
  double Af    = fufif->GetParameter(0);
  double Bf    = fufif->GetParameter(1);
  //double phif  = fufif->GetParameter(2);
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

Double_t TExperiment::GetElectronRange(Double_t PheEnergy){
  return Mixture->GetElectronRange(PheEnergy);
}

Double_t TExperiment::GetPhotoelectricCrossSection(Double_t Energy){
  return Mixture->GetPhotoelectricCrossSection(Energy);
}

Double_t TExperiment::GetAbsorptionLenght(Double_t  Energy){
  return Mixture->GetAbsorptionLenght(Energy);
}

std::vector<TGraph*>  TExperiment::GetAbsProbGraph(Double_t  Energy){
  return Mixture->GetAbsProbGraph(Energy);
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
