#include <iostream>
#include <fstream>
#include "TRandom3.h"
#include "TString.h"
#include "TTreeAnalysis.h"
#include "TExperiment.h"
#include "ResultScan.cpp"
////////////////////////////////////////////
//TRandomGenerator *RandomGenerator = new TRandomGenerator();;
void An(TString Name)
{
 
  TTreeAnalysis A(Name,"EventsTree");
  A.TreeAnalysis();
  A.FillH();
  //  A.Save();
  //  B.Save();
  // A.MuCanvas();
  /* A.CheckCanvas();
  */  
  A.DeleteAll();
  // A.ImagingCanvas();
 }



void FitH(TH1 *His)

{
   TF1 *fun = new TF1("fun","[0]+[1]*(cos(x+[2]))^2",0,2*kPI);
   fun->SetParLimits(2,-kPI,kPI);
   fun->SetParLimits(0,0.,50000);
   fun->SetParLimits(1,0.,50000);
   His->Fit("fun","E","",0.,2*kPI);
   double A    = fun->GetParameter(0);
   double B    = fun->GetParameter(1);
   double phi0 = fun->GetParameter(2);
   double mu   = B / (2.*A+B);
   std::cout<<"mu His =  "<<mu<<std::endl;
   Float_t AError0 = fun->GetParError(0);
   Float_t BError0 = fun->GetParError(1);
   Float_t ModulationFactorError0 = 2.0*(A*BError0+B*AError0)/((B+2*A)*(B+2*A));
   std::cout<<"mu err =  "<<ModulationFactorError0<<std::endl;
}
void DoubleGausFit(TH1 *His)
{
  TF1 *g1 = new TF1("g1","[0]*TMath::Gaus(x,[1],[2]) +  [3]*TMath::Gaus(x,[4],[5])",2000,7000);
  g1->SetParLimits(1,3000,4000);
  g1->SetParLimits(2,10,5000);
  g1->SetParLimits(0,10,5000);
  g1->SetParLimits(3,10,5000);
  g1->SetParLimits(4,5000,6000);
  g1->SetParLimits(5,10,5000);
 
  // g1->Draw();
  His->Fit(g1);
  

}

/*
  void MDP(TString directory)
  {
  TMDP a;
  a.ReadData(directory);
  a.MDPCanvas();
  
  }
*/
void Scan()
{  
  int MixId;
  double Thickness;
  double Pressure;
  double Energy;
  double MuCut;
  double EffCut;
  double  mu;
  double FitProbability;
  double FitProbabilityCut;
  double Eff;
  double Mufirst;
  TRandom3 *rnd = new TRandom3();
  TExperiment Experiment(rnd);
  TTree *BigTree = new TTree("BigTree","BigTree");
  BigTree->Branch("MixID",&MixId,"MixId/I");
  BigTree->Branch("Thickness",&Thickness,"Thickness/D");
  BigTree->Branch("Pressure",&Pressure,"Pressure/D");
  BigTree->Branch("Energy",&Energy,"Energy/D");
  BigTree->Branch("MuCut",&MuCut,"MuCut/D");
  BigTree->Branch("Eff",&Eff,"Eff/D");
  BigTree->Branch("EffCut",&EffCut,"EffCut/D");
  BigTree->Branch("Mu",&mu,"mu/D");
  BigTree->Branch("FitProbability",&FitProbability,"FitProbability/D");
  BigTree->Branch("FitProbabilityCut",&FitProbabilityCut,"FitProbabilityCut/D");
  BigTree->Branch("Mufirst",&Mufirst,"Mufirst/D");
  for(MixId=2;MixId<7;MixId++)
    {
      Experiment.SetMixID(MixId);
      Pressure = 0.5;
      while (Pressure<2.)
	{
	  Thickness = 0.5;
	  Experiment.SetPressure(Pressure);
	  while (Thickness<2.5)
	    {
	      Experiment.SetThickness(Thickness);
	      Energy = 1.5;
	      while (Energy <15)
		{
		  Experiment.SetSource(30.,1.,Energy);
		  Experiment.Generalsetup();
		  // Experiment.EventsTree(10);
		  //Experiment.AnalyzeTree();
		  Experiment.ReadResult(MuCut,Eff,EffCut,mu,FitProbability,FitProbabilityCut,Mufirst);
		  Experiment.Delete();
		  BigTree->Fill();
		  if (Energy<10)Energy += 0.5;
		  else Energy += 1.;
		}
	      Thickness+=0.5;
	    }
	  Pressure+=0.5;
	}
    }
  TFile BigFile("BigFile.root","RECREATE");
  //BigFile->cd();
  BigTree->Write();
  BigFile.Close(); 
}
void ScanW()//Mandato il 5/5 Per analisi MixID low Energy (da 2 a 6)
{
  int MixId;
  double Thickness;
  double Pressure;
  double Energy;
  double MuCut;
  double EffCut;
  double  mu;
  double FitProbability;
  double FitProbabilityCut;
  double Eff;
  double Mufirst;
  TRandom3 *rnd = new TRandom3();
  TExperiment Experiment(rnd);
  TTree *BigTree = new TTree("BigTree","BigTree");//BigTree
  BigTree->Branch("MixID",&MixId,"MixId/I");
  BigTree->Branch("Thickness",&Thickness,"Thickness/D");
  BigTree->Branch("Pressure",&Pressure,"Pressure/D");
  BigTree->Branch("Energy",&Energy,"Energy/D");
  BigTree->Branch("MuCut",&MuCut,"MuCut/D");
  BigTree->Branch("Eff",&Eff,"Eff/D");
  BigTree->Branch("EffCut",&EffCut,"EffCut/D");
  BigTree->Branch("Mu",&mu,"mu/D");
  BigTree->Branch("FitProbability",&FitProbability,"FitProbability/D");
  BigTree->Branch("FitProbabilityCut",&FitProbabilityCut,"FitProbabilityCut/D");
  BigTree->Branch("Mufirst",&Mufirst,"Mufirst/D");
  int MixV[7]={2,3,4,5,6,7,20}; 
  for(int i=0; i<7; i++)//2 7
    {
      MixId = MixV[i];
      Experiment.SetMixID(MixId);
      Pressure = 0.5;//0.5
      while (Pressure<2.5)//2.
	{
	  Thickness = 0.5;//0.5
	  Experiment.SetPressure(Pressure);
	  while (Thickness<2.5)//2.5
	    {
	      Experiment.SetThickness(Thickness);
	      Energy = 1.5;
	      while (Energy <15)//15
		{
		  Experiment.SetSource(30.,1.,Energy);
		  Experiment.Generalsetup();
		  //Experiment.EventsTree(10000);
		  // Experiment.AnalyzeTree();
		  Experiment.ReadResult(MuCut,Eff,EffCut,mu,FitProbability,FitProbabilityCut,Mufirst);
		  Experiment.Delete();
		  BigTree->Fill();
		  if (Energy<10)Energy += 0.5;
		  else Energy += 1.;
		}
	      Thickness+=0.5;
	    }
	  Pressure+=0.5;
	}
    }
  TFile BigFile("BigFile15_11.root","RECREATE");
  //BigFile->cd();
  BigTree->Write();
  BigFile.Close(); 
}
void Scan13()//Mandato il 5/5 Per analisi MixID low Energy (da 2 a 6)
{
  int MixId;
  double Thickness;
  double Pressure;
  double Energy;
  double MuCut;
  double EffCut;
  double  mu;
  double FitProbability;
  double FitProbabilityCut;
  double Eff;
  double Mufirst;
  TRandom3 *rnd = new TRandom3();
  TExperiment Experiment(rnd);
  TTree *BigTree = new TTree("BigTree","BigTree");//BigTree
  BigTree->Branch("MixID",&MixId,"MixId/I");
  BigTree->Branch("Thickness",&Thickness,"Thickness/D");
  BigTree->Branch("Pressure",&Pressure,"Pressure/D");
  BigTree->Branch("Energy",&Energy,"Energy/D");
  BigTree->Branch("MuCut",&MuCut,"MuCut/D");
  BigTree->Branch("Eff",&Eff,"Eff/D");
  BigTree->Branch("EffCut",&EffCut,"EffCut/D");
  BigTree->Branch("Mu",&mu,"mu/D");
  BigTree->Branch("FitProbability",&FitProbability,"FitProbability/D");
  BigTree->Branch("FitProbabilityCut",&FitProbabilityCut,"FitProbabilityCut/D");
  BigTree->Branch("Mufirst",&Mufirst,"Mufirst/D");
  // const int NM = 10;
  const int NM = 12;
  int MixIdV[NM] = {2,3,4,5,6,7,20,10,11,0,1,24};
  for(int i=0;i<NM;i++)//2 7
    {
      Experiment.SetMixID(MixIdV[i]);
      MixId = MixIdV[i];
      if(i>6)
	{
	  Pressure=1;
	  Thickness=1;
	  Experiment.SetPressure(Pressure);
	  Experiment.SetThickness(Thickness);
	  Energy =1.5;
	  while (Energy <15)//15
	    {
	      Experiment.SetSource(30.,1.,Energy);
	      Experiment.Generalsetup();
	      //Experiment.EventsTree(10000);
	      //Experiment.AnalyzeTree();
	      Experiment.ReadResult(MuCut,Eff,EffCut,mu,FitProbability,FitProbabilityCut,Mufirst);
	      Experiment.Delete();
	      BigTree->Fill();
	      if (Energy<10)Energy += 0.5;
	      else Energy += 1.;
	    }
	}
      else
	{
	  Pressure = 0.5;//0.5
	  while (Pressure<2.5)//2.
	    {
	      Thickness = 0.5;//0.5
	      Experiment.SetPressure(Pressure);
	      while (Thickness<2.5)//2.5
		{
		  Experiment.SetThickness(Thickness);
		  Energy = 1.5;
		  while (Energy <15)//15
		    {
		      Experiment.SetSource(30.,1.,Energy);
		      Experiment.Generalsetup();
		      //Experiment.EventsTree(10000);
		      // Experiment.AnalyzeTree();
		      Experiment.ReadResult(MuCut,Eff,EffCut,mu,FitProbability,FitProbabilityCut,Mufirst);
		      Experiment.Delete();
		      BigTree->Fill();
		      if (Energy<10)Energy += 0.5;
		      else Energy += 1.;
		    }
		  Thickness+=0.5;
		}
	      Pressure+=0.5;
	    }
    
	}
      TFile BigFile("BigFile18_11.root","RECREATE");
      //BigFile->cd();
      BigTree->Write();
      BigFile.Close(); 
    }
}
void Scan20_5_05()// Da rifare
{
  int MixId;
  double Thickness;
  double Pressure;
  double Energy;
  double Energymax;
  TRandom3 *rnd = new TRandom3();
  TExperiment Experiment(rnd);
  
  for(MixId=8;MixId<10;MixId++)
 
  {
      Experiment.SetMixID(MixId);
      Pressure =1. ;//1
      while (Pressure<2.5)//
	{
	  Thickness = 0.5;//
	  Experiment.SetPressure(Pressure);
	  while (Thickness<2.5)
	    {
	      Experiment.SetThickness(Thickness);
	      Energy = 4.;//3
	      if (Pressure>1.5) Energymax=31;
	      else Energymax=26;
	      while (Energy <Energymax)//Energymax=30
		{
		  Experiment.SetSource(30.,1.,Energy);
		  Experiment.Generalsetup();
		  //Experiment.EventsTree(10000);
		   Experiment.AnalyzeTree();
		  //Experiment.ReadResult(MuCut,Eff,EffCut,mu,FitProbability,FitProbabilityCut);
		  Experiment.Delete();
		  
		  if (Energy<15)Energy += 1.;
		  else if (Energy>=15)Energy +=5.;
		  
		}
	      Thickness+=0.5;
	    }
	  Pressure+=0.5;
	}
    }
  
}

void Scan23_5_05()//Scan He con nuove drift
{
  int MixId;
  double Thickness;
  double Pressure;
  double Energy;
  TRandom3 *rnd = new TRandom3();
  TExperiment Experiment(rnd);
  
  for(MixId=24;MixId<25;MixId++)//=3 <7
    {
      Experiment.SetMixID(MixId);
      Pressure = 1;//0.5
      while (Pressure<1.5)//2.5
	{
	  Thickness = 1.;//0.5;
	  Experiment.SetPressure(Pressure);
	  while (Thickness<1.5)//2.5
	    {
	      Experiment.SetThickness(Thickness);
	      Energy = 1.5;
	      while (Energy <15)//15 CAMBIARE PER BIGTREE
		{
		  Experiment.SetSource(30.,1.,Energy);
		  Experiment.Generalsetup();
		  //Experiment.EventsTree(10000);
		  Experiment.AnalyzeTree();
		  //Experiment.ReadResult(MuCut,Eff,EffCut,mu,FitProbability,FitProbabilityCut);
		  Experiment.Delete();
		  
		  if (Energy<10)Energy += 0.5;
		  else Energy += 1.;
		}
	      Thickness+=0.5;
	    }
	  Pressure+=0.5;
	}
    }
  
}
void CheckTree()
{
  TRandom3 *rnd = new TRandom3();
  TExperiment Experiment(rnd);
  Experiment.SetMixID(7);
  Experiment.SetPressure(1);
  Experiment.SetThickness(1);
  Experiment.SetSource(30.,1.,6);
  Experiment.Generalsetup();
  Experiment.EventsTree(1000);
  //Experiment.AnalyzeTree();
  
}
void Scan3_10_05()
{
  
  double Thickness;
  double Pressure;
  double Energy;
  TRandom3 *rnd = new TRandom3();
  TExperiment Experiment(rnd);
  //const int NN = 10;//5;
  //int MixId[NN] = {2,3,4,5,6,10,11,7,0,1};0 1 10 11 
  const int NN = 7;
  int MixId[NN] = {20,7,2,3,4,5,6};
  for(int i=0;i<NN;i++)//=0 <NN
    {
      
      Experiment.SetMixID(MixId[i]);
      //Experiment.SetMixID(7);
      Pressure = 0.5;//0.5
      while (Pressure<2.5)//2.
	{
	  Thickness = 0.5;//0.5
	  Experiment.SetPressure(Pressure);
	  while (Thickness<2.5)//2.5
	    {
	      Experiment.SetThickness(Thickness);
	      Energy = 1.5;
	      while (Energy <15)//15
		{
		  Experiment.SetSource(30.,1.,Energy);
		  Experiment.Generalsetup();
		  //Experiment.EventsTree(10000);
		  Experiment.AnalyzeTree();
		  Experiment.Delete();
		  if (Energy<10)Energy += 0.5;
		  else Energy += 1.;
		}
	      Thickness+=0.5;
	    }
	  Pressure+=0.5;
	}
      
      
    }
}

void ScanMix2()//Modifico PixMap
{
  
  double Thickness;
  double Pressure;
  double Energy;
  TRandom3 *rnd = new TRandom3();
  TExperiment Experiment(rnd);
  //const int NN = 4;
  const int NN = 2;
  int MixId[NN]={2,3};
  
  for (int i=0; i<NN;i++)
    {
      Experiment.SetMixID(MixId[i]);
      
      Pressure = 1;//0.5
      Thickness = 1;//0.5;
      Experiment.SetPressure(Pressure);
      Experiment.SetThickness(Thickness);
      Energy = 1.5;
      while (Energy <15)//15 CAMBIARE PER BIGTREE
	{
	  Experiment.SetSource(30.,1.,Energy);
	  Experiment.Generalsetup();
	  Experiment.EventsTree(10000);
	  //Experiment.AnalyzeTree();
      //Experiment.ReadResult(MuCut,Eff,EffCut,mu,FitProbability,FitProbabilityCut);
	  Experiment.Delete();
	  
	  if (Energy<10)Energy += 0.5;
	  else Energy += 1.;
	}
    }
}



//Per fare tree leggibili da Pixie con 50 micron di passo occorre modificare il numero dei pixel riga colonna,riga 42 di TDimension.cpp
void Tree(int MixID,  double Pressure,double Thickness, double Energy,TString Name="PitchTree",double GPitch=50., double Pitch=50.)
{
  TRandom3 *rnd = new TRandom3();
  TExperiment Experiment(rnd);
  Experiment.SetMixID(MixID);
  Experiment.SetThickness(Thickness);
  Experiment.SetPressure(Pressure);
  Experiment.SetSource(30.,1.,Energy);
  Experiment.Generalsetup(GPitch,Pitch);
  Experiment.EventsTree(5000,Name);
  
}
void ScanPixmap()//gia' fatto 0.5 - 1 Ne DME
{
  int MixId;
  double Thickness;
  double Pressure;
  double Energy;
  TRandom3 *rnd = new TRandom3();
  TExperiment Experiment(rnd);
  const int NN = 7;//5;
  int MixIdV[NN] = {2,3,4,5,6,7,20};
  for(int i=0;i<NN;i++)//0 <7
    {
      MixId = MixIdV[i];//CF4
      Experiment.SetMixID(MixId);
      //Thickness = 0.5;
      Thickness = 1;
      Pressure = 1.0;
      Experiment.SetPressure(Pressure);
      // while (Thickness<2.)//2.
      {
	Experiment.SetThickness(Thickness);
	Energy = 1.5;
	while (Energy <15)//15 CAMBIARE PER BIGTREE
	  {
	    Experiment.SetSource(30.,1.,Energy);
	    Experiment.Generalsetup(50.,80.);
	    //Experiment.EventsTree(10000,"Pixmap8");
	    Experiment.AnalyzeTree();
	    //Experiment.ReadResult(MuCut,Eff,EffCut,mu,FitProbability,FitProbabilityCut);
	    Experiment.Delete();
	    
	    if (Energy<10)Energy += 0.5;
	    else Energy += 1.;
	  }
	//  Thickness+=0.5;
	//}
	
      }
    }
}


void ScanW8()//Mandato il 5/5 Per analisi MixID low Energy (da 2 a 6)
{
  int MixId;
  double Thickness;
  double Pressure;
  double Energy;
  double MuCut;
  double EffCut;
  double  mu;
  double FitProbability;
  double FitProbabilityCut;
  double Eff;
  double Mufirst;
  TRandom3 *rnd = new TRandom3();
  TExperiment Experiment(rnd);
  TTree *BigTree = new TTree("Big8","Big8");//BigTree
  BigTree->Branch("MixID",&MixId,"MixId/I");
  BigTree->Branch("Thickness",&Thickness,"Thickness/D");
  BigTree->Branch("Pressure",&Pressure,"Pressure/D");
  BigTree->Branch("Energy",&Energy,"Energy/D");
  BigTree->Branch("MuCut",&MuCut,"MuCut/D");
  BigTree->Branch("Eff",&Eff,"Eff/D");
  BigTree->Branch("EffCut",&EffCut,"EffCut/D");
  BigTree->Branch("Mu",&mu,"mu/D");
  BigTree->Branch("FitProbability",&FitProbability,"FitProbability/D");
  BigTree->Branch("FitProbabilityCut",&FitProbabilityCut,"FitProbabilityCut/D");
  BigTree->Branch("Mufirst",&Mufirst,"Mufirst/D");
  int MixV[7]={2,3,4,5,6,7,20}; 
  for(int i=0;i<7;i++)//0 7
    {
      MixId = MixV[i];
      Experiment.SetMixID(MixId);
      Pressure = 1;
      //if (i<5)
      {
	Thickness = 1;//0.5
	Experiment.SetPressure(Pressure);
	//while (Thickness<2.)//2.5
	{
	  Experiment.SetThickness(Thickness);
	  Energy = 1.5;
	  while (Energy <15)//15
	    {
	      Experiment.SetSource(30.,1.,Energy);
	      Experiment.Generalsetup(50,80);
	      //Experiment.EventsTree(10000);
	      // Experiment.AnalyzeTree();
	      Experiment.ReadResult(MuCut,Eff,EffCut,mu,FitProbability,FitProbabilityCut,Mufirst);
	      Experiment.Delete();
	      BigTree->Fill();
	      if (Energy<10)Energy += 0.5;
	      else Energy += 1.;
	    }
	  //Thickness+=0.5;
	}
	
      }
      /* else
	 {
	 
	 Thickness = 1.0;
	 Experiment.SetPressure(Pressure);
	 Experiment.SetThickness(Thickness);
	 Energy = 1.5;
	 while (Energy <15)//15
	 {
	 Experiment.SetSource(30.,1.,Energy);
	 Experiment.Generalsetup();
	 //Experiment.EventsTree(10000);
	 // Experiment.AnalyzeTree();
	 Experiment.ReadResult(MuCut,Eff,EffCut,mu,FitProbability,FitProbabilityCut,Mufirst);
	 Experiment.Delete();
	 BigTree->Fill();
	 if (Energy<10)Energy += 0.5;
	 else Energy += 1.;
	    }
	    }*/
    }
  TFile BigFile("BigFile8.root","RECREATE");
  //BigFile->cd();
  BigTree->Write();
  BigFile.Close(); 
}
void BasicAnS()
{
  const int N=7;
  double Energy[N] ={2,3,4,5,6,8,10} ;
  for (int i=0;i<N;i++)
    {
      
      TString Name;
      char Dummy[100];
      sprintf(Dummy,"../Pix/M2_8_%.1fkeV",Energy[i]);
      Name =Dummy;
      //Name+=".root";
      std::cout<<Name<<std::endl;
      TTreeAnalysis A(Name,"EventsTree");
      A.TreeAnalysis();
      A.FillH();
      A.DeleteAll();
      
    }
  
}

//////////////////////////////////////////////////
void MakeOneStep(int step=1)//Mandato il 5/5 Per analisi MixID low Energy (da 2 a 6)
{
  int MixId;
  double Thickness;
  double Pressure;
  double Energy;
  double MuCut;
  double EffCut;
  double  mu;
  double FitProbability;
  double FitProbabilityCut;
  double Eff;
  double Mufirst;
  TRandom3 *rnd = new TRandom3();
  TExperiment Experiment(rnd);
  TTree *BigTree = new TTree("BigTree","BigTree");//BigTree
  BigTree->Branch("MixID",&MixId,"MixId/I");
  BigTree->Branch("Thickness",&Thickness,"Thickness/D");
  BigTree->Branch("Pressure",&Pressure,"Pressure/D");
  BigTree->Branch("Energy",&Energy,"Energy/D");
  BigTree->Branch("MuCut",&MuCut,"MuCut/D");
  BigTree->Branch("Eff",&Eff,"Eff/D");
  BigTree->Branch("EffCut",&EffCut,"EffCut/D");
  BigTree->Branch("Mu",&mu,"mu/D");
  BigTree->Branch("FitProbability",&FitProbability,"FitProbability/D");
  BigTree->Branch("FitProbabilityCut",&FitProbabilityCut,"FitProbabilityCut/D");
  BigTree->Branch("Mufirst",&Mufirst,"Mufirst/D");
  static const int NumMixtures=7;
  int MixV[NumMixtures]={2,3,4,5,6,7,20}; 
  for(int i=0; i<NumMixtures; i++)//2 7
    {
      MixId = MixV[i];
      Experiment.SetMixID(MixId);
      Pressure = 0.5;//0.5
      while (Pressure<2.5)//2.
	{
	  Thickness = 0.5;//0.5
	  Experiment.SetPressure(Pressure);
	  while (Thickness<2.5)//2.5
	    {
	      Experiment.SetThickness(Thickness);
	      Energy = 1.5;
	      while (Energy <15)//15
		{
		  Experiment.SetSource(30.,1.,Energy);
		  Experiment.Generalsetup();
		  if(step==1)
		    {
		      Experiment.EventsTree(1000);
		    }
		  else if(step==2)
		    Experiment.AnalyzeTree();
		  else if(step==3)		  
		    {
		      Experiment.ReadResult(MuCut,Eff,EffCut,mu,FitProbability,FitProbabilityCut,Mufirst);
		    }
		  Experiment.Delete();
		  BigTree->Fill();

		  if (Energy<10)Energy += 0.5;
		  else Energy += 1.;
		}
	      Thickness+=0.5;
	    }
	  Pressure+=0.5;
	}
    }
  
  TFile BigFile("BigFile.root","RECREATE");
  //BigFile->cd();
  BigTree->Write();
  BigFile.Close(); 
}




//////////////////////////////////////////////////
void MakeMyScan(int step=1, int mix=5)//Gloria macro for MixID=5 low Energy (da 1.5 a 10)
{
  int MixId;
  double Thickness;
  double Pressure;
  double Energy;
  double MuCut;
  double EffCut;
  double  mu;
  double FitProbability;
  double FitProbabilityCut;
  double Eff;
  double Mufirst;
  TRandom3 *rnd = new TRandom3();
  TExperiment Experiment(rnd);
  TTree *BigTree = new TTree("BigTree","BigTree");//BigTree
  BigTree->Branch("MixID",&MixId,"MixId/I");
  BigTree->Branch("Thickness",&Thickness,"Thickness/D");
  BigTree->Branch("Pressure",&Pressure,"Pressure/D");
  BigTree->Branch("Energy",&Energy,"Energy/D");
  BigTree->Branch("MuCut",&MuCut,"MuCut/D");
  BigTree->Branch("Eff",&Eff,"Eff/D");
  BigTree->Branch("EffCut",&EffCut,"EffCut/D");
  BigTree->Branch("Mu",&mu,"mu/D");
  BigTree->Branch("FitProbability",&FitProbability,"FitProbability/D");
  BigTree->Branch("FitProbabilityCut",&FitProbabilityCut,"FitProbabilityCut/D");
  BigTree->Branch("Mufirst",&Mufirst,"Mufirst/D");
  MixId = mix;
  Pressure = 1.;
  Thickness = 1;
  Experiment.SetMixID(MixId);
  Experiment.SetPressure(Pressure);
  Experiment.SetThickness(Thickness);
  Energy = 1.5;
  while (Energy <11)
    {
      Experiment.SetSource(30.,1.,Energy);
      Experiment.Generalsetup();
      if(step==1)
	{
	  
	  std::cout << "Energy (keV): " << Energy << std::endl;
	  Experiment.EventsTree(50000);
	}
      else if(step==2)
	    Experiment.AnalyzeTree();
      else if(step==3)		  
	{
	  Experiment.ReadResult(MuCut,Eff,EffCut,mu,FitProbability,FitProbabilityCut,Mufirst);
	}
      Experiment.Delete();
      BigTree->Fill();
      
      if (Energy<8)Energy += 0.5;
      else Energy += 1.;
    }
  
  
  TFile BigFile("BigFile.root","RECREATE");
  //BigFile->cd();
  BigTree->Write();
  BigFile.Close(); 
}


void MakeAll(int mix=5)
{
  for(int i=1;i<4;i++)
    //MakeOneStep(i);
    MakeMyScan(i,mix);
    SetOption(mix,1,1,50,"BigFile.root");
    SetEnergy(6);
    MuCanvas(); // Result scan
    ImagingCanvas(); // Result scan

}


void MakeOnlyAnalysis(int mix=5)
{
  for(int i=2;i<4;i++)
    //MakeOneStep(i);
    MakeMyScan(i,mix);
    SetOption(mix,1,1,50,"BigFile.root");
    SetEnergy(6);
    MuCanvas(); // Result scan
    ImagingCanvas(); // Result scan

}

//////////////////////////////////////////////////
void DoOneRun(int MixId,double Pressure, double Thickness,double Energy)// Da rifare
{
  TRandom3 *rnd = new TRandom3();
  TExperiment Experiment(rnd);
  
  Experiment.SetMixID(MixId);
  Experiment.SetPressure(Pressure);
  Experiment.SetThickness(Thickness);
  Experiment.SetSource(25. , 0.1 , Energy);
  Experiment.Generalsetup();
  Experiment.EventsTree(500000);
  Experiment.Delete();  
}

void DoOneAnalysis(int MixId,double Pressure, double Thickness,double Energy)// Da rifare
{
  TRandom3 *rnd = new TRandom3();
  TExperiment Experiment(rnd);
  
  Experiment.SetMixID(MixId);
  Experiment.SetPressure(Pressure);
  Experiment.SetThickness(Thickness);
  Experiment.SetSource(0, 0 , Energy);
  Experiment.Generalsetup();
  Experiment.AnalyzeTree();
}

void DoBoth(int MixId,double Pressure, double Thickness,double Energy)// Da rifare
{
  DoOneRun(MixId, Pressure, Thickness, Energy);
  DoOneAnalysis(MixId, Pressure,Thickness,Energy);
}

void AstroScan(double TOBS=1.0,int N=1)
{  

  for (int I=1;I<=N;I++)
    {
      int MixId;
      double Thickness;
      double Pressure;
      /*  double Energy;
	  double MuCut;
	  double EffCut;
	  double  mu;
	  double FitProbability;
	  double FitProbabilityCut;
	  double Eff;
	  double Mufirst;
      */
      
      TRandom3 *rnd = new TRandom3();
      TString name="";
      if(N>1)
	{
	  char dummy[6];
	  sprintf(dummy,"%06d",I);
	  name+=dummy;
	}
      TExperiment Experiment(rnd,name);
      
      /*
	TTree *BigTree = new TTree("BigTree","BigTree");
    BigTree->Branch("MixID",&MixId,"MixId/I");
    BigTree->Branch("Thickness",&Thickness,"Thickness/D");
    BigTree->Branch("Pressure",&Pressure,"Pressure/D");
    BigTree->Branch("Energy",&Energy,"Energy/D");
    BigTree->Branch("MuCut",&MuCut,"MuCut/D");
    BigTree->Branch("Eff",&Eff,"Eff/D");
    BigTree->Branch("EffCut",&EffCut,"EffCut/D");
    BigTree->Branch("Mu",&mu,"mu/D");
    BigTree->Branch("FitProbability",&FitProbability,"FitProbability/D");
    BigTree->Branch("FitProbabilityCut",&FitProbabilityCut,"FitProbabilityCut/D");
    BigTree->Branch("Mufirst",&Mufirst,"Mufirst/D");
    //////////////////////////////////////////////////
    */
      
      Pressure = 1.0;
      Thickness = 1.0;
      MixId    = 5;
      double EMIN=2.0;
      double EMAX=10.0;
      TMDP *mdp = new TMDP();
      //  TF1* blazars = mdp->Spectrum(2.47,5e-4);
      TF1* blazars = mdp->Spectrum2(2.47,9.9e-3);

      mdp->SetMirror("XEUSNEW.txt");
      TGraph *myMirror = mdp->MirrorGraph(1.0);
      TGraph *myWindow = mdp->GetWindowEfficiency();
      
      TH1D *SourceSpectrum = mdp->Spectrum(blazars,EMIN,EMAX);
      TH1D *SourceSpectrumAeff = mdp->Spectrum(blazars,myMirror,EMIN,EMAX);
      TH1D *SourceSpectrumAeffWind = mdp->Spectrum(blazars,myMirror,myWindow,EMIN,EMAX);
      double R1 = mdp->Rate(SourceSpectrum,EMIN,EMAX);
      double R2 = mdp->Rate(SourceSpectrumAeff,EMIN,EMAX);
      double R3 = mdp->Rate(SourceSpectrumAeffWind,EMIN,EMAX);
      std::cout<<"Source Rate               : "<<R1<<" Tobs : "<<TOBS<<" Nevents: "<<R1*TOBS<<std::endl;
      std::cout<<"Source Rate Mirror        : "<<R2<<" Tobs : "<<TOBS<<" Nevents: "<<R2*TOBS<<std::endl;
      std::cout<<"Source Rate Mirror Window : "<<R3<<" Tobs : "<<TOBS<<" Nevents: "<<R3*TOBS<<std::endl;

      double R=R3;
      
      TCanvas *sourceCNV = new TCanvas("sourceCNV","sourceCNV");
      
      sourceCNV->Divide(2,2);
      sourceCNV->cd(1);
      gPad->SetLogx();
      gPad->SetLogy();
      myWindow->Draw("alp");
      
      sourceCNV->cd(2);
      gPad->SetLogx();
      gPad->SetLogy();
      myMirror->Draw("alp");
      
      sourceCNV->cd(3);
      gPad->SetLogx();
      gPad->SetLogy();
      blazars->Draw();
      
      sourceCNV->cd(4);
      gPad->SetLogx();
      gPad->SetLogy();
      SourceSpectrumAeffWind->Draw();
      
      Experiment.SetMixID(MixId);
      Experiment.SetPressure(Pressure);
      Experiment.SetThickness(Thickness);
      Experiment.SetSource(0.0,0.1,SourceSpectrumAeffWind);
      Experiment.Generalsetup();
      Experiment.SetMaxConvertedEvents(R*TOBS);
      //////////////////////////////////////////////////
      
      Experiment.EventsTree(100000000);
      
      name = "../Trees/";
      name = Experiment.GetName();
      Experiment.Delete();
      std::cout<<name<<std::endl;
    }
}

static const int N=4;
double E[N] = {2.,3., 4. ,10.};

void AstroSplitTree(TString fileName="fileName")
{
  TFile *f = new TFile(fileName,"OPEN");
  TString base_name=fileName.Replace(fileName.Sizeof()-6,6,"");
  std::cout<<base_name<<std::endl;
  TTree *t = (TTree*) f->Get("EventsTree");
  
  char Dummy[100];
  
  
  for(int i=0; i<N-1; i++)
    {
      double  E1 = E[i];
      double  E2 = E[i+1];

      TString OutputFilename=base_name;
      TString OriginalDimensionsFilename=OutputFilename;
      OriginalDimensionsFilename+=".txt";
      
      sprintf(Dummy,"_%dkeV-%dkeV", (int) E1,(int) E2);
      OutputFilename+=Dummy;      
      TString AnalysisFilename=OutputFilename;
      TString DimensionsFilename=OutputFilename;
      OutputFilename+=".root";
      DimensionsFilename+=".txt";
      gSystem->CopyFile(OriginalDimensionsFilename,DimensionsFilename);
      TFile fout(OutputFilename,"RECREATE");
      
      TString sel="";
      sprintf(Dummy,"Energy <= %f && Energy > %f", (float) E2,(float) E1);
      sel+=Dummy;
      
      TTree *selectiedTree = t->CopyTree(sel);
      std::cout<<selectiedTree->GetEntries()<<" between "<<E1<<" and "<<E2<<" keV... saved in  "<<OutputFilename<<std::endl;
      selectiedTree->Write();
      
      fout.Close();
      
      /*
	TTreeAnalysis A(AnalysisFilename,"EventsTree");
	A.TreeAnalysis();
	A.FillH();
	A.DeleteAll();      
      */
    }
}


void AstroProcessTree(TString fileName="fileName",  double ModulationDegree=1.0)
{
  TString base_name=fileName.Replace(fileName.Sizeof()-6,6,"");
  char Dummy[100];

  double Mu,MuErr,Phi,PhiError,CutEfficiency,FitProbabilityCut,Energy;
  double Mu2,MuErr2,Phi2,PhiError2,CutEfficiency2,FitProbabilityCut2;

  TTree *myTREE = new TTree("ResultsTree","results");
  myTREE->Branch("Energy",&Energy,"Energy/D");
  myTREE->Branch("Mu",&Mu,"Mu/D");
  myTREE->Branch("MuErr",&MuErr,"MuErr/D");
  myTREE->Branch("Phi",&Phi,"Phi/D");
  myTREE->Branch("PhiError",&PhiError,"PhiError/D");
  myTREE->Branch("CutEfficiency",&CutEfficiency,"CutEfficiency/D");
  myTREE->Branch("FitProbabilityCut",&FitProbabilityCut,"FitProbabilityCut/D");

  myTREE->Branch("Mu2",&Mu2,"Mu2/D");
  myTREE->Branch("MuErr2",&MuErr2,"MuErr2/D");
  myTREE->Branch("Phi2",&Phi2,"Phi2/D");
  myTREE->Branch("PhiError2",&PhiError2,"PhiError2/D");
  myTREE->Branch("CutEfficiency2",&CutEfficiency2,"CutEfficiency2/D");
  myTREE->Branch("FitProbabilityCut2",&FitProbabilityCut2,"FitProbabilityCut2/D");

  //  static const int N=4;
  //  double E[N] = {2.,3, 4 ,10.};
  //  static const int N=4;
  //double E[N] = {2.,3.,4.,10.};
  
  for(int i=0; i<N-1; i++)
    {
      double  E1 = E[i];
      double  E2 = E[i+1];
      Energy = (E1+E2)/2.;
      TString OutputFilename=base_name;
      sprintf(Dummy,"_%dkeV-%dkeV", (int) E1,(int) E2);
      OutputFilename+=Dummy;      

      TTreeAnalysis A(OutputFilename,"EventsTree");
      A.TreeAnalysis();
      A.FillH(ModulationDegree);
      	
      A.GetResults_SIMPLE(Mu,MuErr,Phi,PhiError,CutEfficiency,FitProbabilityCut);
      A.GetResults_CUT(Mu2,MuErr2,Phi2,PhiError2,CutEfficiency2,FitProbabilityCut2);
      
      A.DeleteAll();      
      myTREE->Fill();
    }
  TFile BigFile("AstroBigFile.root","RECREATE");
  myTREE->Write();
  BigFile.Close(); 
}

void AstroAnalysis(int step=1, double TOBS=100.0)
{
  TString Name = "../Trees/ScanHe_0.2-DME_0.8_1.0cm_1.0atm/EventstreeAstroScan.root";
  if (step==1)
    {
      AstroScan(TOBS);
    }
  else if (step==2)
    AstroSplitTree(Name);
  else if (step==3)
    AstroProcessTree(Name);

  
  

}
