#include "FastSimulator.h"

FastSimulator::FastSimulator(bool draw)
{
  Draw=draw;
  std::cout<<"Create fast simulator"<<std::endl;
  PI=TMath::Pi();
  PolarizationFraction=1.0;
  PolarizationAngle=0.0;

}

void FastSimulator::SetOption(int MixID,double Thickness,double Pressure,int PixmapS,TString BigName)
{
    
  if (PixmapS == 50)
    {
      bigfile = new TFile(BigName);
      bigtree = (TTree*) bigfile->Get("BigTree");
    }
  else
    {
      BigName = "BigFile8.root";
      bigfile = new TFile(BigName);//PixMap80
      bigtree = (TTree*) bigfile->Get("Big8");//PixMap80
    }
  std::cout<<"===>>> Opening .....  "<<BigName<<std::endl; 
  bigtree->SetDirectory(0);
  
  TString Option = "Thickness == ";
  Option+=Thickness;
  Option+=" && Pressure ==";
  Option+=Pressure;
  Option+=" && MixID == ";
  Option+=MixID;

  MIXID= MixID;
  THICKNESS = Thickness;
  PRESSURE = Pressure;
  PIXMAP =PixmapS;

  rnd = new TRandom3();
  Experiment = new TExperiment(rnd);
  Experiment->SetThickness(THICKNESS);
  Experiment->SetPressure(PRESSURE);
  Experiment->SetMixID(MIXID);
}

void FastSimulator::SetEnergy(double Energy,  TString TreeName)
{
  if (Energy !=0)
    {
      Experiment->SetSource(0,0,Energy);
      Experiment->Generalsetup();
      TString Name = Experiment->GetName();
      AnalysisName = "../Trees_scan/";
      AnalysisName += Name;
      char Dummy[100];
      if (PIXMAP == 50)
	{
	  sprintf(Dummy,"/Eventstree_%.1fkeV",Energy); 
	}
      else
	{
	  sprintf(Dummy,"/Pixmap8_%.1fkeV",Energy); 
	}
      AnalysisName +=Dummy ;
      DataName = AnalysisName;
      AnalysisName += "Analysis.root";
      std::cout<<AnalysisName<<std::endl;
      filetree  = new TFile(AnalysisName);
      treemix = (TTree*) filetree->Get("EventsTreef");
      treemix->SetDirectory(0);
    }
  else 
    {
      AnalysisName = TreeName;
      AnalysisName +="Analysis.root"; 
      DataName = TreeName;
      filetree  = new TFile( AnalysisName);
      treemix = (TTree*) filetree->Get("EventsTreef");
      treemix->SetDirectory(0);
    }
  std::cout<<AnalysisName<<" "<<treemix->GetEntries()<<std::endl;  
}


TH1F *FastSimulator::AngularResp()
{
  TString McTreeName = DataName;
  McTreeName +=".root";
  TTree *t = (TTree*)filetree->Get("EventsTreef");
  t->AddFriend("EventsTree",McTreeName);

  std::cout<<McTreeName<< " , " <<AnalysisName<<" N:  "<<t->GetEntries()<<std::endl;
  
  t->Draw("InitialPhi-Theta2>>h2");
  TH1F *h2 = (TH1F*)gDirectory->Get("h2"); 
  h2->GetXaxis()->SetTitle("Angle (rad)");
  h2->SetTitle("Reconstruction Second Step");
  return h2;
}

TF1 *FastSimulator::PhiDistribution()
{
  gDirectory->Delete("PhiDistribution");
  TF1* phiDistribution = new TF1("PhiDistribution","[0] + [1]*pow((cos(x-[2])), 2.0)", 0.0, 2.0*PI);
  phiDistribution->SetParameter(0, (1.0-PolarizationFraction));
  phiDistribution->SetParameter(1, 2.0*PolarizationFraction); // ?
  phiDistribution->SetParameter(2, PolarizationAngle);
  return phiDistribution;
}

void FastSimulator::ProcessEvents(int N)
{
  
  TF1 *phiDistribution = PhiDistribution();
  TH1F *DetMC = AngularResp();

  gDirectory->Delete("Incoming");
  gDirectory->Delete("Reconstructed");

  TH1D *Incoming      = new TH1D("Incoming","Incoming",100,0.0,2.0*PI);
  TH1D *Reconstructed = new TH1D("Reconstructed","Reconstructed",100,0.0,2.0*PI);
  
  for (int i=0 ; i < N ; i++)
   { 
     double Phi_IN  = phiDistribution->GetRandom();
     double Phi_REC = Phi_IN + DetMC->GetRandom();
    
     if (Phi_REC>2.*PI)
       Phi_REC-=2*PI;
     if (Phi_REC<0)
       Phi_REC+=2*PI;
     
     Incoming->Fill(Phi_IN);
     Reconstructed->Fill(Phi_REC);
     
     if (i % 1000000==0)
       { 
	 std::cout<<" ********** Nevent = "<<i<<std::endl;
       }
   }
  
  if(Draw)
    {
      a= new TCanvas("D","D",1200,800);
      a->Divide(1,3);
      a->cd(1);
      Incoming->Draw();
      
      a->cd(2);
      DetMC->Draw();
      
      a->cd(3);
      Reconstructed->Draw();
    }
  gDirectory->Delete("fufi");
  TF1 *fufi1 = new TF1("fufi","[0]+[1]*(cos(x-[2]))^2",0,2*PI);
  fufi1->SetParLimits(2,-PI,PI);
  fufi1->SetParLimits(0,0.,N*1.0);
  fufi1->SetParLimits(1,0.,N*1.0);
  Reconstructed->Fit("fufi","E","",0.,2*PI);  
  double A    = fufi1->GetParameter(0);
  double B    = fufi1->GetParameter(1);
  ModulationFactor   = 100.0 * B / (2.*A+B);
  Float_t DA = fufi1->GetParError(0);
  Float_t DB = fufi1->GetParError(1);

  AngleError = fufi1->GetParError(2)*180.0/TMath::Pi();
  ModulationFactorError = 200.0*(A*DB+B*DA)/((B+2*A)*(B+2*A));
  std::cout<<" Mu = "<<ModulationFactor<<"+-"<<ModulationFactorError <<" Angle Error (gradi)="<<AngleError<<std::endl;
  
} 

void FastSimulator::GetResults(double &Mu,double &MuE, double &AngleErr)
{
  Mu=ModulationFactor;
  MuE=ModulationFactorError;
  AngleErr=AngleError;
}
