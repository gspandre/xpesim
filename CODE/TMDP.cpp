#include "TMDP.h"
#include "TF1.h"

TMDP::TMDP()
{
  MirrorName="XEUS.txt";
  N_bins = 0;
  //Read file Window - Berillium window
  ifstream in;
  in.open("Window.txt");
  Float_t Energy[100], Trasp[100];
  Int_t nlines=0;
  Char_t DummyChar[20];
   for (Int_t j=0; j<3; j++)
     {
       in >> DummyChar;
     }
   while(1) 
     {
       in >> Energy[nlines] >> Trasp[nlines];
       if(!in.good())break;
       // std::cout<<Energy[nlines]<<" Tras =  "<<Trasp[nlines]<<std::endl; 
       nlines++;
     }
   Windowgr = new TGraph(nlines,Energy,Trasp);
   Windowgr->SetFillColor(19);
   Windowgr->SetMarkerColor(4);
   Windowgr->SetMarkerStyle(10);
   TString Title = "Window Transparency";
   Windowgr->GetXaxis()->SetTitle("Energy [keV]");
   //Windowgr->GetYaxis()->SetTitle("Trasnparency [cm^2]");
   Windowgr->GetYaxis()->SetTitle("Transparency ");
   Windowgr->GetYaxis()->SetRangeUser(0,1);
   in.close();
}

double         TMDP::MDP(TGraph *aeff, TF1 *flux, TGraph  *Effgas,double InEne,double FinEne,double T)
{
  double BinEn = 0.01;
  double A,F,M;
  double sqrtR = 0;
  double    MR = 0;
  double     R = 0;
  double mdp = 0;
  double Eff = 0;
  double En = InEne;
  while (En < FinEne) 
    { 
      A   = aeff->Eval(En);    //cm^2
      F   = flux->Eval(En);    //ph/kev/cm^2/s
      M   = Mug->Eval(En)/100.; //Mug al posto di mu
      // Eff = Windowgr->Eval(En) * Effgas->Eval(En)/100.;//good if Effgas (%) now (30/5/2003) with Transparency
      Eff = Effgas->Eval(En)/100.;//Eff ha gia' finestra 22/07/05
      R   +=      F * A * Eff * BinEn;  //ph/s
      MR  +=  M * F * A * Eff * BinEn;  //ph/s
      En += BinEn;
    }
  sqrtR = sqrt(R);
  mdp =  4.24*sqrtR/(MR*sqrt(T));
  //std::cout<<"   mdp = "<<mdp<<std::endl;
  double  MDPpoint=mdp;
  return MDPpoint;
}

TF1 *TMDP::Spectrum(double SpectralIndex, double FluxinCrab)
{
  TF1 *Spectrum = new TF1("Spectrum","[0]*x^(-[1])",0.1,50);
  Spectrum->SetTitle("Spectrum");
  Spectrum->SetMinimum(1e-3);
  std::cout<<"Constant =   "<<EvalFluxParameter(SpectralIndex,FluxinCrab)<<std::endl;
  Spectrum->SetParameters(EvalFluxParameter(SpectralIndex,FluxinCrab), SpectralIndex);
  return Spectrum;
}

TF1 *TMDP::Spectrum2(double SpectralIndex, double Constant)
{
  TF1 *Spectrum = new TF1("Spectrum","[0]*x^(-[1])",0.1,50);
  Spectrum->SetTitle("Spectrum");
  Spectrum->SetMinimum(1e-3);
  std::cout<<"Constant =   "<<Constant<<std::endl;
  Spectrum->SetParameters(Constant, SpectralIndex);
  Spectrum->GetXaxis()->SetTitle("energy[KeV]");
  Spectrum->GetYaxis()->SetTitle("flux [sec-1 KeV-1 cm^-2]");

  return Spectrum;
}

double TMDP::Rate(TF1 *flux, double InEne,double FinEne)// da il rate per la Crab
{
  TGraph *Effgas = EfficiencyGraph;
  TGraph *aeff = MirrorGraph(1.);

  double BinEn = 0.01;
  double En=InEne;
  double TotalEff = 0;
  double Rate=0;
  double SourceRate=0;
  double AfterMirror=0;
  while (En<FinEne)
    {
    
      TotalEff = aeff->Eval(En)* (Effgas->Eval(En))/100.;
      SourceRate += BinEn *flux->Eval(En);//fotoni cm-2 sec-1
      AfterMirror += BinEn * aeff->Eval(En)*flux->Eval(En);//fotoni sec-1 senza gas e berillio
      Rate += BinEn * TotalEff * flux->Eval(En);
      En+=BinEn;
    }
  std::cout<<"Source= "<<SourceRate<<"Source+Ottica="<<AfterMirror<<"Source+Ottica + berillio + Gas= "<<Rate<<std::endl;
  return Rate;//fotoni per secondo nell'intervallo InEne FinEne
}

/// new from Nic!

TH1D *TMDP::Spectrum(TF1 *sourceSpectrum, double emin,double emax)
{
  int N=1000;
  TH1D *obsspectrum = new TH1D("spectrum1","obsspectrum",N,emin,emax);
  for(int i=1;i<=N;i++)
    {
      double energy   = obsspectrum->GetBinLowEdge(i);
      double fi       = sourceSpectrum->Eval(energy); // fotoni cm-2 sec-1  kev-1
      obsspectrum->SetBinContent(i,fi);
    }
  obsspectrum->GetXaxis()->SetTitle("energy[KeV]");
  obsspectrum->GetYaxis()->SetTitle("flux [cm^-2 sec-1 KeV-1]");
  return obsspectrum;
}

TH1D *TMDP::Spectrum(TF1 *sourceSpectrum, TGraph *aeff, double emin,double emax)
{
  int N=1000;
  TH1D *obsspectrum = new TH1D("spectrum2","obsspectrum",N,emin,emax);
  for(int i=1;i<=N;i++)
    {
      double energy   = obsspectrum->GetBinLowEdge(i);
      double fi       = sourceSpectrum->Eval(energy); // fotoni cm-2 sec-1  kev-1
      double trasp    = aeff->Eval(energy);
      obsspectrum->SetBinContent(i,fi*trasp);
    }
  obsspectrum->GetXaxis()->SetTitle("energy[KeV]");
  obsspectrum->GetYaxis()->SetTitle("flux [sec-1 KeV-1]");
  return obsspectrum;
}

TH1D *TMDP::Spectrum(TF1 *sourceSpectrum, TGraph *aeff, TGraph *window, double emin,double emax)
{
  int N=1000;
  TH1D *obsspectrum = new TH1D("spectrum3","obsspectrum",N,emin,emax);
  for(int i=1;i<=N;i++)
    {
      double energy   = obsspectrum->GetBinLowEdge(i);
      double fi       = sourceSpectrum->Eval(energy); // fotoni cm-2 sec-1  kev-1
      double trasp    = aeff->Eval(energy)*window->Eval(energy);
      obsspectrum->SetBinContent(i,fi*trasp);
    }
  obsspectrum->GetXaxis()->SetTitle("energy[KeV]");
  obsspectrum->GetYaxis()->SetTitle("flux [sec-1 KeV-1]");
  return obsspectrum;
}

TH1D *TMDP::Spectrum(TF1 *sourceSpectrum, TGraph *aeff, TGraph *window, TGraph *mixEff, double emin,double emax)
{
  int N=1000;
  TH1D *obsspectrum = new TH1D("spectrum4","obsspectrum",N,emin,emax);
  for(int i=1;i<=N;i++)
    {
      double energy   = obsspectrum->GetBinLowEdge(i);
      double fi       = sourceSpectrum->Eval(energy); // fotoni cm-2 sec-1  kev-1
      double trasp    = aeff->Eval(energy)*window->Eval(energy)*mixEff->Eval(energy);
      obsspectrum->SetBinContent(i,fi*trasp);
    }
  obsspectrum->GetXaxis()->SetTitle("energy[KeV]");
  obsspectrum->GetYaxis()->SetTitle("flux [sec-1 KeV-1]");
  return obsspectrum;
}
//////////////////////////////////////////////////

//////////////////////////////////////////////////
TF1 *TMDP::CrabSpectrum()

{
  TF1 *CrabS =  Spectrum(2.1,1.);
  CrabS->SetLineColor(2);
  CrabS->SetLineWidth(2);
  return CrabS;
}


TF1   *TMDP::HerculesSpectrum()
{
  TFile *f= new TFile("HerSpectrum.root");
  Her     =(TF1*)f->Get("dndeabs");
  return Her;
}
void TMDP::SetMirror(TString MirrorN)
{
  MirrorName = MirrorN;
}

TGraph *TMDP::MirrorGraph(double  WindowEfficiency)
{
  ifstream in;
  //  in.open("axparea.dat");
 
  in.open(MirrorName);   
  Float_t Energy[2000], Area[2000];
  Int_t nlines=0;
  Char_t DummyChar[10];
  
  for (Int_t j=0; j<4; j++)
    in >> DummyChar;
  // cout << DummyChar << endl;
  
  while(1) {
    
    in >> Energy[nlines] >> Area[nlines];
    if(!in.good())break;
    
    if (MirrorName == "XEUS.txt")
      Area[nlines]*=WindowEfficiency * 10000.;
    if (MirrorName == "XEUSNEW.txt")
      {
	Energy[nlines] /= 1000.;
	Area[nlines]*=WindowEfficiency * 10000.;
      }
    if (MirrorName == "jetX.txt")
      {
	Area[nlines]*=WindowEfficiency * 5.;
      }
    else
      Area[nlines]*=WindowEfficiency;
    //std::cout <<nlines << "   " <<  Energy[nlines] << "   " << Area[nlines] << std::endl;
    nlines++;
  }

  printf(" found %d pointsn\n",nlines);
  
  TGraph *gr = new TGraph(nlines,Energy,Area);
  gr->GetXaxis()->SetRange(0,gr->GetXaxis()->FindBin(50));
  gr->SetFillColor(19);
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(21);
  TString Title = "Effective Area";
  char WindowName[150];
  sprintf(WindowName," * Window Efficiency [%.2f]",WindowEfficiency);
  if(WindowEfficiency != 1.0)
    Title+=WindowName;
  gr->SetTitle(Title);
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->CenterTitle();
  gr->GetXaxis()->SetTitle("Energy [keV]");
  gr->GetYaxis()->SetTitle("Area [cm^2]");
  in.close();
  return gr;
}


void     TMDP::DrawFlux(TString Source)
{
  TF1 *Spectrum;
  if (Source =="Crab")
    { 
      Spectrum = CrabSpectrum();
    }
  else 
    {
      Spectrum = HerculesSpectrum();
    }
  Spectrum->GetXaxis()->SetTitle("Energy [keV]");
  Spectrum->GetYaxis()->SetTitle(" ph/sec/keV/cm2 ");
  Spectrum->SetLineWidth(2);
  Spectrum->SetLineColor(4);
  Spectrum->Draw();
}

void TMDP::ReadDatafromGraph(TGraph *MuvsEn, TGraph *EffvsEn,int MIXID,double THICKNESS,double PRESSURE)
{
  //std::cout<<"Go1"<<std::endl;
  double *EnVG = MuvsEn->GetX();
  double *MuVG = MuvsEn->GetY();
  double *EffVG = EffvsEn->GetY();
  N_bins = MuvsEn->GetN();
  //std::cout<<N_bins<<std::endl;
  for (int i=0;i<N_bins;i++)
    {
      EnergyNVector.push_back(EnVG[i]);
      TotalEfficiencyVector.push_back(EffVG[i]);
      MuNVector.push_back(MuVG[i]);
    }
 
  MaximumEnergy = EnergyNVector[N_bins-1];
  MinimumEnergy = EnergyNVector[0];
  Thickness = THICKNESS;
  MixID = MIXID;
  Pressure = PRESSURE;
  FillEfficiencyGraph(); ///<<=========== Add Efficiency Berillium window
  FillMuGraph();
   
}


void TMDP::ReadData(TString FileName)
{
  
  std::ifstream FileIn;
  FileIn.open(FileName,std::ios::in);
  if(!FileIn.is_open())
    {
      std::cout<<"Error opening "<<FileName<<std::endl;
    }
   char dummy[20];
   FileIn >> dummy >> MixtureName;
   for (int i=0;i<7;i++)
    FileIn>>dummy;
   FileIn >> Thickness;
   for (int i=0;i<5;i++)
     FileIn>>dummy;
   FileIn >> MixID;
   FileIn>>dummy;
   FileIn >> N_bins;

   double  Energy;
   double  Mudata;
   double  Efficiency;
   double  CutEfficiency;
   double  TotalEfficiency;
   
   for (int i=0;i<4;i++)
     FileIn>>dummy;
   for (Int_t i=0; i<N_bins; i++)
     {  
       FileIn >> Energy>> Mudata >> Efficiency >> CutEfficiency;
       TotalEfficiency = Efficiency*CutEfficiency/100.;
       EnergyNVector.push_back(Energy);
       MuNVector.push_back(Mudata);
       TotalEfficiencyVector.push_back(TotalEfficiency);
      
	 }
   FileIn.close();
   MinimumEnergy = EnergyNVector[0];
   MaximumEnergy =EnergyNVector[N_bins-1];
   for(int i =0 ; i<N_bins; i++)
     std::cout<<" Energy = "<<EnergyNVector[i]<<" Mu =  "<<MuNVector[i]<<" TotalEff =  "<< TotalEfficiencyVector[i]<<std::endl;
   std::cout<<" Thickness = "<<Thickness<<" MixID ="<<MixID<<std::endl;
   //FillEfficiencyGraph();
   //FillMuGraph();
}

void   TMDP::FillEfficiencyGraph()
{
  const Int_t N = N_bins;
  if (N == 0) std::cout<<" NO DATA in charge "<<std::endl;
  Double_t EnVectorG[N];
  Double_t EfVectorG[N];
  for (int i =0;i<N;i++)
    {
      EnVectorG[i] = EnergyNVector[i];
      EfVectorG[i] = TotalEfficiencyVector[i]* Windowgr->Eval(EnergyNVector[i]);// finestra Berillio dal 22/07/05
    }
  EfficiencyGraph = new TGraph(N,EnVectorG,EfVectorG);
}

void TMDP::FillMuGraph()
{
  const Int_t N = N_bins;
  if (N == 0) std::cout<<" NO DATA in charge "<<std::endl;
  Double_t EnVectorG[N];
  Double_t MuVectorG[N];
  for (int i =0;i<N;i++)
    {
      EnVectorG[i] =  EnergyNVector[i];
      MuVectorG[i] = MuNVector[i];
    }
  Mug= new TGraph(N,EnVectorG,MuVectorG);
}


TGraph *TMDP::FillMDPGraph(TF1 *Source,double Time)
{
  Int_t Ncon = N_bins;
  if (Ncon == 0) std::cout<<" NO DATA in charge "<<std::endl;
  const Int_t N = (Int_t)(MaximumEnergy - MinimumEnergy - 3);//OK jetX 
  TGraph *Mirror;
 
  Mirror = MirrorGraph(1.);
 
   
  double Energypoint[N];
  double MDPpoint[N];

  //Mirror->Draw();
  //  EfficiencyGraph->Draw();
  for (int i=0; i< N; i++)
    {
      Energypoint[i] = 1.*i + MinimumEnergy + 1.;
      // std::cout<< Energypoint[i]<<std::endl;
      MDPpoint[i] = 100*MDP(Mirror,Source,EfficiencyGraph,Energypoint[i]-0.5,Energypoint[i]+0.5,Time);
     
    }
  TGraph *MDPGraph = new TGraph(N,Energypoint,MDPpoint);
  return MDPGraph;
}
 



void    TMDP::MDPCanvas()
{
  
  TCanvas *Ca = new TCanvas("Ca", "Ca",1200, 800);
  Ca->Divide(2,2);
  Ca->cd(1);
  TString Title1 = GetMixtureName();
  EfficiencyGraph->SetTitle(Title1);
  EfficiencyGraph->GetXaxis()->SetTitle("Energy(Kev)");
  EfficiencyGraph->GetYaxis()->SetTitle("Efficiency(%)");
  
  EfficiencyGraph->GetYaxis()->SetTitleColor(2);
  EfficiencyGraph->SetLineColor(2);
  EfficiencyGraph->GetYaxis()->SetLimits(0,100);
  EfficiencyGraph->SetMaximum(100);
  EfficiencyGraph->GetYaxis()->SetLabelSize(0.05);
  EfficiencyGraph->GetYaxis()->SetTitleSize(0.05);
  EfficiencyGraph->GetXaxis()->SetLabelSize(0.05);
  EfficiencyGraph->GetXaxis()->SetTitleSize(0.05);
  EfficiencyGraph->GetXaxis()->SetTitleOffset(0.7);
  EfficiencyGraph->Draw("alp");
  double xmax=EfficiencyGraph->GetXaxis()->GetXmax();
  double ymin=EfficiencyGraph->GetYaxis()->GetXmin();
  double ymax=EfficiencyGraph->GetYaxis()->GetXmax();
  //TGaxis *raxis = new TGaxis(xmax,ymin,xmax,ymax,0,ymax,510,"+L");
  // 
  //

  Mug->GetYaxis()->SetLimits(0,100);
 
  // Mug->SetTitle("Modulation Factor");
  //Mug->GetXaxis()->SetTitle("Energy(Kev)");
  //Mug->GetYaxis()->SetTitle("Modulation Factor(%)");
  Mug->Draw("lp");
  // raxis->Draw();  
  
  /*
    const Int_t N = N_bins;
    if (N == 0) std::cout<<" NO DATA in charge "<<std::endl;
    Double_t EnVectorG[N];
    Double_t MFactor[N];
    Double_t MaximumSF=0.0;
    for (int i =0;i<N;i++)
    {
    EnVectorG[i] =  EnergyNVector[i];
    //      MFactor[i] = (Mug->Eval(EnVectorG[i]))*sqrt((1/100.)*EfficiencyGraph->Eval(EnVectorG[i])) ;
    
    MFactor[i] = MuNVector[i]/100.0*sqrt((1./100.)*TotalEfficiencyVector[i]) ;
    std::cout<<EnVectorG[i]<<" "<<MFactor[i]<<std::endl;
    MaximumSF = TMath::Max(MaximumSF,MFactor[i]);
    }
    
    double ymax2 = 1.1*MaximumSF;
    double Scale = ymax/ymax2;
    
    for (int i =0;i<N;i++) MFactor[i]   = Scale *MFactor[i];
    
    TGraph *MuEp = new TGraph(N,EnVectorG,MFactor);
    MuEp->SetLineColor(3);
    
    MuEp->Draw("l");  
  */
  TGraph *MuEp=MuEpGraph();
  TGaxis *raxis = new TGaxis(xmax,ymin,xmax,ymax,0,0.24,510,"+L");//0.24 = ymax2
  raxis->SetTitle("#mu #sqrt{#epsilon}");
  raxis->SetTextColor(3);
  //
  raxis->SetLabelSize(0.05);
  raxis->SetTitleSize(0.05);
  //
  MuEp->Draw("l");  
  raxis->Draw(); 
  
   
  Ca->cd(2);
  TGraph *MirrorG;
  
  MirrorG = MirrorGraph(1.0);
  
  MirrorG->GetYaxis()->SetLabelSize(0.05);
  MirrorG->GetYaxis()->SetTitleSize(0.05);
  MirrorG->GetXaxis()->SetLabelSize(0.05);
  MirrorG->GetXaxis()->SetTitleSize(0.05);
  MirrorG->GetXaxis()->SetTitleOffset(0.7);
  
  MirrorG->Draw("acp");
  
  Ca->cd(3);
  TCanvas *SS = new TCanvas("SS","SS",1200,800);
  SS->cd();
  TGraph  *MDPher = FillMDPGraph(HerculesSpectrum());
  MDPher->SetMarkerStyle(24);
  MDPher->SetMarkerSize(0.7);
  MDPher->SetTitle("Hercules (black) Crab (red)");
  MDPher->GetXaxis()->SetTitle("Energy(Kev)");
  MDPher->GetYaxis()->SetTitle("MDP(%)");
  MDPher->GetYaxis()->SetLabelSize(0.05);
  MDPher->GetYaxis()->SetTitleSize(0.05);
  MDPher->GetXaxis()->SetLabelSize(0.05);
  MDPher->GetXaxis()->SetTitleSize(0.05);
  MDPher->GetXaxis()->SetTitleOffset(0.7);
  MDPher->Draw("acp");
  TGraph  *MDPcrab = FillMDPGraph(CrabSpectrum());
  MDPcrab->SetMarkerStyle(24);
  MDPcrab ->SetLineColor(2);
  MDPcrab->SetMarkerSize(0.7);
  //MDPcrab->SetTitle("Crab");
  MDPcrab->GetXaxis()->SetTitle("Energy(Kev)");
  MDPcrab->GetYaxis()->SetTitle("MDP(%)");
  MDPcrab->Draw("cp"); 
  TF1 *Crabs  = CrabSpectrum();  
  TF1 *Hers = HerculesSpectrum(); 
  double MDPH = MDP(MirrorG,Hers,EfficiencyGraph,MinimumEnergy,MaximumEnergy,3600.*24.);
  std::cout<<" total MDP =  "<<MDPH;
  char labelC[100];
  char labelH[100];
  sprintf(labelC," Crab MDP integrated (time = 1 day)  = %.3f %%",100*MDP(MirrorG,Crabs,EfficiencyGraph,MinimumEnergy,MaximumEnergy,3600.*24.));
  sprintf(labelH," Her-X1 MDP integrated (time = 1 day) = %.3f %%",100*MDP(MirrorG,Hers,EfficiencyGraph,MinimumEnergy,MaximumEnergy,3600.*24.));
  std::cout<<" MDP 5-12 =  "<<100*MDP(MirrorG,Crabs,EfficiencyGraph,5,12,3600.*24.)<<std::endl;
  std::cout<<" MDP 2-12 =  "<<100*MDP(MirrorG,Crabs,EfficiencyGraph,2,14,3600.*24.)<<std::endl;
  TText *ltC = new TText(3.5, 1.5 ,labelC);//new TText(3.5,0.01,labelC);
  TText *ltH = new TText(3.5, 2.0 ,labelH);//era 0.03
  ltC->SetTextColor(2);
  //ltH->SetTextColor();
  ltC->Draw();
  ltH->Draw();
  
  
  
  Ca->cd(4);
  gPad->SetLogy(); 
  gPad->SetLogx();
  Crabs->GetYaxis()->SetLabelSize(0.05);
  Crabs->GetYaxis()->SetTitleSize(0.05);
  Crabs->GetYaxis()->SetTitle("N photons / (cm^{2} sec keV)");
  Crabs->GetXaxis()->SetLabelSize(0.05);
  Crabs->GetXaxis()->SetTitleSize(0.05);
  Crabs->GetXaxis()->SetTitleOffset(0.07);
  Crabs->GetXaxis()->SetTitle("Energy(keV)");
  Crabs->Draw();
  Hers->Draw("same");
}
/*
  Ca->cd(4);
  TF1 *Crabs  = CrabSpectrum();
  Crabs->Draw();
  Her->Draw("same");
  Compounds[i]->GetCompoundName();
  Compounds[i]->GetFraction()
  GetPressure()
   GetCompounds()  std::vector<TCompound*>
 TPaveText *pt = new TPaveText(0.01,0.01,0.98,0.98);
  pt->AddText(MixtureName);
  long ObservationTime = 3600*24;
  //TString Muint = " 1h MDP  = ";
  char labelC[100];
  char labelH[100];
  sprintf(labelC," Crab   MDP (%d s)  = %.3f %%",ObservationTime,MDP.MDP(MEff,CrabSpectrum, Ef,E[0],E[N-1],ObservationTime)*100.);
  pt->AddText(labelC);
  sprintf(labelH," Her-X1 MDP (%d s) = %.3f %%",ObservationTime,MDP.MDP(MEff,HerSpectrum, Ef,E[0],E[N-1],ObservationTime)*100.);
  pt->AddText(labelH);
  pt->AddText("Gem Pitch =  50 #mum     Pixmap Pitch =  50 #mum");

 TText *ltC = new TText(3.5,0.01,labelC);
  TText *ltH = new TText(3.5,0.02,labelH);//era 0.03
  ltC->SetTextColor(2);
  ltH->SetTextColor(4);

  ltC->Draw();
  ltH->Draw();
  

*/
/*TGaxis *raxis = new TGaxis(xmax,ymin,xmax,ymax,0,100,510,"+L");
  raxis->SetTextColor(3);
  raxis->SetTitle("Mu(%)");
  M->Draw("cp");

  raxis->Draw();*/
void     TMDP::PlotDiffMdp()
{
  gPad->Clear();
  TGraph *MirrorG;
  MirrorG = MirrorGraph(1.0); 
  TGraph  *MDPher = FillMDPGraph(HerculesSpectrum());
  MDPher->SetMarkerStyle(24);
  MDPher->SetMarkerSize(0.7);
  
  // MDPher->SetTitle(Title1);
  MDPher->GetXaxis()->SetTitle("Energy(Kev)");
  MDPher->GetYaxis()->SetTitle("MDP(%)");
  MDPher->GetYaxis()->SetLabelSize(0.05);
  MDPher->GetYaxis()->SetTitleSize(0.05);
  MDPher->GetXaxis()->SetLabelSize(0.05);
  MDPher->GetXaxis()->SetTitleSize(0.05);
  MDPher->GetXaxis()->SetTitleOffset(0.7);
  MDPher->Draw("acp");
  TGraph  *MDPcrab = FillMDPGraph(CrabSpectrum());
  MDPcrab->SetMarkerStyle(24);
  MDPcrab ->SetLineColor(2);
  MDPcrab->SetMarkerSize(0.7);
  //MDPcrab->SetTitle("Crab");
  MDPcrab->GetXaxis()->SetTitle("Energy(Kev)");
  MDPcrab->GetYaxis()->SetTitle("MDP(%)");
  MDPcrab->Draw("cp"); 
  TF1 *Crabs  = CrabSpectrum();  
  TF1 *Hers = HerculesSpectrum(); 
  double MDPH = MDP(MirrorG,Hers,EfficiencyGraph,MinimumEnergy,MaximumEnergy,3600.*24.);
  std::cout<<" total MDP =  "<<MDPH;
  char labelC[100];
  char labelH[100];
  sprintf(labelC," Crab   MDP (1 day)  = %.3f %%",100*MDP(MirrorG,Crabs,EfficiencyGraph,MinimumEnergy,MaximumEnergy,3600.*24.));
  sprintf(labelH," Her-X1 MDP (1 day) = %.3f %%",100*MDP(MirrorG,Hers,EfficiencyGraph,MinimumEnergy,MaximumEnergy,3600.*24.));
  
  TText *ltC = new TText(3.5, 1.5 ,labelC);//new TText(3.5,0.01,labelC);
  TText *ltH = new TText(3.5, 2.0 ,labelH);//era 0.03
  ltC->SetTextColor(2);

  //ltH->SetTextColor();
  
  ltC->Draw();
  ltH->Draw();

}   

void TMDP::CompactPlot()
{
  gPad->Clear();
  TString Title1 = GetMixtureName();
  EfficiencyGraph->SetTitle(Title1);
  EfficiencyGraph->GetXaxis()->SetTitle("Energy(Kev)");
  EfficiencyGraph->GetYaxis()->SetTitle(" (%)");
  //EfficiencyGraph->GetYaxis()->SetTitle("Efficiency (red)  -   #mu (black)    (%)");
  /*
    TPaveLabel *p = new TPaveLabel(-1.4,77,-0.8,98,"#mu   (%)");
    p->SetTextAngle(90);
    p->SetTextSize(0.3);
    p->SetFillColor(10);
    p->SetLineColor(10);
    p->Draw();
    EfficiencyGraph->GetYaxis()->SetTitleOffset(0.9);
  */
  EfficiencyGraph->GetYaxis()->SetTitleSize(0.05);
  //EfficiencyGraph->GetYaxis()->CenterTitle();
  //EfficiencyGraph->GetYaxis()->SetTitleColor(2);
  EfficiencyGraph->SetLineColor(2);
  EfficiencyGraph->SetMarkerColor(2);
  EfficiencyGraph->SetMarkerStyle(8);
  EfficiencyGraph->SetMarkerSize(0.4);
  EfficiencyGraph->GetYaxis()->SetTitleOffset(0.9);
  
  EfficiencyGraph->GetYaxis()->SetLimits(0,100);
  EfficiencyGraph->SetMaximum(100);
  EfficiencyGraph->Draw("alp");
  double xmax=EfficiencyGraph->GetXaxis()->GetXmax();
  double ymin=EfficiencyGraph->GetYaxis()->GetXmin();
  double ymax=EfficiencyGraph->GetYaxis()->GetXmax();
  //TGaxis *raxis = new TGaxis(xmax,ymin,xmax,ymax,0,ymax,510,"+L");

  Mug->GetYaxis()->SetLimits(0,100);
  //Mug->SetTitle("Modulation Factor");
  //Mug->GetXaxis()->SetTitle("Energy(Kev)");
  //Mug->GetYaxis()->SetTitle("Modulation Factor(%)");
  Mug->SetMarkerStyle(23);
  Mug->SetMarkerSize(0.4);
  Mug->Draw("lp");
  
  const Int_t N = N_bins;
  if (N == 0) std::cout<<" NO DATA in charge "<<std::endl;
  Double_t EnVectorG[N];
  Double_t MFactor[N];
  Double_t MaximumMuEp=0.0;
  for (int i =0;i<N;i++)
    {
      EnVectorG[i] =  EnergyNVector[i];
      //  MFactor[i] = (1./100.)*(Mug->Eval(EnVectorG[i]))*sqrt((1/100.)*TotalEfficiencyVector[i]*Windowgr->Eval(EnVectorG[i]));
      MFactor[i] = MuNVector[i]/100.0*sqrt((1./100.)*TotalEfficiencyVector[i]*Windowgr->Eval(EnVectorG[i]));
      // std::cout<<EnVectorG[i]<<" "<<MFactor[i]<<" Eff ="<<EfficiencyGraph->Eval(EnVectorG[i])<<"vect * eff = "<<TotalEfficiencyVector[i]*Windowgr->Eval(EnVectorG[i])<<std::endl;
      MaximumMuEp = TMath::Max(MaximumMuEp,MFactor[i]);
    }
  
  double ymax2 = 0.24;
  double Scale = ymax/ymax2;
  
  for (int i =0;i<N;i++) MFactor[i]   = Scale *MFactor[i];
  
  TGraph *MuEp = new TGraph(N,EnVectorG,MFactor);
  MuEp->SetName("EpMu");
  FindMW(MuEp);
  MuEp->SetLineColor(3);
  MuEp->SetMarkerColor(3);
  MuEp->SetMarkerStyle(7);
    
  MuEp->Draw("pl");  
  
  TGaxis *raxis = new TGaxis(xmax,ymin,xmax,ymax,0,ymax2,510,"+L");  //(xmax,ymin,xmax,ymax,0,ymax2,510,"+L");
  raxis->SetTitle("#mu #sqrt{#epsilon}");
  raxis->SetLabelSize(0.05);
  raxis->SetTitleSize(0.05);
  raxis->CenterTitle();
  raxis->SetTextColor(3);
  raxis->Draw(); 

  TPaveLabel *pl = new TPaveLabel(9.81026,71.5736,14.6655,80.7741,"Modulation Factor","br");
  pl->SetFillColor(0);
  pl->SetLineColor(0);
  pl->SetTextSize(0.551724);
  pl->Draw();
  
  TPaveLabel *pl2 = new TPaveLabel(5.2106,5.90101,8.18122,11.9289,"Efficiency","br");
  pl2->SetFillColor(0);
  pl2->SetLineColor(0);
  pl2->SetTextColor(2);
  pl2->SetTextSize(0.99);
  pl2->Draw();
}


double     TMDP::EvalFluxParameter(double SpectralIndex,double FluxCrab,double MinEnergy, double MaxEnergy)
{
  TF1 IntFlux("IntegralFlux","(1./(2-[0]))*x^(-[0]+2)",0.1,20);
  IntFlux.SetParameter(0,SpectralIndex);
  return 13.4*FluxCrab/(IntFlux.Eval(MaxEnergy) -IntFlux.Eval(MinEnergy)); //7.88381179701537427e+00//12.5*FluxCrab/(IntFlux.Eval(MaxEnergy) -IntFlux.Eval(MinEnergy)); //7.88381179701537427e+00
}
// Aggiungo efficienza finestra
double TMDP::EvalMDP2_10kev(double SpectralIndex,double FluxinCrab,double Time)
{
  TF1 *SourceC = Spectrum(SpectralIndex,FluxinCrab);
  TGraph *MirrorG = MirrorGraph(1.);
  return MDP(MirrorG,SourceC,EfficiencyGraph, 2., 10., Time*3600);

}

void  TMDP::FluxMDPGraph(TString Name)
{

  const int N=10;
  double Flux[N];
  double MDPpointD2[N];
  double MDPpointH2[N];
  double MDPpointD3[N];
  double MDPpointH3[N];
  double MDPpointD0_5[N];
  double MDPpointH0_5[N];
  TCanvas *Canvas = new TCanvas("MDP", "MDP",1200, 800);
  gPad->SetLogy(); 
  gPad->SetLogx(); 
  for (int i=0;i<N;i++)//Flussi 1 giorno 1 ora
    {
      Flux[i] = 0.01*pow(2000./0.01,1.0*i/(N-1));
      MDPpointD0_5[i] = 100.*EvalMDP2_10kev(0.5,Flux[i]/1000.,24.);
      MDPpointD2[i] = 100.*EvalMDP2_10kev(2.01,Flux[i]/1000.,24.);
      MDPpointD3[i] = 100.*EvalMDP2_10kev(3.5,Flux[i]/1000.,24.);
      MDPpointH3[i] = 100.*EvalMDP2_10kev(3.5,Flux[i]/1000.,1.);
      MDPpointH0_5[i] = 100.*EvalMDP2_10kev(0.5,Flux[i]/1000.,1.);
      MDPpointH2[i] = 100.*EvalMDP2_10kev(2.01,Flux[i]/1000.,1.);
    }
  TGraph *GraD2= new TGraph(N,Flux,MDPpointD2);
  TGraph *GraD0_5= new TGraph(N,Flux,MDPpointD0_5);
  TGraph *GraD3= new TGraph(N,Flux,MDPpointD3);
  TGraph *GraH3= new TGraph(N,Flux,MDPpointH3);
  TGraph *GraH0_5= new TGraph(N,Flux,MDPpointH0_5);
  TGraph *GraH2= new TGraph(N,Flux,MDPpointH2);
  GraD2->GetXaxis()->SetRangeUser(0.01,1700);
  GraD2->SetMinimum(0.01);
  GraD2->SetMaximum(100.);
  
  GraD2->SetLineColor(2);
  GraH2->SetLineColor(2);
  GraD3->SetLineColor(3);
  GraH3->SetLineColor(3);
  GraD0_5->SetLineColor(4);
  GraH0_5->SetLineColor(4);
  GraH2->SetLineStyle(3);
  GraH3->SetLineStyle(3);
  GraH0_5->SetLineStyle(3);
  GraD2->GetXaxis()->SetTitle("Flux 2-10keV(mCrab)");
  GraD2->GetYaxis()->SetTitle("Minimum Detectable Polarization(%)");
  GraD2->GetXaxis()->CenterTitle();
  
  GraD2->SetTitle("");
  GraD2->Draw("alp"); 
  AddVariableSource(0.02/1000.,30./1000,1.2,2.7," BL Lac",0,23);
  AddVariableSource(0.02/1000.,5./1000,1.,2.01," FSRQ",0,46);
  AddVariableSource(0.3/1000.,0.75/1000,2.3,2.3," Circinus",1,9);
 
  GraH2->Draw("lp");
  GraD3->Draw("lp");
  GraH3->Draw("lp");
  GraD0_5->Draw("lp");
  GraH0_5->Draw("lp");
  TPaveText *Pave = new TPaveText(0.015,0.07,0.8,0.015);
  TText *t6 = Pave->AddText(Name);
  TText *t1 = Pave->AddText("Dotted line = observation time 1 hour");
  TText *t2 = Pave->AddText("Solid line = observation time 24 hour");
  TText *t3 = Pave->AddText("Spectral index = 0.5");
  TText *t4 = Pave->AddText("Spectral index = 2");
  TText *t5 = Pave->AddText("Spectral index = 3.5");
  t6->SetTextFont(62);
  t1->SetTextFont(72);
  t2->SetTextFont(72);
  t3->SetTextFont(72);
  t4->SetTextFont(72);
  t5->SetTextFont(72);
  t3->SetTextColor(4);
  t4->SetTextColor(2);
  t5->SetTextColor(3);
 
  TF1 *Hers = HerculesSpectrum(); 
  TGraph *MirrorX = MirrorGraph(1.0);
  double mdpherx1= 100*MDP(MirrorX,Hers,EfficiencyGraph,2,10,3600.*24.);
  TMarker *mher = new TMarker(63.3,mdpherx1,30);
  TText *ther = new TText(63.3,mdpherx1," Hercules X-1");
  ther->SetTextSize(0.035);
  ther->SetTextFont(12);
  Canvas->cd();  
  ther->Draw();
  mher->Draw();
  Pave->Draw();
  AddSource(1.,2.01," Crab");
  AddSource(2.7/1000.,2.01," MGC-6-30-15");
  AddSource(7.2/1000.,2.01," Centaurus A");
  AddSource(0.250,2.01," Vela X-1");
  TLine *Limit = new TLine(0.01,0.1,1700,0.1);
  Limit->SetLineColor(50);
  Limit->SetLineStyle(2);
  Limit->Draw();
  TText *Texlim = new TText(0.0266,0.1124,"Present limit of detector systematics");
  Texlim->SetTextFont(12);
  Texlim->SetTextSize(0.03);
  Texlim->SetTextColor(50);
  Texlim->Draw();
}
void TMDP::AddSource(double Flux,double Alpha,TString Name,double Time)
{
  TMarker *ma= new TMarker(1000*Flux,100.*EvalMDP2_10kev(Alpha,Flux,Time),30);
  TText *tex = new TText(1000*Flux,100*EvalMDP2_10kev(Alpha,Flux,Time),Name);
  ma->Draw();
  tex->SetTextSize(0.035);
  tex->SetTextFont(12);
  tex->Draw();
}
void TMDP::AddVariableSource(double Fluxmin,double Fluxmax,double Alphamin,double Alphamax,TString Name,int Marker,int color,double Time)
{
 
  double fluxvar[5] = {1000*Fluxmin,1000*Fluxmin,1000*Fluxmax,1000*Fluxmax,1000*Fluxmin};
  double mdpvar[5];
  mdpvar[0]= 100.*EvalMDP2_10kev(Alphamin,Fluxmin,Time);
  mdpvar[1]= 100.*EvalMDP2_10kev(Alphamax,Fluxmin,Time);
  mdpvar[2]= 100.*EvalMDP2_10kev(Alphamax,Fluxmax,Time);
  mdpvar[3]= 100.*EvalMDP2_10kev(Alphamin,Fluxmax,Time);
  mdpvar[4] = mdpvar[0];
  TText *tex=0; 
  if (Marker == 1)
    {
      AddSource(Fluxmin,Alphamin,Name,Time);
      AddSource(Fluxmax,Alphamax,Name,Time); 
    }
    /*  TLine *line = new TLine(1000*Fluxmin,100.*EvalMDP2_10kev(Alphamin,Fluxmin,Time),1000*Fluxmax,100.*EvalMDP2_10kev(Alphamax,Fluxmax,Time));
      line->Draw();*/
  else 
    {
      double xtext = pow(10.0,(log10(Fluxmax)+log10(Fluxmin))/2.);
      double ytext = pow(10.0,(log10(mdpvar[2])+log10(mdpvar[0]))/2.);

      tex = new TText(1000.0*xtext,ytext/3.,Name);
      tex->SetTextColor(color);
      tex->SetTextFont(12);
      // tex->SetTextSize(0.035);
       tex->SetTextSize(0.05);
  
  }
 
  TPolyLine *Poly = new TPolyLine(5, fluxvar,mdpvar);
  Poly->SetLineColor(color);
  Poly->SetFillColor(color);
  Poly->Draw("f");
  Poly->Draw();
  if (Marker == 0)tex->Draw();
}
TGraph *TMDP::MuEpGraph()
{
  const Int_t N = N_bins;
  if (N == 0) std::cout<<" NO DATA in charge "<<std::endl;
  Double_t EnVectorG[N];
  Double_t MFactor[N];
  for (int i =0;i<N;i++)
    {
      EnVectorG[i] =  EnergyNVector[i];
      //      MFactor[i] = (Mug->Eval(EnVectorG[i]))*sqrt((1/100.)*EfficiencyGraph->Eval(EnVectorG[i])) ;
      
      MFactor[i] = MuNVector[i]/100.0*sqrt((1./100.)*TotalEfficiencyVector[i]*Windowgr->Eval(EnergyNVector[i]));
      //std::cout<<EnVectorG[i]<<" "<<MFactor[i]<<std::endl;
    }
  
  for (int i =0;i<N;i++) MFactor[i]   =100/0.24*MFactor[i]; //100/0.24 *MFactor[i];
  
  TGraph *MuEpGr = new TGraph(N,EnVectorG,MFactor);
  MuEpGr->SetLineColor(3);
 
  return MuEpGr;
}
void TMDP::FindMW(TGraph *Graph)//Trova Max, x del Max, intervallo in X a meta' altezza
{
  GraphMax = 0.;
  GraphWsx = 0.;
  GraphWdx = 0.; 
  GraphMaxX = 0.;
  double *x=Graph->GetX();
  double *y=Graph->GetY();
  int N = Graph->GetN();
 
  for(int i=0;i<N;i++)
    {
      if (y[i]>GraphMax) 
	{
	  GraphMax = y[i]; 
	  GraphMaxX = x[i];
	}
    }
  double xvals =x[0];
  while (Graph->Eval(xvals)<GraphMax/2) 
    xvals+=0.2;
  GraphWsx = xvals;
  double xvald = x[N-1]- 0.1;
  while (Graph->Eval(xvald)<GraphMax/2)
    xvald-=0.1;
  GraphWdx = xvald;
}

double               TMDP::Rate(double InEne,double FinEne,TString Source)// da il rate per la Crab
{
  TGraph *Effgas = EfficiencyGraph;
  TGraph *aeff = MirrorGraph(1.);
  TF1 *flux=0;
  if (Source =="CRAB")
    {
      flux = CrabSpectrum();
      std::cout<<" CRAB "<<std::endl;
    }
  else if (Source =="CYGNUS")
    {
      flux = Spectrum(2.78,0.787); //Cygnus x-1
      std::cout<<" CYGNUS"<<std::endl;
    }
  else std::cout<<" Source ?  Select CRAB or CYGNUS"<<std::endl;
  double BinEn = 0.01;
  double En=InEne;
  double TotalEff = 0;
  double Rate=0;
  double SourceRate=0;
  double AfterMirror=0;
  while (En<FinEne)
    {
    
      TotalEff = aeff->Eval(En)* (Effgas->Eval(En))/100.;
      SourceRate += BinEn *flux->Eval(En);//fotoni cm-2 sec-1
      AfterMirror += BinEn * aeff->Eval(En)*flux->Eval(En);//fotoni sec-1 senza gas e berillio
      Rate += BinEn * TotalEff * flux->Eval(En);
      En+=BinEn;
    }
  std::cout<<"Source= "<<SourceRate<<"Source+Ottica="<<AfterMirror<<"Source+Ottica + berillio + Gas= "<<Rate<<std::endl;
  return Rate;//fotoni per secondo nell'intervallo InEne FinEne
}

double               TMDP::Rate(TH1D *spectrum, double InEne,double FinEne)// da il rate per la Crab
{

  double Rate=0;
  int    N    = spectrum->GetNbinsX();
  int    i    = 1;
  double energy= spectrum->GetBinLowEdge(1);

  while(energy<InEne)
    {
      i++;
      energy= spectrum->GetBinLowEdge(i);
    }
  while (energy<FinEne)
    {
      Rate+=spectrum->GetBinContent(i)*spectrum->GetBinWidth(i);
      i++;
      energy= spectrum->GetBinLowEdge(i);
    }
  std::cout<<" Rate [ "<<InEne<<" , "<<FinEne<<" ] = "<<Rate<<" ph/s "<<std::endl;
  return Rate;//fotoni per secondo nell'intervallo InEne FinEne
}
