#include <iostream>
#include <fstream>
#include "TRandom3.h"
#include "TString.h"
#include "TMDP.h"
#include "TTree.h"
#include "TFile.h"
#include "TPaveText.h"
#include "TText.h"
#include  "TExperiment.h"
#include "TTreeAnalysis.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TGasMixture.h"
#include "TGraphErrors.h"
#include "TMath.h"

TFile *file;
TFile *filetree;
TFile *filedata;
TTree *tree;
TTree *treemix;
TTree *treedata;
TString AnalysisName,Name,DataName;

TString Option;
int MIXID;
double THICKNESS;
double PRESSURE;
TExperiment *Experiment;
TRandom3 *rnd;
//Settaggio miscela, Spessore, Pressione
void SetOption(int MixID,double Thickness,double Pressure,int PixmapS = 50,TString BigName="BigFile.root")
{
  if (PixmapS == 50)
    {
      file = new TFile(BigName);
      tree = (TTree*) file->Get("BigTree");
    }
  else
    {
      BigName = "BigFile8.root";
      file = new TFile(BigName);//PixMap80
      tree = (TTree*) file->Get("Big8");//PixMap80
    }
  std::cout<<"===>>> Opening .....  "<<BigName<<std::endl; 
  tree->SetDirectory(0);
  //  file->Close();

  Option = "Thickness == ";
  Option+=Thickness;
  Option+=" && Pressure ==";
  Option+=Pressure;
  Option+=" && MixID == ";
  Option+=MixID;
  MIXID= MixID;
  THICKNESS = Thickness;
  PRESSURE = Pressure;
  std::cout<<Option<<std::endl;
  rnd= new TRandom3();
  Experiment = new TExperiment(rnd);
  Experiment->SetThickness(THICKNESS);
  Experiment->SetPressure(PRESSURE);
  Experiment->SetMixID(MIXID);
   
}

//Fornisce il grafico del Fattore di Modulazione (con tagli):Energia
TGraph *GetMuvsEn()
{
  tree->Draw("MuCut:Energy",Option);
  TGraph *MuVEn = (TGraph*)gPad->FindObject("Graph");
  MuVEn->SetTitle("MuvEn");
  return (TGraph*)MuVEn->Clone();
}

////Fornisce il grafico dell'efficienza (con tagli):Energia
TGraph *GetEffvsEn()
{
  tree->Draw("EffCut:Energy",Option);
  TGraph *EffVEn = (TGraph*)gPad->FindObject("Graph");
  EffVEn->SetTitle("EffvsEn");
  return (TGraph*)EffVEn->Clone();
}

////Fornisce il grafico della MDP crab e Her X-1:Energia
void ConfrontoAnalisi(TString File1, TString File2, double Pressure=1, double Thickness=1)
{
  TCanvas *con=new TCanvas("con","con",1200,800);
  con->Divide(3,3);
  int MixV[7]={2,3,4,5,6,7,20};
  for (int i=0; i<7;i++)
  {
    con->cd(i+1);
    gPad->SetGridx();
    gPad->SetGridy();
    SetOption(MixV[i],Thickness,Pressure,50, File1);
    TGraph *g1 =GetMuvsEn();
    g1->Draw("alp");
    SetOption(MixV[i],Thickness,Pressure,50, File2);
    TGraph *g2 =GetMuvsEn();
    g2->SetMarkerColor(2);
    g2->Draw("lp");
  }
}

void DiffMdp()
{
  TMDP a;
  TGraph *MuvsEn = GetMuvsEn();
  TGraph *EffvsEn = GetEffvsEn();
  a.ReadDatafromGraph(MuvsEn,EffvsEn,2,1.0,1.0);
  Experiment->Generalsetup();
  TString Ti= Experiment->GetName();
  TCanvas *c = new TCanvas("c",Ti,1200,800);
  a.PlotDiffMdp();
  std::cout<<" Ti = "<<Ti<<" Exp->" <<Experiment->GetName() <<std::endl;
  TText *TT = new TText(3.5,2.5,Ti);
  TT->Draw();
}

//Fornisce il grafico del Fattore di Modulazione (con tagli):Energia, Efficienza(con finestra Berillio), Fattore di Qualita' (Mu * sqrt(effic.))
void CCompactPlot() 
{
 
  TMDP a;
  TGraph *MuvsEn = GetMuvsEn();
  TGraph *EffvsEn = GetEffvsEn();
  a.ReadDatafromGraph(MuvsEn,EffvsEn,5,1.0,1.0);
  a.CompactPlot();

}

void FactorGraph(int Mixid)//
{
  TCanvas *prov=new TCanvas("PRO","PRO",1200,800);
  //
  const int NP = 4;
  const int NT = 4;
  double Pressure[NP]={0.5,1.,1.5,2.};
  double Thickness[NT]={0.5,1.,1.5,2.};
  const int N = NP*NT;
  double Max[N];
  double MaxX[N];
  double Widthsx[N];
  double Widthdx[N];
  TString SettingName[N];
  double TotalMDP[N];
  double Index[N];
  double Mu2keV[N];
  double Mu6keV[N];
  for (int i=0; i<NP;i++)
    {
      for (int j=0; j<NT;j++)
	{
	   TMDP a;
	  SetOption(Mixid,Thickness[j],Pressure[i]);
	  TGraph *MuvsEn = GetMuvsEn();
	  Mu2keV[(NP)*i+j]=MuvsEn->GetY()[1];
	  Mu6keV[(NP)*i+j]=MuvsEn->GetY()[9];
	  std::cout<<"Energy Mu2 ="<<MuvsEn->GetX()[1]<<std::endl;
	  std::cout<<"Energy Mu6 ="<<MuvsEn->GetX()[9]<<std::endl;
	  TGraph *EffvsEn = GetEffvsEn();
	  a.ReadDatafromGraph(MuvsEn,EffvsEn,2,1.0,1.0);
	  a.FindMW(a.MuEpGraph());
	  Max[(NP)*i+j]=a.GetGraphMax()*0.0024;
	  MaxX[(NP)*i+j]=a.GetGraphMaxX();
	  Widthsx[(NP)*i+j]=a.GetGraphWsx();
	  Widthdx[(NP)*i+j]=a.GetGraphWdx();
	  std::cout<<"Max =  "<<Max[(NP)*i+j]<<" Energy max ="<<MaxX[(NP)*i+j]<<" En sx = "<<MaxX[(NP)*i+j] - Widthsx[(NP)*i+j]<<" En dx = "<<MaxX[(NP)*i+j] + Widthdx[(NP)*i+j]<<std::endl;
	  Experiment->Generalsetup();
	  SettingName[(NP)*i+j] = Experiment->GetName();
	  TotalMDP[(NP)*i+j]=100.*a.EvalMDP2_10kev(2.01,1.,24);//24 ora
	  Index[(NP)*i+j] =(NP)*i+j; 
	  std::cout<<Index[(NP)*i+j]<<"Totalmdp =   "<<TotalMDP[(NP)*i+j]<<std::endl;
	}
    }
   TCanvas *can = new TCanvas("Factor","Factor",1200,800);
   TH2D *graph1 = new TH2D("FQ","FQ",N,1,N+1,10,1.5,5);
   TAxis *ax1 = graph1->GetXaxis();
   for(int i=1;i<=N;i++)
     {
       ax1->SetBinLabel(i,SettingName[i-1]);
       graph1->SetBinContent(i,MaxX[i-1],Max[i-1]);
     }
   can->SetTopMargin(0.0417311);
   can->SetBottomMargin(0.451314);
   graph1->LabelsOption("v","X");
   graph1->SetStats(false);
   graph1->Draw("LEGO");
   TCanvas *ca2 = new TCanvas("TotalMDP", "TotalMDP",6,20,1199,799);
   TH1D *graph2 = new TH1D("","",N,1,N+1);
   TAxis *ax = graph2->GetXaxis();
   for(int i=1;i<=N;i++)
     {
       ax->SetBinLabel(i,SettingName[i-1]);
       graph2->SetBinContent(i,TotalMDP[i-1]);
     }
   //  gStyle->SetOptStat(0);
   ca2->SetRightMargin(0.0327056);
   ca2->SetTopMargin(0.0417311);
   ca2->SetBottomMargin(0.451314);
   graph2->LabelsOption("v","X");
   graph2->SetMaximum(0.06);
   // graph2->GetYaxis()->SetNdivisions(520);
   
   graph2->SetMinimum(0.02);
   graph2->SetStats(false);
   graph2->Draw();
   for (int i = 0; i<N ; i++)
   
     {
       
       std::cout<<SettingName[i]<<std::endl;;
       std::cout<<"Max = "<<Max[i]<<" Energy max ="<<MaxX[i]<<" Total MDP 2-10="<<TotalMDP[i]<<" Mu2keV= "<< Mu2keV[i]<<" Mu6= "<<Mu6keV[i] <<std::endl;
     }
   //for (int i = 0; i<N ; i++)
   std::cout<<" \\hline "<<std::endl;
   std::cout<<"  Pr  Sp&Max FdQ (Energia)& Mu $2$ keV&  Mu $6$ keV & MDP (2-10) \\\\ "<<std::endl;
   for (int i=0; i<16; i++)
     {
       char Text[400];
       sprintf(Text,"%.1f %.1f & %.3f ( %.1f )  & %.2f  & %.2f & %.3f \\\\  ",Pressure[i/4],Thickness[i- (i/4)*4],Max[i],MaxX[i],Mu2keV[i],Mu6keV[i],TotalMDP[i]);
       //std::cout<<" & "<<Max[i]<<"("<<MaxX[i]<<") &"<<Mu2keV[i]<<" & "<<Mu6keV[i]<<" & "<<TotalMDP[i]<<" \\" <<"\\"<<std::endl;
       std::cout<<" \\hline "<<std::endl;
       TString TextS = "";
       TextS+=Text;
       std::cout<<TextS<<std::endl;
     }
}
void FactorGraph1_1(int Mixid,int PixmapS = 50)
{
  TCanvas *prov=new TCanvas("PRO","PRO",1200,800);
  //
  
  double Max;
  double MaxX;
  double Widthsx;
  double Widthdx;
  TString SettingName;
  double TotalMDP;
  TMDP a;
  SetOption(Mixid,1,1,PixmapS);
  TGraph *MuvsEn = GetMuvsEn();
  TGraph *EffvsEn = GetEffvsEn();
  a.ReadDatafromGraph(MuvsEn,EffvsEn,2,1.0,1.0);
  a.FindMW(a.MuEpGraph());
  Max = a.GetGraphMax()*0.0024;
  MaxX = a.GetGraphMaxX();
  Widthsx = a.GetGraphWsx();
  Widthdx = a.GetGraphWdx();
  std::cout<<"Max =  "<<Max<<" Energy max ="<<MaxX<<" En sx = "<<Widthsx<<" En dx = "<<Widthdx<<std::endl;
  Experiment->Generalsetup();
  SettingName = Experiment->GetName();
  TotalMDP=100.*a.EvalMDP2_10kev(2.01,1.,24);
  std::cout<<" TotalMDP 2-10 =  "<<TotalMDP<<std::endl;

}

//CompactPlot per alcune miscele
void VarMixture(double Th=1,double Pr =1)
{
  //FactorGraph(Mixid);
  const int NT = 7;
  int Mixid[NT]={2,3,4,5,6,7,20};//7,10,11,0,1
  double MaxMuEp[NT];
  TCanvas *MdpA = new TCanvas("MixtureScan","MixtureScan",1200,800);
  MdpA->SetTitle(Name);
  MdpA->Divide(3,3);
  // TMDP a;
  for (int j=0; j<NT;j++)
    {
      SetOption(Mixid[j],Th,Pr);
      MdpA->cd(j+1);
     
      //CCompactPlot();
      TMDP a;
      TGraph *MuvsEn = GetMuvsEn();
      TGraph *EffvsEn = GetEffvsEn();
      a.ReadDatafromGraph(MuvsEn,EffvsEn,2,1.0,1.0);
      a.CompactPlot();
      MaxMuEp[j]=a.GetMaximumMuEp();
      
      MdpA->cd(j+1);
      Experiment->Generalsetup();
      Name=Experiment->GetName();
      TText *TT = new TText(4.5,90,Name);
	
      TT->Draw();
    }

}
//CompactPlot per una miscelaq variando Spessore
void MixtureAnT(int Mixid)
{
  FactorGraph(Mixid);
  const int NT = 4;
  double Thickness[NT]={0.5,1.,1.5,2.};//{0.5,1.,1.5,2.}
  TCanvas *MdpA = new TCanvas("MixtureScan","MixtureScan",1200,800);
  MdpA->SetTitle(Name);
  MdpA->Divide(2,2);
  for (int j=0; j<NT;j++)
    {
      SetOption(Mixid,Thickness[j],1);
      MdpA->cd(j+1);
     
      CCompactPlot();
      MdpA->cd(j+1);
      Experiment->Generalsetup();
      Name=Experiment->GetName();
      TText *TT = new TText(4.5,90,Name);
	
      TT->Draw();
    }

}
//CompactPlot per una miscelaq variando Pressione
void MixtureAnP(int Mixid)
{
  FactorGraph(Mixid);
  const int NP = 4;
  double Pressure[NP]={0.5,1.,1.5,2.};//{1.,1.5,2.}
  TCanvas *MdpA = new TCanvas("MixtureScan","MixtureScan",1200,800);
  MdpA->SetTitle(Name);
  MdpA->Divide(2,2);
  for (int j=0; j<NP;j++)
    {
      SetOption(Mixid,1,Pressure[j]);
      MdpA->cd(j+1);
     
      CCompactPlot();
      MdpA->cd(j+1);
      Experiment->Generalsetup();
      Name=Experiment->GetName();
      TText *TT = new TText(4.5,90,Name);
	
      TT->Draw();
    }

}

//CompactPlot per una miscela variando Spessore e pressione
void MixtureAn(int Mixid)
{
  //SetOption(Mixid,1.,1.);
  //Name=Experiment->GetName();
  //FactorGraph(Mixid);
  const int NP = 4;
  const int NT = 4;
  double Pressure[NP]={0.5,1.,1.5,2.};//{1.,1.5,2.}
  double Thickness[NT]={0.5,1.,1.5,2.};//{0.5,1.,1.5,2.}
  TCanvas *MdpA = new TCanvas("MixtureScan","MixtureScan",1200,800);
  MdpA->SetTitle(Name);
  MdpA->Divide(NP,NT);
  for (int i=0; i<NP;i++)
    {
      for (int j=0; j<NT;j++)
	{
	  SetOption(Mixid,Thickness[j],Pressure[i]);
	  MdpA->cd((NP)*i+j+1);
	  std::cout<<"i= "<<i<<" j = "<<j<<" cd =  "<<(NP)*i+j+1<<std::endl;
	  /* char text[100];
	  sprintf(text,"Thickness  =%.1f Pressure = %.1f",Thickness[j],Pressure[i]);
	  TText *T = new TText(6.5,5,text);
	  */
	  CCompactPlot();
	  MdpA->cd((NP)*i+j+1);
	  Experiment->Generalsetup();
	  Name=Experiment->GetName();
	  TText *TT = new TText(4.5,90,Name);
	  //T->Draw();
	  TT->Draw();
	}
    }
}
//Scelta la miscela(con SetOption), fornisce grafico CompactPlot+XEUS+FlussoCrab e Her+Mdp differenziale
void Mdp(TString Mir = "XEUS.txt")
{
  TMDP a;
  TGraph *MuvsEn = GetMuvsEn();
  TGraph *EffvsEn = GetEffvsEn();
  a.ReadDatafromGraph(MuvsEn,EffvsEn,2,1.0,1.0);
  a.SetMirror(Mir);
  a.MDPCanvas();
 
  Experiment->Generalsetup();
  // a.FluxMDPGraph(Experiment->GetName());
  
}
//Grafico MDP per varie sorgenti
void SourceMdp(TString Mir = "XEUS.txt")
{
  
  TMDP a;
  a.SetMirror(Mir);
  TCanvas *pro2 = new TCanvas("","",1200,800);
  TGraph *MuvsEn = GetMuvsEn();
  TGraph *EffvsEn = GetEffvsEn();
  a.ReadDatafromGraph(MuvsEn,EffvsEn,2,1.0,1.0);
  a.MDPCanvas();
  Experiment->Generalsetup();
  a.FluxMDPGraph(Experiment->GetName());

}
//Settaggio Energia
void SetEnergy(double Energy,  TString TreeName = "Eventstree",int  PixmapS = 50)
{
  
  if (Energy !=0)
    {
      Experiment->SetSource(0,0,Energy);
      Experiment->Generalsetup();
      Name = Experiment->GetName();
      AnalysisName = "../Trees/";
      AnalysisName += Name;
      char Dummy[100];
      if (PixmapS == 50)
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
}
// Settati Miscela (con SetOption) e Energia (SetEnergy) visualizza canvas Fattore di Modulazione
void  MuCanvas()
{
  TCanvas *Mod = new TCanvas("Mod",Name,1200,800);
  Mod->Divide(3,2); 
  Mod->cd(2);
  TH1D *ThirdMomHi = (TH1D*) filetree->Get("ThirdMomHi");
  ThirdMomHi->Draw();

  Mod->cd(5);
  TH1D *SecondRatioHi =(TH1D*) filetree->Get("SecondRatioHi");
  SecondRatioHi->Draw();
  double SecondRatioMean = SecondRatioHi->GetMean();
  std::cout<<" Second Ratio Mean = "<<SecondRatioMean<<std::endl;
  
  Mod->cd(1);
  TH1D *ModulationCutH = (TH1D*) filetree->Get("ModulationCut");
  TF1 *fuficut = (TF1*) ModulationCutH->GetListOfFunctions()->FindObject("fufi");
  double A2    = fuficut->GetParameter(0);
  double B2    = fuficut->GetParameter(1);
  double phi02 = fuficut->GetParameter(2);
  ModulationCutH->Draw();
  double mucut   = 100.0 * B2 / (2.*A2+B2); 
  

  
  Mod->cd(3); 
  TH1D *ModulationSimpleH = (TH1D*)filetree->Get("ModulationSimple");
  ModulationSimpleH->Draw("");
  TF1 *fufi = (TF1*) ModulationSimpleH->GetListOfFunctions()->FindObject("fufi");
  double A    = fufi->GetParameter(0);
  double B    = fufi->GetParameter(1);
  double phi = fufi->GetParameter(2);
  
  double mu2   = 100.*B / (2.*A+B);
  std::cout<<"mu2 Second= "<<mu2<<std::endl;
  
  Mod->cd(6);
  TString Mufit2 = "#mu2 =";
  Mufit2 +=mu2;
  TString Err2 = "#mu2 error = ";
  TString Stat = "N statistic = ";
  Stat +=  ModulationSimpleH->GetEntries();
  Float_t AError2 = fufi->GetParError(0);
  Float_t BError2 = fufi->GetParError(1);
  Float_t ModulationFactorError2 = 2.0*(A*BError2+B*AError2)/((B+2*A)*(B+2*A));
  Float_t AngleError2 = fufi->GetParError(2)*180.0/TMath::Pi();
  Err2+= ModulationFactorError2;
  TPaveText *pt2=new TPaveText(0.01,0.01,0.98,0.98);
  TText *tt1 =pt2->AddText(Mufit2);
  TText *tt2 =pt2->AddText(Err2);
  TText *tt3 = pt2->AddText(Stat);
  pt2->Draw();
  
  Mod->cd(4);
  TString MufitC = "#mu (Cut) =";
  MufitC +=mucut;
  TString ErrC = "#mu (Cut) error = ";
  TString StatC = "N statistic = ";
  StatC +=   ModulationCutH->GetEntries();
  Float_t AError = fuficut->GetParError(0);
  Float_t BError = fuficut->GetParError(1);
  Float_t ModulationFactorError = 2.0*(A2*BError+B2*AError)/((B2+2*A2)*(B2+2*A2));
  Float_t AngleError = fuficut->GetParError(2)*180.0/TMath::Pi();
  ErrC+= ModulationFactorError;
  TPaveText *ptc=new TPaveText(0.01,0.01,0.98,0.98);
  TText *tt1c =ptc->AddText(MufitC);
  TText *tt2c =ptc->AddText(ErrC);
  TText *tt3c = ptc->AddText(StatC);
  ptc->Draw();
  
  TCanvas *Co= new TCanvas("Mod Factor","Mod Factor",1200,800);
  Co->Divide(2,2);
  Co->cd(1);
  TCanvas *SS = new TCanvas("SS","SS",1200,800);
  SS->cd();
  ModulationCutH->Draw();
  std::cout<<"Mu= "<<mucut<<"Delta Mu= "<<ModulationFactorError<<std::endl;
  Co->cd(2);
  ModulationSimpleH->Draw();
  Co->cd(4);
  TH1D *ModulationFirst = (TH1D*)filetree->Get("ModulationFirstStep");
  ModulationFirst->Draw();
  TF1 *fufifirst = (TF1*) ModulationFirst->GetListOfFunctions()->FindObject("fufi");
  double Afirst    = fufifirst->GetParameter(0);
  double Bfirst    = fufifirst->GetParameter(1);
  double phifirst = fufifirst->GetParameter(2);
  double mufirst   = 100.0 * Bfirst / (2.*Afirst+Bfirst);
  TString Mufirst = "#mu first =";
  Mufirst +=mufirst;
  Co->cd(3);
  TPaveText *pt4=new TPaveText(0.01,0.01,0.98,0.98);
  TText *ttcut =pt4->AddText(MufitC);
  TText *tt2uncut =pt4->AddText(Mufit2);
  TText *ttfirst = pt4->AddText(Mufirst);
  pt4->SetTextFont(12);
  pt4->SetTextSize(0.07);
  pt4->Draw();
}
void MUC()
{
  TCanvas *Mod2 = new TCanvas("Mod2",Name,1200,800);
  TH1D *ModulationCutH = (TH1D*) filetree->Get("ModulationCut");
  TF1 *fuficut = (TF1*) ModulationCutH->GetListOfFunctions()->FindObject("fufi");
  double A2    = fuficut->GetParameter(0);
  double B2    = fuficut->GetParameter(1);
  double phi02 = fuficut->GetParameter(2);
  ModulationCutH->Draw();
  double mucut   = 100.0 * B2 / (2.*A2+B2); 
   

}
// Settati Miscela (con SetOption) e Energia (SetEnergy) visualizza canvas Impatto ricostruito, Mu primo passo, Pulse Height
void CheckCanvas()
{
  TCanvas *ChCanvas = new TCanvas("ChCanvas",Name,1200,800);
  ChCanvas->Divide(3,2); 
  
  ChCanvas->cd(1);
  treemix->Draw(" (ImpactY-RealImpactY) * (BaricenterX-RealImpactX)/(sqrt(pow(RealImpactX- BaricenterX,2)+pow(RealImpactY- BaricenterY,2))) - (ImpactX-RealImpactX) * (BaricenterY-RealImpactY)/(sqrt(pow(RealImpactX- BaricenterX,2)+pow(RealImpactY- BaricenterY,2))): (ImpactX-RealImpactX) * (BaricenterX-RealImpactX) /(sqrt(pow(RealImpactX- BaricenterX,2)+pow(RealImpactY- BaricenterY,2))) + (ImpactY-RealImpactY) * (BaricenterY-RealImpactY)/(sqrt(pow(RealImpactX- BaricenterX,2)+pow(RealImpactY- BaricenterY,2))) >> Hbid(1000,-300,800,1000,-400,400)");
  
  ChCanvas->cd(2);
  //t->Draw("Theta2:InitialPhi");
  treemix->Draw("  (ImpactX-RealImpactX) * (BaricenterX-RealImpactX) /(sqrt(pow(RealImpactX- BaricenterX,2)+pow(RealImpactY- BaricenterY,2))) + (ImpactY-RealImpactY) * (BaricenterY-RealImpactY)/(sqrt(pow(RealImpactX- BaricenterX,2)+pow(RealImpactY- BaricenterY,2)))");
  ChCanvas->cd(3);
  treemix->Draw("BaricenterY:BaricenterX");

  ChCanvas->cd(4);
  treemix->Draw("Pulse");
  
  ChCanvas->cd(5);
  treemix->Draw("AnClusterdim");
  
  ChCanvas->cd(6);
  treemix->Draw("ThetaI>>FirstStep");//,"Auger==0");
  TH1D *NoA = (TH1D*)gDirectory->FindObject("FirstStep");
  TF1 *fu = new TF1("fu","[0]+[1]*(cos(x+[2]))^2",0,2*kPI);
  fu->SetParLimits(2,-kPI,kPI);
  fu->SetParLimits(0,0.,50000);
  fu->SetParLimits(1,0.,50000);
  NoA->Fit("fu","E","",0.,2*kPI);
  double AA    = fu->GetParameter(0);
  double BB    = fu->GetParameter(1);
  double phip = fu->GetParameter(2);
  double mum   = BB / (2.*AA+BB);
  Float_t AAError = fu->GetParError(0);
  Float_t BBError = fu->GetParError(1);
  Float_t ModulationFactorErrorI = 2.0*(AA*BBError+BB*AAError)/((BB+2*AA)*(BB+2*AA));
  std::cout<<" Mu = "<<mum<<" Error = "<<ModulationFactorErrorI<<std::endl;
}
void DrawString(TString Variable)
{
  treemix->Draw(Variable); 
}
// Settati Miscela (con SetOption) e Energia (SetEnergy) visualizza canvas capacita' di Imaging
void  ImagingCanvas()
{
  TCanvas *ImagC = new TCanvas("ImagC",Name,1200,800);
  ImagC->Divide(2,2);
  
  ImagC->cd(1);
  treemix->Draw("ImpactY:ImpactX");//>>hi2(1000,100,12000,1000,2000,10000)");
 
  ImagC->cd(2);
  treemix->Draw("ImpactX");
 
  ImagC->cd(3);
  treemix->Draw("BaricenterY:BaricenterX>>hb2(1000,100,12000,100,2000,10000)");
 
  ImagC->cd(4);
  treemix->Draw("RealImpactY:RealImpactX>>hre(1000,100,12000,1000,2000,10000)");

  TCanvas *Res = new TCanvas("Res",Name,1200,800);
  Res->Divide(2,2);
  
  Res->cd(1);
  treemix->Draw("ImpactY");
  
  Res->cd(2);
  treemix->Draw("RealImpactY >>hreal");
  
  Res->cd(3);
  treemix->Draw("ImpactX-RealImpactX>>hres");
    
  
  Res->cd(4);
  treemix->Draw("ImpactY-RealImpactY");
  
  TCanvas *ImpX =new TCanvas("Im","Im",1200,800);
  ImpX->Divide(2,2);
  ImpX->cd(1);
  treemix->Draw("abs(ImpactX-RealImpactX)");
  ImpX->cd(2);
  treemix->Draw("abs(BaricenterX-RealImpactX)");
}


//Settati Miscela (con SetOption) e Energia (SetEnergy) visualizza traccia numero n
void ClusterviewInCharge(int n)
{
  Name+=".root";
  TTreeAnalysis A(DataName,"EventsTree"); 
  A.ClusterAnalysis(n);
}

void RawClusterview(int n)
{
  Name+=".root";
  TTreeAnalysis A(DataName,"EventsTree"); 
  std::vector<ADC> Vect= A.ReadCluster(n);
  // TH1D *H = new TH1D("Raw","Raw",20000,0,22000);
  //const int N=Vect.size();
  //const int N=22000;
  const int N=Vect.size();
  double Channel[N];
  double DigiCharge[N];
  
  /*std::vector<ADC>::iterator pos;
  for (int i=10000; i<12000; i++)
    {
      for (pos= Vect.begin(); pos!=Vect.end();++pos)
	{
	  if (i == (*pos).Channel )
	    {
	      Channel[i] = (*pos).Channel;
	      DigiCharge[i] = (*pos).DigiCharge;
	      std::cout<<Channel[i]<<"   "<<DigiCharge[i]<<std::endl;
	    }
	  else 
	    {
	      Channel[i]=i;
	      DigiCharge[i] =1400;
	    }
	}
	}*/
   for (int i=0; i<N; i++)
     {
       Channel[i] = Vect[i].Channel;
       DigiCharge[i] = Vect[i].DigiCharge-1400;
       std::cout<<Channel[i]<<"   "<<DigiCharge[i]<<std::endl;
     }
  TGraph *g = new TGraph(N,Channel,DigiCharge);
  
  g->Draw("alp");
}
void ClusterComp(TString InPut1,TString InPut2,int n)
{
  TCanvas *Confr = new TCanvas("Confr","Confr",1200,800);
  //Confr->Divide(4,4);
  Confr->Divide(2,2);
  //  for (int i = 1; i<9; i++)
  for (int i = 1; i<3; i++)
    {
      SetEnergy(0,InPut1);
      TTreeAnalysis A(DataName,"EventsTree"); 
      A.ClusterViewComp(n+i); 
      // TVirtualPad *p1 = gPad;
      TCanvas *p1 = (TCanvas*)gPad;
      Confr->cd(i);
      p1->DrawClonePad();
    }
  //for (int i = 9; i<17; i++)
  for (int i = 3; i<5; i++)
    {
      SetEnergy(0,InPut2);
      TTreeAnalysis A(DataName,"EventsTree"); 
      A.ClusterViewComp(i-2+n); 
      // TVirtualPad *p1 = gPad;
      TCanvas *p1 = (TCanvas*)gPad;
      Confr->cd(i);
      p1->DrawClonePad();
    }
  /* TCanvas *Confr2 = new TCanvas("Confr2","Confr2",1200,800);
  SetEnergy(0,InPut1);
  TTreeAnalysis A(DataName,"EventsTree"); 
  A.ClusterViewComp(n+2);
  TCanvas *p2 = (TCanvas*)gPad;
  p2->DrawClonePad();
  TCanvas *Confr3 = new TCanvas("Confr3","Confr3",1200,800);
  SetEnergy(0,InPut2);
  TTreeAnalysis B(DataName,"EventsTree"); 
  B.ClusterViewComp(n+1);
  TCanvas *p3 = (TCanvas*)gPad;
  p3->DrawClonePad();*/
}
/*{
  Name+=".root";
  TTreeAnalysis A(DataName,"EventsTree"); 
  TCanvas *Confr = new TCanvas("Confr","Confr",1200,800);
  Confr->Divide(4,2);
  A.ClusterViewComp(n); 
  for (int i = 1; i<9; i++)
    {
      A.ClusterViewComp(i); 
      // TVirtualPad *p1 = gPad;
      TCanvas *p1 = (TCanvas*)gPad;
      Confr->cd(i);
      p1->DrawClonePad();
    }
 
    }*/

TH1D  *GetClusterH()
{
  
  TH1D *G = (TH1D*) filetree->Get("ClusterDim");
  // G->Draw();
  return G;
}
// restituisce il Rate per Crab(vedi TMDP.cpp), per una miscela selezionata tra 2 estremi di energia
double   BinRate(int MixID,double Thickness,double Pressure,double Initialbin,double Finalbin,TString MirrorN="XEUSNEW.txt")
{
  TMDP aa;
  aa.SetMirror(MirrorN);
  SetOption(MixID,Thickness,Pressure);
  TGraph *MuvsEn = GetMuvsEn();
  TGraph *EffvsEn = GetEffvsEn();
  aa.ReadDatafromGraph(MuvsEn,EffvsEn,2,1.0,1.0);
  //double TotalRate = aa.Rate(2.,10.);
  return aa.Rate(Initialbin,Finalbin,"CYGNUS");
}

void NichDistribution(double R,double Tau, TString Title)
{
  TCanvas *ca = new TCanvas("ca","ca",1200,800);
  ca->Divide(2,2);
  TH1D *His = new TH1D("","",1000,0,1000);
  TH1I *Poiss = new TH1I("","",20,0,20);
  double Time =0;
  int Counts =0;
  for (int i=0;i<10000;i++)
    {
      double csi = rnd->Uniform();
      double dt = -1./R * log(1. - csi);
      His->Fill(dt);
      Time += dt;
      if (Time<= Tau) Counts++;
      else 
	{
	  Time = 0;
	  Poiss->Fill(Counts);
	  Counts = 0;
	}
  
    }
  ca->cd(1);
  His->GetXaxis()->SetTitle("Time interval(#musec)");
  double TotRate = 1e3/His->GetMean();
  char TotRateString[100];
  sprintf(TotRateString,"Total Rate (for Crab 2-10 keV) = %.1f (kHz)", TotRate);
  TText *totratetext =new  TText(40,300,TotRateString);
  His->SetTitle(Title);
  His->Draw();
  totratetext->Draw();
  ca->cd(2);
  char titleNumber[100];
  sprintf(titleNumber,"Cluster's Number - Hold Time (#musec)= %.1f",Tau);
  Poiss->GetXaxis()->SetTitle(titleNumber);
  Poiss->SetTitle(Title);
  Poiss->Draw();
  //TF1 *Poisson1 = new TF1("Poisson1","[0]*TMath::PoissonI(x,[1])",0,20);
  //  Poiss->Fit("Poisson1");
}

void ClusterH(int MixID,double Thickness,double Pressure,double Tau)//legge tutti gli istogrammi del cluster number e fa la somma pesando con probabilita' di rivelazione(in funzione dell'energia)
{
  SetOption(MixID, Thickness, Pressure);
  TH1D *Hist = new TH1D("Hist","Hist",100,0,300);
  double Energy[18]={2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7.,7.5,8,8.5,9,9.5,10,10.5};
  double TotalRate = BinRate(MixID,Thickness,Pressure,2.,10.);
  double WeightV[18];
  for (int i=0; i<17;i++)
    {
      double Weight = BinRate(MixID,Thickness,Pressure,Energy[i],Energy[i+1])/TotalRate;
      Experiment->Generalsetup();
      TH1D *h = Experiment->GetClusterH(Energy[i]);//Pixel in Cluster
      Hist->Add(h,Weight);
      h->Draw();
      gPad->Update();
     
      WeightV[i] = Weight;
    }
  
  TCanvas *aa = new TCanvas("Cluster","Cluster",1200,800);
  Hist->GetXaxis()->SetTitle("Number of Pixel");
  TString Name = Experiment->GetName();
  Hist->SetTitle(Name);

  TText *t = new TText(150,1000,"Source = Crab (2 - 10 Kev)");
 
  Hist->Draw();
  t->Draw();
  NichDistribution(TotalRate*1e-6,Tau,Name);
  double Sum=0;
  for (int i=0;i<12;i++)
    {
      
      Sum +=WeightV[i];
      std::cout<<" "<<Sum<<std::endl;
    }
 
  std::cout<<Sum<<std::endl;
}
//Dopo SetOption e SetEnergy,  Apre il Tree dei dati MC e dell'analisi e visualizza funzione di risposta angolare
void ViewMC_Analys()
{
  TString McTreeName = DataName;
  McTreeName +=".root";
  TTree *t = (TTree*)filetree->Get("EventsTreef");
  t->AddFriend("EventsTree",McTreeName);
  std::cout<<McTreeName<< "  " <<AnalysisName<<std::endl;
  TCanvas *TreeMC_An = new TCanvas("Emission Direction Reconstruction","Emission Direction Reconstruction",1200,800);
  TreeMC_An->Divide(2,2);
  TreeMC_An->cd(1);
  t->StartViewer();
  t->Draw("AnClusterdim>>h");
  TH1D *h1 = (TH1D*)gDirectory->Get("h");
  double CluMean = h1->GetMean();
  double CluRMS= h1->GetRMS();
  TString Select = "AnClusterdim< ";
  double CluMin = CluMean - CluRMS;
  double CluMax = CluMean + CluRMS;
  Select += CluMax;
  Select+=" && AnClusterdim >";
  Select+= CluMin;
  TreeMC_An->cd(2);
  t->Draw("sqrt((RealImpactX-BaricenterX)^2 + (RealImpactY-BaricenterY)^2)/sqrt(MaxS)",Select);
  //delete h1;
  TreeMC_An->cd(3);
  t->Draw("sqrt((RealImpactX-BaricenterX)^2 + (RealImpactY-BaricenterY)^2)/sqrt(MaxS)","");
  /*t->Draw("InitialPhi:Theta2");
  TreeMC_An->cd(2);
  t->Draw("InitialPhi:ThetaI");
  TreeMC_An->cd(3);
  t->Draw("InitialPhi-Theta2");
  TreeMC_An->cd(4);
  t->Draw("InitialPhi-ThetaI");
  TCanvas *TreeMC_An2 = new TCanvas("Emission Direction Reconstruction2","Emission Direction Reconstruction2",1200,800);
  TreeMC_An2->Divide(2,2);
  TreeMC_An2->cd(1);
  t->Draw("abs(InitialPhi-Theta2)","");
  TreeMC_An2->cd(2);
  t->Draw("abs(InitialPhi-ThetaI)","");
  TreeMC_An2->cd(3);
  t->Draw("Theta2:ThetaI");*/
  /*  TreeMC_An->cd(1);
  t->Draw("InitialPhi-ThetaI>>h1");
  TH1D *h1 = (TH1D*)gDirectory->Get("h1");
  h1->GetXaxis()->SetTitle("Angle(rad)");
  h1->SetTitle("Reconstruction First Step");
  TreeMC_An->cd(2);
  t->Draw("InitialPhi-Theta2>>h2");
  TH1D *h2 = (TH1D*)gDirectory->Get("h2"); 
  h2->GetXaxis()->SetTitle("Angle(rad)");
  h2->SetTitle("Reconstruction Second Step");
  TreeMC_An->cd(3);
  t->Draw("abs(InitialPhi-ThetaI)>>h3","");
  TH1D *h3 = (TH1D*)gDirectory->Get("h3");
  h3->GetXaxis()->SetTitle("Angle(rad)");
  h3->SetTitle("Reconstruction First Step");
  TreeMC_An->cd(4);
  t->Draw("abs(InitialPhi-Theta2)>>h4");
  TH1D *h4 = (TH1D*)gDirectory->Get("h4");
  h4->GetXaxis()->SetTitle("Angle(rad)");
  h4->SetTitle("Reconstruction Second Step");
  std::cout<<Name<<std::endl;
  */}
//TH1F *f=(TH1F*)gPad->FindObject("h2")
void Convolution(int run,TH1F *DetMC,double Degree, double Theta)
{
  // double kPI = TMath::Pi(); 
  double PolarizationDegree = 1.;
  double PolarizationAngle = 30.;
  TF1  *PhiDistribution=new TF1("PhiDistribution","[0] + [1]*pow((cos(x-[2])), 2.0)", 0.0, 2.0*kPI);
  PhiDistribution->SetParameter(0, (1-Degree));
  PhiDistribution->SetParameter(1, 2*Degree); // ?
  PhiDistribution->SetParameter(2, Theta);
  TF1 *DetectorF = new TF1("DetectorF ","TMath::Gaus(x,[0],[1])",-kPI,kPI);
  TF1 *DetectorF2 = new TF1("DetectorF ","-x^2+5",-kPI,kPI);
  DetectorF->SetParameters(0,0.25);
  
  TH1D *PhiRecH = new TH1D("Real","Real",100,0.0,2.0*kPI);

  TH1D *PhiRecH4 = new TH1D("His4MC","His4MC",100,0.0,2.0*kPI);
  for (int i=0 ; i < run ; i++)
   { 
     double Phi = PhiDistribution->GetRandom();
     
     double PhiRec4 = Phi + DetMC->GetRandom();
     
     if (PhiRec4>2.*kPI)
       PhiRec4-=2*kPI;
     if (PhiRec4<0)
       PhiRec4+=2*kPI;
     PhiRecH->Fill(Phi);

     PhiRecH4->Fill(PhiRec4);
     if (i % 1000000==0)std::cout<<" i ="<<i<<std::endl;
   } 
  TCanvas *a = new TCanvas("D","D",1200,800);
  a->Divide(2,2);
  a->cd(1);
  PhiRecH->Draw("ap");
  TF1 *fufi1 = new TF1("fufi","[0]+[1]*(cos(x-[2]))^2",0,2*kPI);
  fufi1->SetParLimits(2,-kPI,kPI);
  fufi1->SetParLimits(0,0.,50000);
  fufi1->SetParLimits(1,0.,50000);
  PhiRecH->Fit("fufi","E","",0.,2*kPI);  
  a->cd(2);

  a->cd(3);

  a->cd(4);
  PhiRecH4->Draw("ap");
  PhiRecH4->Fit("fufi","E","",0.,2*kPI);  
  double A    = fufi1->GetParameter(0);
  double B    = fufi1->GetParameter(1);
  double  mu   = 100.0 * B / (2.*A+B);
  Float_t AError2 = fufi1->GetParError(0);
  Float_t BError2 = fufi1->GetParError(1);
  Float_t AngleError2 = fufi1->GetParError(2)*180.0/TMath::Pi();
  Float_t ModulationFactorError2 = 100*2.0*(A*BError2+B*AError2)/((B+2*A)*(B+2*A));
  std::cout<<" Mu = "<<mu<<"+-"<<ModulationFactorError2 <<" Angle Error="<<AngleError2<<std::endl;
} 
    
TH1F *GetHistogram()
{
 TH1F *H = (TH1F*)gPad->FindObject("htemp"); 
 return H;
}
void OpenBig(TString FileBig = "BigFile.root")
{
  std::cout<<FileBig<<std::endl;
  file = new TFile(FileBig);
 
  tree = (TTree*) file->Get("BigTree");
  tree->SetDirectory(0);
  tree->StartViewer();

}
void BigDraw(TString Select, TString FileBig = "BigFile.root")
{
  OpenBig(FileBig);
  TCanvas *ca = new TCanvas("ca","ca",1200,800);
  ca->Divide(2,2);
  ca->cd(1);
  tree->Draw("sqrt(EffCut)*MuCut",Select);
  TH1* H = (TH1*)gPad->FindObject("htemp");
  double MaxMuep=0;
  int Nbins=H->GetNbinsX();
  int i=Nbins;
  double Control=0;
  while (Control==0)
    {
      Control = H->GetBinContent(i);
      MaxMuep=H->GetBinCenter(i);
      i-=1;
    }
  std::cout<<"Max MuEf  = "<<MaxMuep<<"   "<<std::endl;
  Select+="&& sqrt(EffCut)*MuCut >";
  double MuepCut = MaxMuep - 10;
  Select+= MuepCut;
  ca->cd(2);
  tree->Draw("sqrt(Eff)*Mu",Select);
  tree->Scan("sqrt(EffCut)*MuCut:MuCut:Energy:MixId:Thickness:Pressure",Select);
}
void BigMaxMu(TString Select, TString FileBig = "BigFile.root")
{
  delete tree;
  OpenBig(FileBig);
  TCanvas *ca = new TCanvas("ca","ca",1200,800);
  ca->Divide(2,2);
  ca->cd(1);
  tree->Draw("Mu",Select);
  TH1* H = (TH1*)gPad->FindObject("htemp");
  double MaxMu=0;
  int Nbins=H->GetNbinsX();
  int i=Nbins;
  double Control=0;
  while (Control==0)
    {
      Control = H->GetBinContent(i);
      MaxMu=H->GetBinCenter(i);
      i-=1;
    }
  std::cout<<"Max Mu  = "<<MaxMu<<"   "<<std::endl;
  Select+="&& Mu>";
  double MuepCut = MaxMu - 10;
  Select+= MuepCut;
  ca->cd(2);
  tree->Draw("Mu",Select);
  tree->Scan("Mu:Mufirst:Energy:MixId:Thickness:Pressure",Select);
}
/*void OpenPixmap5_8()
{
  TCanvas *PixmapC = new TCanvas("PixmapC","PixmapC",1200,800);
  PixmapC->Divide(2,2);
  file = new TFile("BigFile.root");
  tree = (TTree*) file->Get("BigTree");
  tree->Draw("Mu:Energy","Pressure==1 && Thickness==0.5 && MixId==2");
  file->Close();
  file = new TFile("BigFile8.root");
  tree = (TTree*) file->Get("BigTree");
  TTree *Pix8 = tree->CloneTree();
   
}
*/
void ImagingComp(int MixID,double Thickness,double Pressure)
{
  double Energy[18] = {1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10};
  double Xrms[18];
  double Yrms[18];
  SetOption(MixID,Thickness,Pressure);
  for (int i = 0; i<18;i++)
    {
      SetEnergy(Energy[i]);
      treemix->Draw("ImpactX-RealImpactX>>hx");
      TH1D *hx = (TH1D*)gPad->FindObject("hx");
      treemix->Draw("ImpactY-RealImpactY>>hy");
      TH1D *hy = (TH1D*)gPad->FindObject("hy");
      Xrms[i] = hx->GetRMS(); 
      Yrms[i] = hy->GetRMS();
      gPad->Clear();
    }
  TGraph *ImagX = new TGraph(18,Energy,Xrms);
  TGraph *ImagY = new TGraph(18,Energy,Yrms);
  TCanvas *ca = new TCanvas("ca","ca",1200,800);
  ImagX->GetYaxis()->SetTitle("RMS X (#mum)");
  ImagX->GetXaxis()->SetTitle("Energy (keV)");
  ImagX->GetYaxis()->SetRangeUser(0,200);
  ImagX->SetMarkerStyle(2);
  ImagX->SetMarkerColor(2);
  ImagX->Draw("ap");
  TCanvas *ca2 = new TCanvas("ca2","ca2",1200,800);
  ImagY->GetYaxis()->SetTitle("RMS Y (#mum)");
  ImagY->GetXaxis()->SetTitle("Energy (keV)");
  ImagX->GetYaxis()->SetRangeUser(0,200);
  ImagY->SetMarkerStyle(2);
  ImagY->SetMarkerColor(2);
  ImagY->Draw("ap");
}

/*void BitRateO(int MixID,double Thickness,double Pressure,double Tau)
{
  SetOption(MixID, Thickness, Pressure);
  TH1D *Hist = new TH1D("Hist","Hist",100,0,300);
  TH1D *HPixel = new TH1D("Hist","Hist",100,0,300);
  TH1I *Poiss = new TH1I("","",20,0,20);
  double Energy[18]={2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7.,7.5,8,8.5,9,9.5,10,10.5};
  double TotalRate = BinRate(MixID,Thickness,Pressure,2.,10.);
  double TotalTime = -1;
  double PixelsinTimeInt = 0;
  for (int i=0; i<17;i++)
    {
      double Weight = BinRate(MixID,Thickness,Pressure,Energy[i],Energy[i+1])/TotalRate;
      Experiment->Generalsetup();
      TH1D *h = Experiment->GetClusterH(Energy[i]);
      Hist->Add(h,Weight);
      h->Draw();
      gPad->Update();
    }
    double Time =0;
    int Counts =0;
    double DeadTime = 0;
    for (int i=0;i<10000;i++)
    {
      double dt = 0;
      while (Time <= DeadTime)
	{
	  double csi = rnd->Uniform();
	  dt = -1.0e6/ TotalRate * log(1. - csi);//tempo tra un fotone e il successivo
	  TotalTime += dt;//tempo totale misura
	  Time += dt;
	}
      Time = 0;
      while (Time  < Tau)
	{
	  double csi = rnd->Uniform();
	  dt = -1.0e6/ TotalRate * log(1. - csi);
	  TotalTime += dt;
	  Time += dt;
	  Counts++;
	}
      // std::cout<<Counts<<std::endl;
      double TotalPixel = 0;
      for(int NT =0; NT <Counts;NT++)
	{
	  double NPixel =  Hist->GetRandom(); 
	  TotalPixel +=  NPixel;
	}
      PixelsinTimeInt += TotalPixel;
      HPixel->Fill(TotalPixel);
      DeadTime = TotalPixel / 10.e6;
      Time = 0;
      Poiss->Fill(Counts);
      Counts = 0;
    }
  TCanvas *ca= new TCanvas("Cluster","Cluster",1200,800);
  ca->Divide(2,2);
  ca->cd(1);
  HPixel->Draw();
  ca->cd(2);
  Hist->Draw();
  ca->cd(3);
  Poiss->Draw();
  std::cout<<" Total Time = "<<TotalTime<<" N tot Pixel ( = N byte)"<< PixelsinTimeInt<<
    " tot byte / Tot Time( microS ) = "<<PixelsinTimeInt/TotalTime<<  std::endl;
}*/


void BitRate(int MixID,double Thickness,double Pressure,double Tau)
{
  SetOption(MixID, Thickness, Pressure);
  // TH1D *Hist = new TH1D("Hist","Hist",100,0,500);
  //TH1D *HPixel = new TH1D("PHist","PHist",200,0,1000);
  TH1D *Hist = new TH1D("Window Dimension","Window Dimension",1000,0,10000);
  TH1D *Histc = new TH1D("Cluster Dimension","Cluster Dimension",100,0,1000);
  TH1D *HPixel = new TH1D("Eff. Window","Eff. Window",2000,0,20000);
  TH1I *Poiss = new TH1I("","",20,0,20);
  
 
  
  double Energy[18]={2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7.,7.5,8,8.5,9,9.5,10,10.5};
  double TotalRate = BinRate(MixID,Thickness,Pressure,2.,10.);
  double TotalTime = -1;
  double PixelsinTimeInt = 0;
  double NPhotonLoss =0;
  double NPhotonDet =0;
  for (int i=0; i<17;i++)
    {
      double Weight = BinRate(MixID,Thickness,Pressure,Energy[i],Energy[i+1])/TotalRate;
      Experiment->Generalsetup();
     
      TH1D *h = Experiment->GetRecordDimH(Energy[i]);
      std::cout<<" Pixel mean  "<<h->GetMean()<<" at Energy =  "<<Energy[i]<<std::endl;
      Hist->Add(h,Weight);
      
      h->Draw();
      
      gPad->Update();
    
    }
  for (int i=0; i<17;i++)
  { 
    double Weight = BinRate(MixID,Thickness,Pressure,Energy[i],Energy[i+1])/TotalRate;
    Experiment->Generalsetup();
    TH1D *hc = Experiment->GetClusterH(Energy[i]);
    Histc->Add(hc,Weight);
    hc->Draw();
    gPad->Update();
  }
int Counts =1;
  double DeathTime = 0;
  int Mode=0;
  double TotDeath=0; ;
  double T0 =0;
  for (int i=0 ; i<100000; i++)
    
    {
      double dt = 0;
     
      double csi = rnd->Uniform();
      dt = -1.0e6/ TotalRate * log(1. - csi);
      TotalTime+=dt;
      double T1 = T0+Tau;
      //
      if (i<100)
	{
	  std::cout<<"Hold end ="<<T1<<"     Time = "<<TotalTime<<std::endl;
	}
      //
      if (TotalTime<T1)
	{
	  Counts++;
	}
      else 
	{
	  
	  double TotalPixel = 0;
	  for(int NT =0; NT <Counts;NT++)
	    {
	      double NPixel =  Hist->GetRandom(); 
	      TotalPixel +=  NPixel;
	    }
	  PixelsinTimeInt += TotalPixel;
	  HPixel->Fill(TotalPixel);
	  DeathTime = TotalPixel / 10.;
	  TotDeath+=DeathTime;
	  Poiss->Fill(Counts);
	  NPhotonDet += Counts;
	  while (TotalTime<T1+DeathTime)
	    {
	      // if (i<1000)std::cout<<TotalTime<<" = Tot Death ="<<T1+DeathTime<<std::endl;
	      csi = rnd->Uniform();
	      dt = -1.0e6/ TotalRate * log(1. - csi);
	      TotalTime += dt;
	      NPhotonLoss++;
	    }
	  //
	  if (i<100)
	    {
	      std::cout<<"  Number of track in  Hold="<<Counts<<"     Death time="<<DeathTime<<std::endl;
	      std::cout<<" new Event at Time= "<<TotalTime <<std::endl;
	    }
	  //
	  T0 = TotalTime;
	  Counts = 1;
	}
       
    }
  TCanvas *ca= new TCanvas("Cluster","Cluster",1200,800);
  ca->Divide(2,2);
  ca->cd(1);
  //HPixel->Draw();
  Histc->Draw();
  ca->cd(2);
  Hist->Draw();
  ca->cd(3);
  Poiss->Draw();
  TCanvas *ct = new TCanvas ("ct","ct");
  TMDP aa;
  SetOption(MixID,Thickness,Pressure);
  TGraph *MuvsEn2 = GetMuvsEn();
  TGraph *EffvsEn2 = GetEffvsEn();
  aa.ReadDatafromGraph(MuvsEn2,EffvsEn2,2,1.0,1.0);
  TGraph *Eff =aa.GetEfficiencyGraph(); 
  ca->cd(4);
  //  Eff->Draw("ap");
  std::cout<<" ******************"<<std::endl;
  std::cout<<" Total Time ="<<TotalTime<<"   N tot Pixel (= N byte) ="<< PixelsinTimeInt<<"     Megabyte per Second ="<<1e6/(1024*1024)*PixelsinTimeInt/TotalTime<<"    Total  Time Death= "<<TotDeath<< std::endl;
  std::cout<<" Total Death Time / Time (%) ="<<100*TotDeath/TotalTime<<std::endl;
  std::cout<<"--------Number of Detected Photons = "<<NPhotonDet<<"----Number of UnDetected Photons ="<<NPhotonLoss <<" Detection (%)= "<<100 * NPhotonDet/(NPhotonDet+NPhotonLoss) <<std::endl;
  TText *t1=new TText(0.1,0.9,"Photon Rate= 00 kHz (Crab 2-10 keV) Mixture ---");
  TText *t2=new TText(0.1,0.8,"Trigger Level = ---- electrons, Window Width = --");
  TText *t3=new TText(0.1,0.7,"Clock 10 MHz, Death Time = 00.0%");
  TText *t4=new TText(0.1,0.6,"Effective Photon Rate = -- kHz");
  TText *t5=new TText(0.1,0.5,"Rate = -- Megabyte / sec");
  t1->Draw();
  t2->Draw();
  t3->Draw();
  t4->Draw();
  t5->Draw();
}
void ModulC()
{
  TH1D *ModulationCutH = (TH1D*) filetree->Get("ModulationCut");
  TF1 *fuficut = (TF1*) ModulationCutH->GetListOfFunctions()->FindObject("fufi");
  double A2    = fuficut->GetParameter(0);
  double B2    = fuficut->GetParameter(1);
  double phi02 = fuficut->GetParameter(2);
  Float_t AError = fuficut->GetParError(0);
  Float_t BError = fuficut->GetParError(1);
  Float_t ModulationFactorError = 2.0*(A2*BError+B2*AError)/((B2+2*A2)*(B2+2*A2));
  ModulationCutH->SetTitle("Modulation Factor");
  ModulationCutH->GetXaxis()->SetTitle("#phi (rad)");
  ModulationCutH->Draw();
  double mucut   = 100.0 * B2 / (2.*A2+B2); 
  
  char MufitC[100];
  sprintf(MufitC,"#mu = %.1f   #Delta #mu = %.2f ",mucut,ModulationFactorError);
  TLatex *t = new TLatex(3,20,MufitC);
  t->Draw();
}
TH1F *AngularResp()
{
  TString McTreeName = DataName;
  McTreeName +=".root";
  TTree *t = (TTree*)filetree->Get("EventsTreef");
  t->AddFriend("EventsTree",McTreeName);
  std::cout<<McTreeName<< "  " <<AnalysisName<<std::endl;
  t->Draw("InitialPhi-Theta2>>h2");
  TH1F *h2 = (TH1F*)gDirectory->Get("h2"); 
  h2->GetXaxis()->SetTitle("Angle (rad)");
  h2->SetTitle("Reconstruction Second Step");
  return h2;
}
//confronto fattori di Modulazione tra Pixmap 50 e 80
void PixmapEval()
{
  TCanvas *PEC = new TCanvas("Pix","Pix",1200,800);
  PEC->Divide(3,2);
  TCanvas *temp=new TCanvas("t","t",800,400);
  for(int Mixid = 2; Mixid<8; Mixid++)
    {
      SetOption(Mixid,1,1);
      temp->cd();
      TGraph *G5=GetMuvsEn();
      PEC->cd(Mixid-1);
      G5->Draw("alp");
      SetOption(Mixid,1,1,80);
      temp->cd();
      TGraph *G8=GetMuvsEn();
      G8->SetLineColor(2);
      PEC->cd(Mixid-1);
      G8->Draw("lp");
    }
}
void PixmapMix(int Mixid,double Th =1, double Pr=1)
{
  TCanvas *PEC = new TCanvas("Pix","Pix",1200,800);
  TCanvas *temp=new TCanvas("t","t",800,400);
  SetOption(Mixid,Th,Pr,50);
  temp->cd();
  TGraph *G5=GetMuvsEn();
  
  double *En = G5->GetX();
  double *Mu5 =  G5->GetY();
  PEC->cd();
  G5->Draw("alp");
  std::cout<<" Mu5 1 =  "<<Mu5[1]<<std::endl;
  SetOption(Mixid,Th,Pr,80);
  temp->cd();
  TGraph *G8=GetMuvsEn();
  double *Mu8 = G8->GetY();
  G8->SetLineColor(2);
  PEC->cd();
  G8->Draw("lp");
  TH1D *DeltaH = new TH1D("Delta","Delta",50,-10,10);
  TH1D *DeltaH2 = new TH1D("DeltaH2","DeltaH2",100,-100,100);
  //  for (int i=0;i<22;i++)
  for (int i=1;i<14;i++)//tolto il primo valore, a 1.5 keV, poco significativo egli ultimi, alle alte energie, poca differenza
    {   
      DeltaH->Fill(Mu5[i]-Mu8[i]);
      double MaxM = TMath::Max(Mu5[i],Mu8[i]);
      DeltaH2->Fill( 100.*(Mu5[i]-Mu8[i]) /(1.0 * MaxM));
      std::cout<<" %Gain =   "<<100.*(Mu5[i]- Mu8[i]) /(1.0 * MaxM)  <<std::endl;
    }
  TCanvas *c2 = new TCanvas("c2","delta",1200,800);
  c2->cd();
  DeltaH->Draw();
  TCanvas *c3 = new TCanvas("c3","delta%",1200,800);
  c3->cd();
  DeltaH2->Draw();
}
//TH1F *f=(TH1F*)gPad->FindObject("h2")
void Convolution2(int run,TH1F *DetMC,double Degree, double Theta)
{
  // double kPI = TMath::Pi(); 
  double PolarizationDegree = 1.;
  double PolarizationAngle = 30.;
  TF1  *PhiDistribution=new TF1("PhiDistribution","[0] + [1]*pow((cos(x-[2])), 2.0)", 0.0, 2.0*kPI);
  PhiDistribution->SetParameter(0, (1-Degree));
  PhiDistribution->SetParameter(1, 2*Degree); // ?
  PhiDistribution->SetParameter(2, Theta);
  
  TH1D *PhiRecH = new TH1D("Real","Real",100,0.0,2.0*kPI);
  TH1D *PhiRecH4 = new TH1D("His4MC","His4MC",100,0.0,2.0*kPI);
  for (int i=0 ; i < run ; i++)
   { 
     double Phi = PhiDistribution->GetRandom();
     double PhiRec4 = Phi + DetMC->GetRandom();
    
     if (PhiRec4>2.*kPI)
       PhiRec4-=2*kPI;
     if (PhiRec4<0)
       PhiRec4+=2*kPI;
     PhiRecH->Fill(Phi);
     
     PhiRecH4->Fill(PhiRec4);
     if (i % 1000000==0)
	{ 
	  std::cout<<" ********** Nevent = "<<i<<std::endl;
	}
    
   } 
  TCanvas *a = new TCanvas("D","D",1200,800);
  a->Divide(2,2);
  a->cd(1);
  PhiRecH->Draw("ap");
  TF1 *fufi1 = new TF1("fufi","[0]+[1]*(cos(x-[2]))^2",0,2*kPI);
  fufi1->SetParLimits(2,-kPI,kPI);
  fufi1->SetParLimits(0,0.,run*1.0);
  fufi1->SetParLimits(1,0.,run*1.0);
  PhiRecH->Fit("fufi","E","",0.,2*kPI);  
  a->cd(2);
  a->cd(3);
  a->cd(4);
  PhiRecH4->Draw("ap");
  PhiRecH4->Fit("fufi","E","",0.,2*kPI);  
  double A    = fufi1->GetParameter(0);
  double B    = fufi1->GetParameter(1);
  double  mu   = 100.0 * B / (2.*A+B);
  Float_t AError2 = fufi1->GetParError(0);
  Float_t BError2 = fufi1->GetParError(1);
  Float_t AngleError2 = fufi1->GetParError(2)*180.0/TMath::Pi();
  Float_t ModulationFactorError2 = 100*2.0*(A*BError2+B*AError2)/((B+2*A)*(B+2*A));
  std::cout<<" Mu = "<<mu<<"+-"<<ModulationFactorError2 <<" Angle Error (gradi)="<<AngleError2<<std::endl;
} 

// Modello della variazione dell'angolo di polarizzazione di Cygnus X-1 in funzione dell'energia
TGraph *CignModel(Bool_t ReturnRad=true)//Radianti
{ 
  ifstream in;
  in.open("Tre.txt");
  Float_t E[2000], Th[2000];
  Int_t nlines=0;
   while(1) {
    
    in >> E[nlines] >> Th[nlines];
    if(!in.good())break;
    if (ReturnRad==true)
      Th[nlines] = TMath::Pi()*Th[nlines]/ 180. ;
    nlines++;
   }
   TGraph *AngModel=new TGraph(nlines,E,Th);
   AngModel->Draw("alp");
   return AngModel;

}
TH1F *AngleVar(TGraph *ThetaG, double E,int run)
{
  TF1  *PhiDistribution=new TF1("PhiDistribution","[0] + [1]*pow((cos(x-[2])), 2.0)", 0.0, 2.0*kPI);
  PhiDistribution->SetParameter(0, (1-0.05));
  PhiDistribution->SetParameter(1, 2*0.05); // ?
 
 
  TH1F *PhiIN = new TH1F("Real","Real",100,0.0,2.0*kPI);
  for (int i=0 ; i < run ; i++)
    {
      double En=rnd->Gaus(E,0.2*En);
      double Theta= ThetaG->Eval(En);
      PhiDistribution->SetParameter(2, Theta);
      double Phi = PhiDistribution->GetRandom();
      PhiIN->Fill(Phi);
    }
  PhiIN->Draw();
  /*TF1 *fufi1 = new TF1("fufi","[0]+[1]*(cos(x-[2]))^2",0,2*kPI);
  fufi1->SetParLimits(2,-kPI,kPI);
  fufi1->SetParLimits(0,0.,run*1.0);
  fufi1->SetParLimits(1,0.,run*1.0);
  PhiIN->Fit("fufi","E","",0.,2*kPI);*/
  return PhiIN;
}
//Fit MC osservando Cygnus X-1 con la miscela selezionata in SetOption() tenendo conto della variabilita' dell'angolo di polarizzazione in accordo con il modello del grafico CignModel()  
void Convolsp(int run,double En)
{
  double polarizationDegree=0.1;
  SetOption(5,1,1,50);
  SetEnergy(En);
  TH1F *DetMC=AngularResp();
  TGraph *ACign=CignModel();
  double SigmaThp = ACign->Eval(En+0.2*En);
  double SigmaThm = ACign->Eval(En-0.2*En);
  double ThetaM=ACign->Eval(En);
  double Sigma=(abs(SigmaThp-ThetaM) + abs(SigmaThm -ThetaM) )/2.;
  //TH1F *PhiD=AngleVar(ACign,En,20000);
  TF1  *PhiDistribution=new TF1("PhiDistribution","[0] + [1]*pow((cos(x-[2])), 2.0)", 0.0, 2.0*kPI);

  PhiDistribution->SetParameter(0, (1-polarizationDegree));
  PhiDistribution->SetParameter(1, 2*polarizationDegree);
  PhiDistribution->SetParameter(2, ThetaM);
  std::cout<<"Sigma"<<Sigma<<"Theta"<<ThetaM<<std::endl;
  TH1D *PhiRecH = new TH1D("Real","Real",100,0.0,2.0*kPI);
  TH1D *PhiRecH4 = new TH1D("His4MC","His4MC",100,0.0,2.0*kPI);
  for (int i=0 ; i < run ; i++)
   { 
     double Phi = PhiDistribution->GetRandom();
     //Phi+=rnd->Gaus(0,Sigma);
     double PhiRec4 = Phi + DetMC->GetRandom();//+rnd->Gaus(0,Sigma);
     PhiRec4+=rnd->Gaus(0,Sigma);
     if (PhiRec4>2.*kPI)
       PhiRec4-=2*kPI;
     if (PhiRec4<0)
       PhiRec4+=2*kPI;
     PhiRecH->Fill(Phi);
     
     PhiRecH4->Fill(PhiRec4);
     if (i % 1000000==0)
	{ 
	  std::cout<<" ********** Nevent = "<<i<<std::endl;
	}
    
   } 
  TCanvas *a = new TCanvas("D","D",1200,800);
  a->Divide(2,2);
  a->cd(1);
  PhiRecH->Draw("ap");
  TF1 *fufi1 = new TF1("fufi","[0]+[1]*(cos(x-[2]))^2",0,2*kPI);
  fufi1->SetParLimits(2,-kPI,kPI);
  fufi1->SetParLimits(0,0.,run*1.0);
  fufi1->SetParLimits(1,0.,run*1.0);
  PhiRecH->Fit("fufi","E","",0.,2*kPI);  
  a->cd(2);
  DetMC->Draw();
  a->cd(3);
  ACign->Draw("alp");
  a->cd(4);
  PhiRecH4->Draw("ap");
  PhiRecH4->Fit("fufi","E","",0.,2*kPI);  
  double A    = fufi1->GetParameter(0);
  double B    = fufi1->GetParameter(1);
  double  mu   = 100.0 * B / (2.*A+B);
  Float_t AError2 = fufi1->GetParError(0);
  Float_t BError2 = fufi1->GetParError(1);
  Float_t AngleError2 = fufi1->GetParError(2)*180.0/TMath::Pi();
  Float_t ModulationFactorError2 = 100*2.0*(A*BError2+B*AError2)/((B+2*A)*(B+2*A));
  std::cout<<" Mu = "<<mu<<"+-"<<ModulationFactorError2 <<" Angle Error (gradi)="<<AngleError2<<"Angle (gradi) "<<fufi1->GetParameter(2)*180.0/TMath::Pi()<<std::endl;
} 


void ConnorsRec()
{
  Bool_t grad=false;
  TGraph *ACign=CignModel(grad);
  double En[5]={2.5,4,6,9};
  double Ener[5]={0.5,0.2*3.5,0.2*6,0.2*9};
  double Theta[5]={-1.129,-0.943,-0.79,-0.622};
  double Thetaerr[5]={0.002,0.0035,0.0081,0.027};
  
  for (int i=0;i<4;i++)
    {
      Theta[i]*=180./ TMath::Pi();
      Thetaerr[i]*=180./ TMath::Pi();
      std::cout<<" "<<Theta[i]<<"  "<<Thetaerr[i]<<std::endl;
    }
  TGraphErrors *Ar = new TGraphErrors(5,En,Theta,Ener,Thetaerr);
  TCanvas *can=new TCanvas("can","can",1200,800);
 
  //TCanvas *b=new TCanvas("b","b",1200,600);
  ACign->GetXaxis()->SetTitle("Energy (keV)");
  ACign->GetYaxis()->SetTitle("Polarization Angle");
  ACign->Draw("al");
  ACign->GetYaxis()->LabelsOption("<");
  ACign->Draw("al");
  Ar->Draw("p");
}
void ConnorsRecjetx()
{
  Bool_t grad=false;
  TGraph *ACign=CignModel(grad);
  double En[5]={2.5,4,6,9};
  double Ener[5]={0.5,0.2*3.5,0.2*6,0.2*9};
  double Theta[5]={-64.675,-54.69,-44.79,-34.69};
  double Thetaerr[5]={0.171,0.29,0.61,1.99};
  TGraphErrors *Ar = new TGraphErrors(5,En,Theta,Ener,Thetaerr);
  TCanvas *can=new TCanvas("can","can",1200,800);
 
  //TCanvas *b=new TCanvas("b","b",1200,600);
  
  ACign->Draw("al");
  ACign->GetYaxis()->LabelsOption("<");
  ACign->GetXaxis()->SetTitle("Energy (keV)");
  ACign->GetYaxis()->SetTitle("Polarization Angle");
  ACign->Draw("al");
  Ar->Draw("p");
}
