#include "TTreeAnalysis.h"

static const double PI = TMath::Pi();


TTreeAnalysis:: TTreeAnalysis (TString FileName,TString DataTreeName)
{
  keyName = FileName;
  TString FName = FileName;  
  FName += ".root";
  TString DimName = FileName;
  DimName += ".txt";
  TFile *f= new TFile(FName);
  t = (TTree*)f->Get(DataTreeName);
  //  t->SetDirectory(0);
  //  f->Close();
  Nentries=(int)t->GetEntries();
  rnd = new TRandom();
  Dimension = new TDimension();
  Dimension->ReadDimensionFromFile(DimName);//+.txt
  Pitch = Dimension->GetPitch();
  std::cout<<"Pitch : "<<Pitch<<std::endl;
  std::pair<int,int>  Pi=Dimension->GetPixel();
  int IMax = Pi.second;
  int JMax = Pi.first;
  XShift = ((JMax+3.)*Pitch*sqrt(3.)/4.)-Pitch/sqrt(3.)+7.;
  //XShift = ((JMax+3.)*Pitch*sqrt(3.)/4.)-Pitch/sqrt(3.)+7.-3507;;
  YShift = Pitch*(IMax+2.)/2.+25.;
  //YShift = Pitch*(IMax+2.)/2.+25.-3525;
  
  ElADC =  Dimension->GetElectronsADC();
  Avalance = Dimension->GetGain();
  std::cout<<" ElADC: "<<ElADC<<std::endl;
  std::cout<<" Avalance: "<<Avalance<<std::endl;
  MixID = 0;/////////////////////////////////
  
  t->SetBranchAddress("InitialPhi",&Phi);
  t->SetBranchAddress("MixID",&MixID);
  t->SetBranchAddress("Clusterdim",&Clusterdim);
  t->SetBranchAddress("Channel",Channel);
  t->SetBranchAddress("Charge",Charge);
  t->SetBranchAddress("DigiCharge",DigiCharge);
  
  t->SetBranchAddress("XInitial",&XIn);
  t->SetBranchAddress("YInitial",&YIn);
 
  std::cout<<" Number of Entries "<<t->GetEntries()<<std::endl;
    
  //  t->Scan();  
  t->GetEntry(0);

  std::cout<<MixID<<" "<<ElADC<<" "<<Avalance<<std::endl;
  Mixture = new TGasMixture(MixID, rnd,Dimension);//Settaggio Mixture
  TDimension *DimensionNorm = new TDimension();
  DimensionNorm->ReadDimensionFromFile("DimensionNorm.txt");
  MixtureRif = new TGasMixture(2, rnd,DimensionNorm);
 
  myDetector = new TDetector(Mixture,rnd,0,Dimension);
  Readout = TReadout(Dimension);
  WIonization = Mixture->GetWIonization();
  Pressure = Dimension->GetPressure();
  //
  //delete f;
  keyName += "Analysis.root";
  UpdateFile = new TFile(keyName,"RECREATE");
  // F.SetOption("RECREATE"); 
  // F.SetName(keyName);
  // F = new TFile(keyName,"RECREATE");
  //
} 

std::vector<ADC> TTreeAnalysis::ReadCluster(int NCluster)
{
  ClusterPulse = 0;
  Clr = 0;
  t->GetEntry(NCluster);
  std::vector<ADC> Digi;
  for(int ch=0;ch<Clusterdim;ch++)
    {
      ADC Pix;  // ADC: structure (Channel,Charge,DigiCharge) defined in TDetector.h 
      Pix.Channel=Channel[ch];
      Pix.Charge=Charge[ch];
      if (DigiCharge[ch]<1394) // Con piedistallo a 1400 ADC DigiCharge<1394 vuol dire soglia 6 ADC counts ovvero 3 sigma noise (= 2 ADC counts)
  	{
	  Pix.DigiCharge=DigiCharge[ch];
	
	  Clr += 1;
	}
      else Pix.DigiCharge=1400;
      Digi.push_back(Pix);
      //ClusterPulse += 1400 - DigiCharge[ch];// Charge[ch]/ElADC;   
      ClusterPulse += 1400 - Pix.DigiCharge;// Charge[ch]/ElADC;   
    }
  if(Clr>2)
    return Digi;
  else{
    std::cout<<" TOO SMALL "<<std::endl;
    Digi.clear();
    return Digi;
  }
  
}


void TTreeAnalysis::ClusterAnalysis(int NCluster)
{
  std::vector<ADC> Digi = ReadCluster(NCluster);
  std::cout << "========>>> Draw cluster with method myDrawSignal of TDetector.cpp" << std::endl;
  myDetector->myDrawSignal(Digi,Readout);  // disegna i pixels del cluster (v. TDetector.cpp) se il segnale (ped sub) e' > 0 (NO taglio a 3 sigma)
  int Clusterdim = Digi.size();
  SetClusterPar(Clusterdim);
  
  std::cout<<"===>> Params: Small radius /large radius / Weight: "<< a <<"-- / -- " << b << "-- / -- "<<c<<std::endl;
  TCluster Cluster(a,b,c);
  //SquareCheck
  std::vector<double> limits = Cluster.ClusterRecord (Digi,Readout,2000,10);
  //std::cout<<" Total Record  = "<<limits[4]<<" Raw Clustr Dim= "<<Clusterdim<<std::endl; // limits[4] dovrebbe essere la window size. (NdG 13/12/06)
  TBox *Box = new TBox(limits[0],limits[1],limits[2],limits[3]);
  Box->SetFillStyle(3004);
  Box->SetFillColor(5);
  Box->Draw();
  //
  std::vector<std::pair<double,double> > TriggerPos = Cluster.GetTriggeredPositionXY();
  //std::vector<std::pair<double,double> >:: iterator iter;
  const int DimTr = TriggerPos.size();
  vector <double> XTr;
  vector <double> YTr;
  for (int i= 0;i<DimTr;i++)
    {
      //std::cout<<" X = "<<TriggerPos[i].first<<" Y ="<< TriggerPos[i].second<<std::endl;   // Pixel che triggerano????  (NdG 13/12/06)

      XTr[i] = TriggerPos[i].first + Pitch/(sqrt(3.));
      YTr[i] = TriggerPos[i].second + Pitch;
     }
  TGraph *TrigPosGr = new TGraph(DimTr,&XTr[0],&YTr[0]);
  TrigPosGr->SetMarkerStyle(3);
  TrigPosGr->SetMarkerSize(2.0);
  TrigPosGr->SetMarkerColor(4);
  std::vector<Pixel> D = Cluster.ADCtoPixel(Digi,Readout);
  Cluster.Analysis(D);  //// <<<< ============================================
  double ImpactRecX = Cluster.GetImpactX();
  double ImpactRecY =Cluster.GetImpactY(); 
  double X= Cluster.GetBaricenterX();
  double Y= Cluster.GetBaricenterY();
  double Theta2 = Cluster.GetTheta2();
  double Theta = Cluster.GetTheta();
  double TM = Cluster.GetThirdMomentum();
  double MaxS =Cluster.GetMaxMom();
  double MinS = Cluster.GetMinM();
  double RealImpactX =XIn*10000+XShift;
  double RealImpactY =YIn*10000+YShift;
 
  std::pair<double,double>  realp = Cluster.Froc(RealImpactX,RealImpactY,Phi);
  std::pair<double,double>  analysisp = Cluster.Froc( ImpactRecX, ImpactRecY,Theta2);   //( ImpactRecX, ImpactRecX,Theta2);

  TArrow  *real=new  TArrow(RealImpactX,RealImpactY,realp.first,realp.second);
  TArrow  *analys=new  TArrow(ImpactRecX,ImpactRecY,analysisp.first,analysisp.second);
  TMarker *ImpactRec = new TMarker(ImpactRecX,ImpactRecY,28);
  TMarker *Baric = new TMarker(X,Y,28); 
  TMarker *RealI = new TMarker(RealImpactX,RealImpactY,28);
  ImpactRec->SetMarkerColor(2); //ROSSO IMPATTO RICOSTRUITO
  ImpactRec->SetMarkerSize(2);
  ImpactRec->Draw();
  Baric->SetMarkerColor(4);     //BLU BARICENTRO
  Baric->SetMarkerSize(2);
  Baric->Draw();
  RealI->SetMarkerColor(3);     //VERDE  IMPATTO VERO
  RealI->SetMarkerSize(2);
  RealI->Draw();
  real->SetLineColor(3);
  real->Draw();
  analys->SetLineColor(2);
  analys->Draw();
  //
  
  std::pair<double,double>  MaxAxdx = Cluster.Froc(X,Y,Theta);
  std::pair<double,double>  MinAxdx = Cluster.Froc(X,Y,Theta+PI/2.);
  std::pair<double,double>  MaxAxsx = Cluster.Froc(X,Y,Theta+PI);
  std::pair<double,double>  MinAxsx = Cluster.Froc(X,Y,Theta+PI*3./2);
  TLine *MinAx =new TLine(MinAxsx.first,MinAxsx.second,MinAxdx.first,MinAxdx.second);
  TLine *MaxAx =new TLine(MaxAxsx.first, MaxAxsx.second, MaxAxdx.first, MaxAxdx.second);
  MinAx->Draw();
  MaxAx->Draw();
  TText *texMa = new TText(X+200,Y+200,"Major Axis");
  TText *texma = new TText(X+250,Y+250,"Minor Axis");
 
  texMa->SetTextFont(42);
  texMa->SetTextSize(0.03);
  texMa->Draw();
  texma->SetTextFont(42);
  texma->SetTextSize(0.03);
  // texma->Draw();
  std::cout<<"==>> Theta0(rad) = " << Theta << " -> (gradi) = " << Theta*360/PI <<std::endl;
  std::cout<<"==>> Shape =" << MaxS/MinS << " abs(Mom III) =" << abs(TM) << std::endl;
  //
  //Modifiche grafiche
  double LargRad = b * sqrt(MaxS);
  double SmallRad = a * sqrt(MaxS);
  TEllipse *SREll = new TEllipse(X,Y,SmallRad);
  TEllipse *LREll = new TEllipse(X,Y,LargRad);
  SREll->SetLineColor(4);
  LREll->SetLineColor(7);
  SREll->Draw();
  LREll->Draw();
  //std::cout<<"  sqrt II =  "<<sqrt(MaxS)<<std::endl;
  
  //Fine modifiche
  
  TH1D *XProfile = new TH1D("XProfile","XProfile",20,-500,800);
  TH1D *YProfile = new TH1D("YProfile","YProfile",20,-600,600);
  
  //double RealB=sqrt(pow(RealImpactX- X,2.)+pow(RealImpactY- Y,2.));
  //double Cos = (X-RealImpactX)/RealB;
  //double Sin = (Y-RealImpactY)/RealB;
  
  //std::vector<Pixel> VectROT = Cluster.RotateVector(D,Cos,Sin,RealImpactX,RealImpactY);
  double Cos = cos(Theta2);
  double Sin = sin(Theta2);
  std::vector<Pixel> VectROT = Cluster.RotateVector(D,Cos,Sin,RealImpactX,RealImpactY);
  Cluster.ClusterProfile(VectROT, XProfile , YProfile);
  Cluster.Analysis(VectROT);
  double Xbar = Cluster.GetBaricenterX();
  //double Ybar = Cluster.GetBaricenterY();
  double sp = a*sqrt(Cluster.GetMaxMom()) + Xbar;
  double sp2 = -(a*sqrt(Cluster.GetMaxMom())) + Xbar;
  double lp = b*sqrt(Cluster.GetMaxMom()) + Xbar;
  double lp2 = -(b*sqrt(Cluster.GetMaxMom())) + Xbar;
  double xrecimp = Cluster.GetImpactX();
  double yrecimp = Cluster.GetImpactY();
  
  std::cout << " Impact X (MC ref.) = "<< RealImpactX  << " Impact Y (MC ref.) = "<< RealImpactY<< " Theta2 = "<< Theta2 << "Phi = "<< Phi << std::endl;  
  std::cout << " ImpactRec X = "<< ImpactRecX  << " ImpactRec Y = "<< ImpactRecY  << std::endl;
  std::cout << " Baricenter X = " << X <<" Baricenter Y = "<< Y <<std::endl;  
  std::cout <<"Shift X : " << XShift <<"  Shift Y : " << YShift << std::endl; 
    //std::cout << " Froc fin X =  "<<analysisp.first << " Froc fin Y = "<<analysisp.second << std::endl;  
  
  std::cout <<"Inpact X (Pixie ref, mm): "<< XIn <<" Inpact Y (Pixie ref, mm): "<< YIn  << std::endl;  

  TMarker *px = new TMarker(xrecimp,0.,28);
  TMarker *py = new TMarker(yrecimp,0.,28);//
  TLine *smallx = new TLine(sp,0.,sp,20.);
  TLine *smallx2 = new TLine(sp2,0.,sp2,20.);
  TLine *largex = new TLine(lp,0.,lp,20.);
  TLine *largex2 = new TLine(lp2,0.,lp2,20.);
  TCanvas *Profile = new TCanvas("Profile","Profile",1200,800);
  Profile->Divide(2,0);
  Profile->cd(1);
  XProfile->Draw();
  px->SetMarkerColor(2);
  px->Draw();
  smallx->SetLineColor(3);//verde
  smallx2->SetLineColor(3);//verde
  smallx->Draw();
  smallx2->Draw();
  largex->SetLineColor(4);//blu
  largex2->SetLineColor(4);
  largex->Draw();
  largex2->Draw();
  Profile->cd(2);
  YProfile->Draw();
  py->SetMarkerColor(2);
  py->Draw();
 
}


//PER GLORIA: tu devi prendere i valori con Pitch==80 (gli altri sono per la Pixmap da 50),Clusterdim e' il numero di cluster in cui e'arrivato qualcosa, anche 1 solo elettrone, quindi vanno ritarati i valori contando che molti pixel qui contati, vengono tagliati con i piedistalli.E' un lavorino che mi tocchera' prossimamente
void            TTreeAnalysis::SetFixedClusterPar(double ain,double bin, double cin)
{
  a=ain;
  b=bin;
  c=cin;

}
void            TTreeAnalysis::SetClusterPar(int Clu)//int Clusterdim sostituisco Clusterdim con Clu:
{
  double Erec = 0.4 + ClusterPulse * (ElADC/Avalance) * WIonization/1000.;
  if (Erec<2.5) {a=0.4;b=1.;}
  else if (Erec<3){a=0.9; b=1.5;}
  else if(Erec<4){a=1.1; b=2.0;}
  else if(Erec<5){a=1.3; b=2.3;}
  else if(Erec<6){a=1.5; b=2.5;}
  else if(Erec<7){a=1.55; b=2.7;}
  else if(Erec<8){a=1.65; b=2.8;}
  else if(Erec<9){a=1.7; b=3.0;}
  else if(Erec<12){a=2; b=4;}
  else {a=2.5; b=4;}
  double Dim=0;
  if(Pitch==80) Dim = Clu*1.6;
  else Dim = Clu;
  Dim = Clu;
  
  //  int ParCor = 1;
  double MixtureMod = 0.15*Mixture->GetElectronRange(Erec)/MixtureRif->GetElectronRange(Erec);
  
  //a*=1.1;Dati confronto Pixmap con fattore 1.1
  //b*=1.1;
  //if(ParCor){
    a+=MixtureMod;
    b+=MixtureMod;
    //}
  
  if (a>2.5)a=2.5;
  if (Dim>60)
    c = 2*(Dim-30);
  else if (Dim<20)
  c= 5*(30-Dim);
  else c=50;// c = 10000*9*Mixture->GetElasticMeanFreePath(Erec, "MOTT");
  if(c>200) c=200;
  //c=50;
}



void TTreeAnalysis::TreeAnalysis( Bool_t ParSelectFix )
{ 
  NStatistic = 0;
  std::cout<<Nentries<<std::endl;
  treef = new TTree("EventsTreef","analysis");
  //  treef->SetDirectory(0);
  double BaricenterX;
  double BaricenterY;
  double ImpactX;
  double ImpactY;
  double Theta2;
  double ThetaI;
  double RealX =0;
  double RealY =0;
  double ThirdMom;
  double MaxS;
  double ReImBar;
  double MinMom;
  int Clusterr;
  double Pulse;
  int RecordClusterdim;
  int RawClusterdim;
  int Recorddim;
  int SqDim;
  double ChargeSq[1000];
  treef->Branch("AnClusterdim",&Clusterr,"Clusterr/I");
  treef->Branch("BaricenterY",&BaricenterY,"BaricenterY/D");
  treef->Branch("BaricenterX",&BaricenterX,"BaricenterX/D");
  treef->Branch("ImpactX",&ImpactX,"ImpactX/D");
  treef->Branch("ImpactY",&ImpactY,"ImpactY/D");
  treef->Branch("RealImpactX",&RealX,"RealImpactX/D");
  treef->Branch("RealImpactY",&RealY,"RealImpactY/D");
  treef->Branch("Theta2",&Theta2,"Theta2/D");
  treef->Branch("ThirdMom",&ThirdMom,"ThirdMom/D");
  treef->Branch("MaxS", &MaxS, "MaxS/D");
  treef->Branch("ReImBar", &ReImBar, "ReImBar/D");
  treef->Branch("MinMom",&MinMom,"MonMom/D");
  treef->Branch("ThetaI",&ThetaI,"ThetaI/D");
  treef->Branch("Pulse",&Pulse,"Pulse/D");

  treef->Branch("RecordClusterdim",&RecordClusterdim,"RecordClusterdim/I");
  treef->Branch("RawClusterdim",&RawClusterdim,"RawClusterdim/I");
  treef->Branch("Recorddim",&Recorddim,"Recorddim/I");
  
  treef->Branch("SqDim",&SqDim,"SqDim/I");
  treef->Branch("ChargeSquare",ChargeSq,"ChargeSq[SqDim]/D");
  
  std::vector<ADC> Digi;
  for (int i=0;i<Nentries;i++) 
    {
      Digi = ReadCluster(i);
      if(Digi.size()>0)
	{
	  Clusterr =  GetClr();
	  RealX = XIn*10000+XShift;
	  RealY = YIn*10000+YShift;
	  Pulse = ClusterPulse;
	  if ( ParSelectFix == false )
	    {
	      SetClusterPar(Clusterr);
	    }
	  
	  TCluster Cluster(a,b,c);
	  std::vector<Pixel> D = Cluster.ADCtoPixel (Digi,Readout);// trasforma vettore(chan,charge,digicharge) in (posX, posY,charge)
	 
	  std::vector<double> RangeV = Cluster.ClusterRecord (Digi,Readout,5500,10);
	  std::vector<Pixel> RecordD = Cluster.ADCtoPixelinRange(Digi,RangeV,Readout);
	  RecordClusterdim =  RecordD.size();
	  RawClusterdim = D.size();
	  Recorddim = (int) RangeV[4];
	  std::vector<double> SquareV = Cluster.GetSquareCharge();
	  SqDim = SquareV.size();
	 
	  for (int is=0; is<SqDim ; is++)
	  {
	    ChargeSq[is] = SquareV[is];
	  }
	  Cluster.Analysis(D);
	  BaricenterX = Cluster.GetBaricenterX();
	  BaricenterY = Cluster.GetBaricenterY();
	  ImpactX = Cluster.GetImpactX();
	  ImpactY = Cluster.GetImpactY();
	  Theta2 = Cluster.GetTheta2();
	  ThirdMom = Cluster.GetThirdMomentum();
	  MaxS = Cluster.GetMaxMom(); 	
	  MinMom = Cluster. GetMinM();
	  ThetaI = Cluster.GetTheta();
	  
	  NStatistic += Cluster.GetStatistic();
	  
	  //double TM = Cluster.GetThirdMomentum();
	  ReImBar = sqrt(pow(RealX- ImpactX,2.)+pow(RealY- ImpactY,2.));
	  //double RealB=sqrt(pow(RealX- BaricenterX,2.)+pow(RealY- BaricenterY,2.));
	  //double ImB=sqrt(pow(ImpactX- BaricenterX,2.)+pow(ImpactY- BaricenterY,2.));
	  
	  //double CosBarReal = (BaricenterX-RealX)/RealB;
	  //double SinBarReal = (BaricenterY-RealY)/RealB;
	  //double xImpactMod = (ImpactX-RealX) * CosBarReal + (ImpactY-RealY) * SinBarReal; 
	  //double yImpactMod = (ImpactY-RealY) * CosBarReal - (ImpactX-RealX) * SinBarReal; 
	  
	  treef->Fill();
    
	}
    }
  std::cout<<" NStatistic  = "<<NStatistic<<std::endl;
}  


void TTreeAnalysis::FillH(double ModulationDegree)
{  
  gDirectory->Delete("ThirdMomHi");
  gDirectory->Delete("SecondRatioHi");
  gDirectory->Delete("ClusterDim");
  gDirectory->Delete("RecordDimH");
  gDirectory->Delete("RawClusterDimH");
  gDirectory->Delete("RecordClusterDimH");
 
  gDirectory->Delete("ModulationSimpleH");
  gDirectory->Delete("ModulationCutH");
  gDirectory->Delete("ModulationFirstStep");
  
 

  treef->Draw("MaxS/MinMom>>SecondRatioHi","MaxS/MinMom<100. && MaxS/MinMom>0.");
  treef->Draw("abs(ThirdMom)>>ThirdMomHi","abs(ThirdMom)<1.0e+8");
  treef->Draw("AnClusterdim>>ClusterDim(100,0,500)");
  treef->Draw("Recorddim>>RecordDimH(1000,0,10000)");
  treef->Draw("RawClusterdim>>RawClusterDimH(100,0,500)");
  treef->Draw("RecordClusterdim>>RecordClusterDimH(500,0,500)");
  
  RawClusterDimH = (TH1D*)gDirectory->Get("RawClusterDimH");
  RecordDimH = (TH1D*)gDirectory->Get("RecordDimH");
  RecordClusterDimH = (TH1D*)gDirectory->Get("RecordClusterDimH");
  ClusterDim = (TH1D*)gDirectory->Get("ClusterDim");    
  ThirdMomHi = (TH1D*)gDirectory->Get("ThirdMomHi");
  SecondRatioHi = (TH1D*)gDirectory->Get("SecondRatioHi");

  
  double ThirdMomMean = ThirdMomHi->GetMean();
  double SecondRatioMean = SecondRatioHi->GetMean();

  std::cout<<ThirdMomMean<<" "<<SecondRatioMean <<std::endl;
  
  char myCut[100];

  sprintf(myCut,"(MaxS/MinMom > %f) || (abs(ThirdMom) > %f)",SecondRatioMean/2.,ThirdMomMean/2.);
  std::cout<<myCut<<std::endl;
  //treef->Draw("Theta2>>ModulationSimple(80, 0., 6.28)");
  //treef->Draw("Theta2>>ModulationCut(80, 0., 6.28)",myCut);

  treef->Draw("Theta2>>ModulationSimple(720)");
  treef->Draw("Theta2>>ModulationCut(720)",myCut);
  //
  treef->Draw("ThetaI>>ModulationFirstStep(720)");
  //
  ModulationCutH = (TH1D*)gDirectory->FindObject("ModulationCut");
  ModulationSimpleH = (TH1D*)gDirectory->FindObject("ModulationSimple");
  //
  ModulationFirstStep = (TH1D*)gDirectory->FindObject("ModulationFirstStep");
  //
  if (ModulationDegree<1.0)
    {
      float n_non_polar =  (1.0-ModulationDegree)/ModulationDegree;

      long Nsimple0 = (long) ModulationSimpleH->GetEntries();
      long Ncut0    = (long) ModulationCutH->GetEntries();
      long Nfirst0  = (long) ModulationFirstStep->GetEntries();

      long Nsimple = (long) n_non_polar*Nsimple0;
      long Ncut    = (long) n_non_polar*Ncut0;
      long Nfirst  = (long) n_non_polar*Nfirst0;

      for(long  i = 0; i < Nsimple; i++) ModulationSimpleH->Fill(rnd->Uniform(0,2*PI));
      for(long  i = 0; i < Ncut;    i++) ModulationCutH->Fill(rnd->Uniform(0,2*PI));
      for(long  i = 0; i < Nfirst;    i++) ModulationFirstStep->Fill(rnd->Uniform(0,2*PI));
      std::cout<<"Increasing the statistics..."<<Nsimple0<<" -> "<<Nsimple<<std::endl;
      std::cout<<"Increasing the statistics..."<<Ncut0<<" -> "<<Ncut<<std::endl;
      std::cout<<"Increasing the statistics..."<<Nfirst0<<" -> "<<Nfirst<<std::endl;
      
    }
  
  fufi1 = new TF1("fufi","[0]+[1]*(cos(x+[2]))^2",0,2*PI);
  fufi1->SetParLimits(2,-PI,PI);
  fufi1->SetParLimits(0,0.,50000);
  fufi1->SetParLimits(1,0.,50000);

  ModulationSimpleH->Fit("fufi","E","",0.,2*PI);  
  ModulationCutH->Fit("fufi","E","",0.,2*PI);  
  ModulationFirstStep->Fit("fufi","E","",0.,2*PI);
 
  
  UpdateFile->cd();
  treef->Write();
  ClusterDim->Write();
  ThirdMomHi->Write();
  SecondRatioHi->Write();
  ModulationCutH->Write();
  ModulationSimpleH->Write();
  RecordClusterDimH->Write();
  RawClusterDimH->Write();
  RecordDimH->Write();
  ModulationFirstStep->Write();
}

void TTreeAnalysis::GetResults_SIMPLE(double &Mu, double &MuErr, double &Phi, double &PhiError, double &CutEfficiency, double &FitProbabilityCut)
{
  TF1 *fuficut = (TF1*) ModulationSimpleH->GetListOfFunctions()->FindObject("fufi");
  double SelectedEvents = ModulationSimpleH->GetEntries();
  
  double A2    = fuficut->GetParameter(0);
  double B2    = fuficut->GetParameter(1);
  double phi02 = fuficut->GetParameter(2);
  double Deltaphi02 = fuficut->GetParError(2);
  double TotalEvents    = ModulationSimpleH->GetEntries();

  double AError2 = fuficut->GetParError(0);
  double BError2 = fuficut->GetParError(1);
  
  std::cout<<"SIMPLE: "<<A2<<" "<<AError2<<" "<<B2<<" "<<BError2<<std::endl;

  Mu   = 100.0 * B2 / (2.*A2+B2);
  Phi = phi02*180./PI;
  PhiError = Deltaphi02*180./PI;
  MuErr= 200.0*(A2*BError2+B2*AError2)/((B2+2*A2)*(B2+2*A2));
  
  CutEfficiency=SelectedEvents/TotalEvents;
  FitProbabilityCut = fuficut->GetProb();
  std::cout<<"SIMPLE: Mu: "<<Mu<<" Mu err: "<<MuErr<<" Phi: "<<Phi<<" PhiError: "<<PhiError<<" CutEfficiency: "<< CutEfficiency <<" FitProbabilityCut: "<< FitProbabilityCut <<std::endl;
}

void TTreeAnalysis::GetResults_CUT(double &Mu, double &MuErr, double &Phi, double &PhiError, double &CutEfficiency, double &FitProbabilityCut)
{
  TF1 *fuficut = (TF1*) ModulationCutH->GetListOfFunctions()->FindObject("fufi");
  double SelectedEvents = ModulationCutH->GetEntries();
  
  double A2    = fuficut->GetParameter(0);
  double B2    = fuficut->GetParameter(1);
  double phi02 = fuficut->GetParameter(2);

  double Deltaphi02 = fuficut->GetParError(2);
  double TotalEvents    = ModulationSimpleH->GetEntries();

  double AError2 = fuficut->GetParError(0);
  double BError2 = fuficut->GetParError(1);

  std::cout<<"CUT: "<<A2<<" "<<AError2<<" "<<B2<<" "<<BError2<<std::endl;

  Mu   = 100.0 * B2 / (2.*A2+B2);
  Phi = phi02*180./PI;
  PhiError = Deltaphi02*180./PI;
  MuErr= 200.0*(A2*BError2+B2*AError2)/((B2+2*A2)*(B2+2*A2));
  
  CutEfficiency=SelectedEvents/TotalEvents;
  FitProbabilityCut = fuficut->GetProb();
  std::cout<<"Mu: "<<Mu<<" Mu err: "<<MuErr<<" Phi: "<<Phi<<" PhiError: "<<PhiError<<" CutEfficiency: "<< CutEfficiency <<" FitProbabilityCut: "<< FitProbabilityCut <<std::endl;
}

void TTreeAnalysis::Save()
{
  /*
    TFile UpdateFile(keyName,"RECREATE");
    treef->Write();
    ThirdMomHi->Write();
    SecondRatioHi->Write();
    ModulationCutH->Write();
    //  ModulationSimpleH->Write();
    UpdateFile.Close();
  */
}

void TTreeAnalysis::DeleteAll()
{
  UpdateFile->Close();
  
  //  delete t;
  //  delete treef;
  //  delete fufi1;

}

void TTreeAnalysis::ClusterViewComp(int NCluster)
{
  std::vector<ADC> Digi = ReadCluster(NCluster);
  myDetector->myDrawSignal(Digi,Readout);
  int Clusterdim = Digi.size();
  SetClusterPar(Clusterdim);
  //a = 200;
  //b= 1.3;
  //c= 2.9;
  std::cout<<" Parameters "<<a<<"--"<<b<<"--"<<c<<std::endl;
  TCluster Cluster(a,b,c);
  std::vector<Pixel> D = Cluster.ADCtoPixel(Digi,Readout);
 
  Cluster.Analysis(D);
  double ImpactRecX = Cluster.GetImpactX();
  double ImpactRecY =Cluster.GetImpactY(); 
  double X= Cluster.GetBaricenterX();
  double Y= Cluster.GetBaricenterY();
  double Theta2 = Cluster.GetTheta2();
  double Theta = Cluster.GetTheta();
  //double TM = Cluster.GetThirdMomentum();
  double MaxS =Cluster.GetMaxMom();
  //double MinS = Cluster.GetMinM();
  double RealImpactX =XIn*10000+XShift;
  double RealImpactY =YIn*10000+YShift;//0;
 
  std::pair<double,double>  realp = Cluster.Froc(RealImpactX,RealImpactY,Phi);
  std::pair<double,double>  analysisp = Cluster.Froc( ImpactRecX, ImpactRecY,Theta2);   //( ImpactRecX, ImpactRecX,Theta2);
  std::cout<< " Froc fin X =  "<<analysisp.first<<"Froc fin Y = "<<analysisp.second<<std::endl;  
  TArrow  *real=new  TArrow(RealImpactX,RealImpactY,realp.first,realp.second);
  TArrow  *analys=new  TArrow(ImpactRecX,ImpactRecY,analysisp.first,analysisp.second);
  TMarker *ImpactRec = new TMarker(ImpactRecX,ImpactRecY,28);
  TMarker *Baric = new TMarker(X,Y,28); 
  TMarker *RealI = new TMarker(RealImpactX,RealImpactY,28);
  ImpactRec->SetMarkerColor(2);//ROSSO RICOSTRUITO
  ImpactRec->SetMarkerSize(2);
  ImpactRec->Draw();
  Baric->SetMarkerColor(4);
  Baric->SetMarkerSize(2);
  Baric->Draw();
  RealI->SetMarkerColor(3);//VERDE VERO IMPATTO
  RealI->SetMarkerSize(2);
  RealI->Draw();
  real->SetLineColor(3);
  real->Draw();
  analys->SetLineColor(2);
  analys->Draw();
  //
  
  std::pair<double,double>  MaxAxdx = Cluster.Froc(X,Y,Theta);
  std::pair<double,double>  MinAxdx = Cluster.Froc(X,Y,Theta+PI/2.);
  std::pair<double,double>  MaxAxsx = Cluster.Froc(X,Y,Theta+PI);
  std::pair<double,double>  MinAxsx = Cluster.Froc(X,Y,Theta+PI*3./2);
  TLine *MinAx =new TLine(MinAxsx.first,MinAxsx.second,MinAxdx.first,MinAxdx.second);
  TLine *MaxAx =new TLine(MaxAxsx.first, MaxAxsx.second, MaxAxdx.first, MaxAxdx.second);
  MinAx->Draw();
  MaxAx->Draw();
  TText *texMa = new TText(X+200,Y+200,"Major Axis");
  TText *texma = new TText(X+250,Y+250,"Minor Axis");
  
  texMa->SetTextFont(42);
  texMa->SetTextSize(0.03);
  texMa->Draw();
  texma->SetTextFont(42);
  texma->SetTextSize(0.03);

  std::cout<<"theta0 ="<<Theta<<"gradi = "<<Theta*360/PI<<std::endl;
  //
  //Modifiche grafiche
  double LargRad = b * sqrt(MaxS);
  double SmallRad = a * sqrt(MaxS);
  TEllipse *SREll = new TEllipse(X,Y,SmallRad);
  TEllipse *LREll = new TEllipse(X,Y,LargRad);
  SREll->SetLineColor(4);
  LREll->SetLineColor(7);
  std::cout<<"  sqrt II =  "<<sqrt(MaxS)<<std::endl;
}

