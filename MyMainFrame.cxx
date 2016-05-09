#include "MyMainFrame.h"

MyMainFrame::MyMainFrame(const TGWindow *p,UInt_t w,UInt_t h):TGMainFrame(p, w, h){
  // Create a main frame

  fMenuBarLayout = new TGLayoutHints(kLHintsTop | kLHintsExpandX);
  fMenuBarItemLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft);
  fMenuBarHelpLayout = new TGLayoutHints(kLHintsTop | kLHintsRight);
    
  fMenuFile = new TGPopupMenu(gClient->GetRoot());
  fMenuFile->AddEntry("&Open...", M_FILE_OPEN);
  fMenuFile->AddEntry("&Save", M_FILE_SAVE);
  fMenuFile->AddEntry("S&ave as...", M_FILE_SAVEAS);
  fMenuFile->AddSeparator();
  fMenuFile->AddEntry("&Print", M_FILE_PRINT);
  fMenuFile->AddEntry("P&rint setup...", M_FILE_PRINTSETUP);
  fMenuFile->AddSeparator();
  fMenuFile->AddEntry("E&xit", M_FILE_EXIT);

  fMenuFile->Connect("Activated(Int_t)", "MyMainFrame", this, "HandleMenu(Int_t)");
  //fMenuFile->Connect("PoppedUp()", "MyMainFrame", this, "HandlePopup()");
  //fMenuFile->Connect("PoppedDown()", "MyMainFrame", this, "HandlePopdown()");
  fMenuFile->HideEntry(M_FILE_SAVEAS);

  fMenuHelp = new TGPopupMenu(gClient->GetRoot());
  fMenuHelp->AddEntry("&Contents", M_HELP_CONTENTS);
  fMenuHelp->AddSeparator();
  fMenuHelp->AddEntry("&Annotation", M_HELP_ABOUT);
  fMenuHelp->Connect("Activated(Int_t)", "MyMainFrame", this, "HandleMenu(Int_t)");

  fMenuBar = new TGMenuBar(this, 1, 1, kHorizontalFrame);
  fMenuBar->AddPopup("&File", fMenuFile, fMenuBarItemLayout);
  fMenuBar->AddPopup("&Help", fMenuHelp, fMenuBarHelpLayout);
  this->AddFrame(fMenuBar, fMenuBarLayout);

  fMainFrameHor = new TGCompositeFrame(this, 150, 150, kHorizontalFrame|kRaisedFrame);
  fRightFrame = new TGCompositeFrame(fMainFrameHor, 100, 100, kVerticalFrame | kLHintsRight);
  fRightVertFrame = new TGCompositeFrame(fRightFrame, 100, 100, kVerticalFrame);
  fLeftFrame = new TGCompositeFrame(fMainFrameHor, 100, 100, kVerticalFrame | kRaisedFrame);
  fLeftVertFrame = new TGCompositeFrame(fLeftFrame, 100, 100, kVerticalFrame);

 // Create widget for input
  fLeftFrame->AddFrame(fLeftVertFrame, new TGLayoutHints(kLHintsLeft,2,0,4,2));
  fMainFrameHor->AddFrame(fLeftFrame, new TGLayoutHints(kLHintsLeft,5,2,5,5));
 
  string tarray[] = {"Events:","Energy (keV) from:","Energy (keV) to:","Gas Id:","Pressure (atm):","Thickness (cm):", "Source:", 
		     "Pol. angle (deg):", "Pol. degree:"};
  string Invarray[] = {"1000","2","12","5","1.0","1.0","monochromatic","30.","1."};

  if ((!gSystem->AccessPathName("info.dat", kFileExists)))ReadConfigFile();
  else {
     for (int j = 0; j < 9; j++) {
       varray[j] = Invarray[j];
     }
  }
 
  // Left frame
  // 2 column, n rows for Input data
  fF1 = new TGGroupFrame(fLeftVertFrame, "Input", kVerticalFrame);  
  fF1->SetTitlePos(TGGroupFrame::kRight); // right aligned
  fLeftVertFrame->AddFrame(fF1, new TGLayoutHints(kLHintsTop| kLHintsExpandX| kLHintsLeft, 5, 5, 5, 5));

  tlo = new TGTableLayout(fF1, 10, 2);
  fF1->SetLayoutManager(tlo);

  for (int j = 0; j < 9; j++) {
    TString buff = tarray[j];
    tloh = new TGTableLayoutHints(0,1,j,j+1,
                                    kLHintsExpandX|kLHintsExpandY |
                                  kLHintsShrinkX|kLHintsShrinkY );
                       
    fF1->AddFrame(new TGLabel(fF1, new TGHotString(buff)),tloh);     

    tbuf = new TGTextBuffer(10);
    TString vbuff = varray[j];
    tbuf->AddText(0, vbuff);
    entry.push_back(tbuf);
    tent  = new TGTextEntry(fF1, tbuf); 
    tent->Resize(145, tent->GetDefaultHeight());
    //tent->SetFont("-adobe-courier-bold-r-*-*-12-*-*-*-*-*-iso8859-1");
    tloh = new TGTableLayoutHints(1,2,j,j+1,
                                    kLHintsExpandX|kLHintsExpandY |
				  kLHintsShrinkX|kLHintsShrinkY );
                                
    fF1->AddFrame(tent,tloh);
  }
  
  TString buff = "Window:";
  tloh = new TGTableLayoutHints(0,1,9,10,
                                    kLHintsExpandX|kLHintsExpandY |
				kLHintsShrinkX|kLHintsShrinkY );                       
  fF1->AddFrame(new TGLabel(fF1, new TGHotString(buff)),tloh); 
  fCombo = new TGComboBox(fF1, 110);
  tloh = new TGTableLayoutHints(1,2,9,10,
                                    kLHintsExpandX|kLHintsExpandY |
				kLHintsShrinkX|kLHintsShrinkY );

  TString win[] = {"Poly 1um - Al 0.06mu - grid", "Be 50um", "Be 150um"};
  for (int i = 0; i < 3; i++) {
    fCombo->AddEntry(win[i], i+1);
  }
  //printf("%c",230);
  fCombo->Resize(145, 20);
  fCombo->Select(2);
  fF1->AddFrame(fCombo,tloh);
  fF1->Layout();
  fCombo->Connect("Selected(int)", "MyMainFrame", this, "SelectCombo(int)"); 

  Pixel_t green;
  fClient->GetColorByName("green", green);
  fF3 = new TGGroupFrame(fLeftVertFrame, "Run Simulation", kVerticalFrame);
  fLeftVertFrame->AddFrame(fF3, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5));
  fF3->SetLayoutManager(new TGMatrixLayout(fF3, 0, 2, 10));

  fF3->AddFrame(runMC = new TGTextButton(fF3,"Single run")); 
  runMC->ChangeBackground(green);
  runMC->Connect("Clicked()","MyMainFrame",this,"runSimulation()");
  runMC->Resize(80, runMC->GetDefaultHeight());  

  fF3->AddFrame(runScan = new TGTextButton(fF3,"Energy scan")); 
  runScan->Connect("Clicked()","MyMainFrame",this,"runEnergyScan()");
  runScan->Resize(80, runScan->GetDefaultHeight());
  fF3->Resize();

  int entrySize = 110;
  fF2 = new TGGroupFrame(fLeftVertFrame, "Draw", kVerticalFrame);
  fLeftVertFrame->AddFrame(fF2, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5));
  fF2->SetLayoutManager(new TGMatrixLayout(fF2, 0, 3, 10));

  fF2->AddFrame(drawEff = new TGTextButton(fF2,"Gas Efficiency"));
  drawEff->Connect("Clicked()","MyMainFrame",this,"DrawEff()");
  drawEff->Resize(entrySize, drawEff->GetDefaultHeight());  

  fF2->AddFrame(drawRange = new TGTextButton(fF2,"Elec. Range"));
  drawRange->Connect("Clicked()","MyMainFrame",this,"DrawElectronRange()");
  drawRange->Resize(entrySize, drawRange->GetDefaultHeight()); 

  fF2->AddFrame(drawPheX = new TGTextButton(fF2,"PhEl. Xsection"));
  drawPheX->Connect("Clicked()","MyMainFrame",this,"DrawPhElX()");
  drawPheX->Resize(entrySize, drawPheX->GetDefaultHeight()); 
	     
  fF2->AddFrame(drawAbsL = new TGTextButton(fF2,"Absorption Lenght"));
  drawAbsL->Connect("Clicked()","MyMainFrame",this,"DrawAbsLenght()");
  drawAbsL->Resize(entrySize, drawAbsL->GetDefaultHeight()); 

  fF2->AddFrame(drawAbsPr = new TGTextButton(fF2,"Absorption Prob."));
  drawAbsPr->Connect("Clicked()","MyMainFrame",this,"DrawAbsProb()");
  drawAbsPr->Resize(entrySize, drawAbsPr->GetDefaultHeight()); 
  fF2->Resize();

  // Progress bar for event generation and analysis  
  fHProg = new TGHProgressBar(fLeftVertFrame, TGProgressBar::kFancy, 250);
  fHProg->SetBarColor("purple");
  fHProg->ShowPosition();
  fLeftVertFrame->AddFrame(fHProg, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 5, 5,  5, 10));
  
  // Right frame
  fEcanvas = new TRootEmbeddedCanvas("Ecanvas",fRightVertFrame,610,600);
  fRightVertFrame->AddFrame(fEcanvas, new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,2,0,0,1));
 
  fCanvas = fEcanvas->GetCanvas();
  fCanvas->Connect("ProcessedEvent(int,int,int,TObject*)", "MyMainFrame", 
		this, "ExecEvent(int,int,int,TObject*)");

  status = new TGStatusBar(fRightVertFrame,510,100);
  status->SetParts(3);
  fRightVertFrame->AddFrame(status, new TGLayoutHints(kLHintsLeft,2,2,0,5));
  fRightFrame->AddFrame(fRightVertFrame, new TGLayoutHints(kLHintsLeft,2,0,4,2));

  fMainFrameHor->AddFrame(fRightFrame, new TGLayoutHints(kLHintsExpandX,2,0,2,2));

  // Create a horizontal for messages
  hframe = new TGCompositeFrame(this,1,1, kHorizontalFrame | kRaisedFrame);
  fMsg = new TGTextEntry(hframe, fTbMsg = new TGTextBuffer(100));
  fTbMsg->AddText(0, "Message window");
  hframe->AddFrame(fMsg, new TGLayoutHints(kLHintsExpandX|kLHintsCenterX| kLHintsCenterY,10,10,0,0));


  this->AddFrame(fMainFrameHor, new TGLayoutHints(kLHintsCenterX|kLHintsExpandX ,2,2,2,2));
  this->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX|kLHintsExpandX,2,2,2,2));
  // Set a name to the main frame
  this->SetWindowName("Pixie simulation GUI");
  // Map all subwindows of main frame
  this->MapSubwindows();
  this->Resize();
  // Map main frame
  this->MapWindow();

  Experiment = 0;
  rnd = 0;
  trasparency = 0;
  SelectCombo(1);
}


void MyMainFrame::runSimulation(){
  ReadInput();
  InitExp();
  int print = 1;
  cout << "Initialize experiment ==>>  Source:" << Source << " Energy: " << EnI << " keV" << endl;
  if(Source == "monochromatic")Experiment->SetSource(polAngle,polDeg,EnI);
  cout << "Num Events: " << Events <<endl;
  Experiment->EventsTree(Events);
  fCanvas = fEcanvas->GetCanvas();
  fCanvas->Clear();
  fCanvas->Divide(2,2);
  fCanvas->cd(1);
  (Experiment->pheRange)->Draw(); 
  fCanvas->cd(2);
  (Experiment->augRange)->Draw();
  fCanvas->cd(3);
  (Experiment->pheTrueRange)->Draw();  
  fCanvas->cd(4);
  (Experiment->augTrueRange)->Draw(); 
  fCanvas->Update();
 
  if(print) {
    TRegexp re("Scan");
    TString name = Experiment->GetNameforFile();
    name(re) = "Erange-";
    name += "_";
    name += EnI;
    name += "keVphoton.png";
    TImage *img = TImage::Create();
    img->FromPad(fCanvas);
    img->WriteImage(name);
    delete img;
    cout << name << endl;
    cout << "Photoelectron practical range: " << Experiment->pheRange->GetMean() << endl;
    cout << "Auger electron practical range: " << Experiment->augRange->GetMean() << endl;
    cout << "Photoelectron total range: " << Experiment->pheTrueRange->GetMean() << endl;
    cout << "Auger electron total range: " << Experiment->augTrueRange->GetMean() << endl;
  }
}

void MyMainFrame::runEnergyScan(){
  ReadInput();
  InitExp();
  fHProg->Reset();
  fHProg->SetRange(EnI,EnF+1.);
  if(Source == "monochromatic") cout << "Initialize experiment ==>>  Source:" << Source << endl;
  float ScanEnergy = EnI;
  float DeltaEn1 = 0.5;
  float DeltaEn2 = 1.0;
  float Delta;
  while(ScanEnergy <= EnF){
    if(Source == "monochromatic")cout << " Energy: " << ScanEnergy << " keV" << endl;
    Experiment->SetSource(polAngle,polDeg,ScanEnergy);
    Experiment->EventsTree(Events);
    if (ScanEnergy < 10.) Delta = DeltaEn1;
    else Delta = DeltaEn2;
    ScanEnergy+=Delta;
    fHProg->Increment(Delta);
  }
}


void MyMainFrame::SelectCombo(int id){
  TString Dummy;
  ifstream in0;
  ifstream in1;
  ifstream in2;
  switch (id) {
  case 1: 
    Energy.clear();
    Trasp.clear();
    in0.open("WINDOWS/Window_Poly1um_Al006um_Grid018.txt");
    for (int j=0; j<3; j++)in0 >> Dummy;
    while(1){
      in0 >> En >> Tr;
      if(!in0.good())break;
      Energy.push_back(En);
      Trasp.push_back(Tr);
     }
    nlines = (int)Energy.size();
    fCanvas = fEcanvas->GetCanvas();
    fCanvas->Clear();
    if(trasparency) delete trasparency;
    trasparency = 0;
    trasparency = new TGraph(nlines,&Energy[0],&Trasp[0]);
    trasparency->SetMarkerStyle(2);
    trasparency->SetMarkerSize(0.7);
    trasparency->SetMarkerColor(2);
    trasparency->GetYaxis()->SetTitleSize(0.035);
    trasparency->GetXaxis()->SetTitleSize(0.035);
    trasparency->GetYaxis()->SetTitleOffset(1.3);
    trasparency->GetXaxis()->SetTitleOffset(1.2);
    trasparency->GetYaxis()->SetLabelSize(0.035);
    trasparency->GetXaxis()->SetLabelSize(0.035);
    trasparency->GetXaxis()->SetTitle("Electron energy (keV)");
    trasparency->GetYaxis()->SetTitle("Window Trasmission");
    trasparency->SetTitle("Polyimide 1 micron- Al 0.06 micron Grid 0.18 micron");
    fCanvas->cd();
    trasparency->Draw("APL");
    fCanvas->Update();
    in0.close();
    break;

  case 2: 
    Energy.clear();
    Trasp.clear();

    in1.open("WINDOWS/Window_Be50um.txt");
    while(1){
      in1 >> En >> Tr;
      if(!in1.good())break;
      Energy.push_back(En);
      Trasp.push_back(Tr);
    }
    nlines = (int)Energy.size();
    fCanvas = fEcanvas->GetCanvas();
    fCanvas->Clear();
    if(trasparency) delete trasparency;
    trasparency = 0;
    trasparency = new TGraph(nlines,&Energy[0],&Trasp[0]);
    trasparency->SetMarkerStyle(2);
    trasparency->SetMarkerSize(0.7);
    trasparency->SetMarkerColor(2);
    trasparency->GetYaxis()->SetTitleSize(0.035);
    trasparency->GetXaxis()->SetTitleSize(0.035);
    trasparency->GetYaxis()->SetTitleOffset(1.3);
    trasparency->GetXaxis()->SetTitleOffset(1.2);
    trasparency->GetYaxis()->SetLabelSize(0.035);
    trasparency->GetXaxis()->SetLabelSize(0.035);
    trasparency->GetXaxis()->SetTitle("Electron energy (keV)");
    trasparency->GetYaxis()->SetTitle("Window Trasmission");
    trasparency->SetTitle("Be 50 micron ");
    fCanvas->cd();
    trasparency->Draw("APL");
    fCanvas->Update();
    in1.close();
    break;

  case 3: 
    Energy.clear();
    Trasp.clear();
    in2.open("WINDOWS/Window_Be150um.txt");
    while(1){
      in2 >> En >> Tr;
      if(!in2.good())break;
      Energy.push_back(En);
      Trasp.push_back(Tr);
    }
    nlines = (int)Energy.size();
    fCanvas = fEcanvas->GetCanvas();
    fCanvas->Clear();
    if(trasparency) delete trasparency;
    trasparency = 0;
    trasparency = new TGraph(nlines,&Energy[0],&Trasp[0]);
    trasparency->SetMarkerStyle(2);
    trasparency->SetMarkerSize(0.7);
    trasparency->SetMarkerColor(2);
    trasparency->GetYaxis()->SetTitleSize(0.035);
    trasparency->GetXaxis()->SetTitleSize(0.035);
    trasparency->GetYaxis()->SetTitleOffset(1.3);
    trasparency->GetXaxis()->SetTitleOffset(1.2);
    trasparency->GetYaxis()->SetLabelSize(0.035);
    trasparency->GetXaxis()->SetLabelSize(0.035);
    trasparency->GetXaxis()->SetTitle("Electron energy (keV)");
    trasparency->GetYaxis()->SetTitle("Window Trasmission");
    trasparency->SetTitle("Be 150 micron ");
    fCanvas->cd();
    trasparency->Draw("APL");
    fCanvas->Update();
    in2.close();
    break;
  }
}


void MyMainFrame::ExecEvent(int event,int x,int y,TObject* selected){
  const Int_t kTMAX=256;
  static char atext[kTMAX];  
  double xc = gPad->PadtoX(gPad->AbsPixeltoX(x));
  double yc = gPad->PadtoY(gPad->AbsPixeltoY(y));
  sprintf(atext, "x: %f - y: %f", xc, yc); 
  status->SetText(atext,0); 
}


TGraph *MyMainFrame::GasEff(double E1,double E2)
{
  double Eff;
  int npts;
  std::vector<Double_t> Efficiency;
  std::vector<Double_t> Energy;
  npts = 0;

  double En = E1;
  while (En <E2){
    Eff = Experiment->GetEfficiencyMixture(En);
    Efficiency.push_back(Eff);
    Energy.push_back(En);
    En+=0.5;
    npts++;
  }
  TGraph *Effgr = new TGraph(npts,&Energy[0],&Efficiency[0]);
  Effgr->SetMarkerStyle(2);
  Effgr->SetMarkerSize(0.7);
  Effgr->SetMarkerColor(2);
  Effgr->GetYaxis()->SetTitleSize(0.035);
  Effgr->GetXaxis()->SetTitleSize(0.035);
  Effgr->GetYaxis()->SetTitleOffset(1.4);
  Effgr->GetXaxis()->SetTitleOffset(1.2);
  Effgr->GetYaxis()->SetLabelSize(0.035);
  Effgr->GetXaxis()->SetLabelSize(0.035);
  Effgr->GetXaxis()->SetTitle("Electron energy (keV)");
  Effgr->GetYaxis()->SetTitle("Gas Efficiency");
  WriteInfoBar();
  return Effgr;
}

void MyMainFrame::DrawEff()
{
  ReadInput();
  InitExp();
 
  gPad->SetLogy(1);
  gPad->SetLogx(0);
  fCanvas = fEcanvas->GetCanvas();
  fCanvas->Clear();
  
  TGraph *eff = GasEff(EnI,EnF);
  fCanvas->cd();
  eff->Draw("APL");
  fCanvas->Update();
}

void MyMainFrame::DrawElectronRange(){
  ReadInput();
  InitExp();  
  gPad->SetLogy(1);
  gPad->SetLogx(0);
  fCanvas = fEcanvas->GetCanvas();
  fCanvas->Clear();
  int npts;
  std::vector<Double_t> Range;
  std::vector<Double_t> Energy;
  npts = 0;
  double En = EnI;
  while (En <EnF){
    Energy.push_back(En);
    Range.push_back(Experiment->GetElectronRange(En));
    En+=0.5;
    npts++;
  }
  TGraph *Erange = new TGraph(npts,&Energy[0],&Range[0]);
  Erange->SetMarkerStyle(2);
  Erange->SetMarkerSize(0.7);
  Erange->SetMarkerColor(2);
  Erange->GetYaxis()->SetTitleSize(0.035);
  Erange->GetXaxis()->SetTitleSize(0.035);
  Erange->GetYaxis()->SetTitleOffset(1.3);
  Erange->GetXaxis()->SetTitleOffset(1.2);
  Erange->GetYaxis()->SetLabelSize(0.035);
  Erange->GetXaxis()->SetLabelSize(0.035);
  Erange->GetXaxis()->SetTitle("Electron energy (keV)");
  Erange->GetYaxis()->SetTitle("Pratical range (mm)");
  fCanvas->cd();
  Erange->Draw("APL");
  fCanvas->Update();
  WriteInfoBar();
}

void MyMainFrame::DrawPhElX(){
  ReadInput();
  InitExp();  
  gPad->SetLogy(1);
  gPad->SetLogx(0);
  fCanvas = fEcanvas->GetCanvas();
  fCanvas->Clear();
  int npts;
  std::vector<Double_t> PhElXsec;
  std::vector<Double_t> Energy;
  npts = 0;
  double En = EnI;
  while (En <EnF){
    Energy.push_back(En);
    PhElXsec.push_back(Experiment->GetPhotoelectricCrossSection(En));
    En+=0.5;
    npts++;
  }
 TGraph *Xsection = new TGraph(npts,&Energy[0],&PhElXsec[0]);
 Xsection->SetMarkerStyle(2);
 Xsection->SetMarkerSize(0.7);
 Xsection->SetMarkerColor(2);
 Xsection->GetYaxis()->SetTitleSize(0.035);
 Xsection->GetXaxis()->SetTitleSize(0.035);
 Xsection->GetYaxis()->SetTitleOffset(1.3);
 Xsection->GetXaxis()->SetTitleOffset(1.2);
 Xsection->GetYaxis()->SetLabelSize(0.035);
 Xsection->GetXaxis()->SetLabelSize(0.035);
 Xsection->GetXaxis()->SetTitle("Photon energy (keV)");
 Xsection->GetYaxis()->SetTitle("Photoelectric cross section (cm^2/g)");
 fCanvas->cd();
 Xsection->Draw("APL");
 fCanvas->Update();
 WriteInfoBar();
}

void MyMainFrame::DrawAbsLenght(){
  ReadInput();
  InitExp(); 
  gPad->SetLogy(1);
  gPad->SetLogx(0);
  fCanvas = fEcanvas->GetCanvas();
  fCanvas->Clear();
  int npts;
  std::vector<Double_t> AbsLen;
  std::vector<Double_t> Energy;
  npts = 0;
  double En = EnI;
  while (En <EnF){
    Energy.push_back(En);
    AbsLen.push_back(Experiment->GetAbsorptionLenght(En));
    En+=0.5;
    npts++;
  }
  TGraph *AbsLenght = new TGraph(npts,&Energy[0],&AbsLen[0]);
  AbsLenght->SetMarkerStyle(2);
  AbsLenght->SetMarkerSize(0.7);
  AbsLenght->SetMarkerColor(2);
  AbsLenght->GetYaxis()->SetTitleSize(0.035);
  AbsLenght->GetXaxis()->SetTitleSize(0.035);
  AbsLenght->GetYaxis()->SetTitleOffset(1.3);
  AbsLenght->GetXaxis()->SetTitleOffset(1.2);
  AbsLenght->GetYaxis()->SetLabelSize(0.035);
  AbsLenght->GetXaxis()->SetLabelSize(0.035);
  AbsLenght->GetXaxis()->SetTitle("Photon energy (keV)");
  AbsLenght->GetYaxis()->SetTitle("Absorption lenght (cm)");
  fCanvas->cd();
  AbsLenght->Draw("APL");
  fCanvas->Update();
  WriteInfoBar();

  fHProg->Reset();
  int cnt = 0;
  int inc = 1;
  while (cnt < 100){   
    if (cnt < 100) {
     cnt += inc;
     fHProg->Increment(inc);
   }
  }
}

void MyMainFrame::DrawAbsProb(){
  ReadInput();
  InitExp(); 
  gPad->SetLogy(1);
  gPad->SetLogx(1);
  fCanvas = fEcanvas->GetCanvas();
  fCanvas->Clear();
  std::vector<TGraph*>  AbsProb = Experiment->GetAbsProbGraph(EnI);
  AbsProb[0]->GetYaxis()->SetTitleSize(0.035);
  AbsProb[0]->GetXaxis()->SetTitleSize(0.035);
  AbsProb[0]->GetYaxis()->SetTitleOffset(1.3);
  AbsProb[0]->GetXaxis()->SetTitleOffset(1.2);
  AbsProb[0]->GetYaxis()->SetLabelSize(0.035);
  AbsProb[0]->GetXaxis()->SetLabelSize(0.035);
  AbsProb[0]->Draw("apl");
  int size = AbsProb.size();
  for(int i=1; i<size;i++) {
    AbsProb[i]->Draw("pl");
   }
  fCanvas->Update();
  WriteInfoBar();
}

void MyMainFrame::WriteInfoBar(){
  fTbMsg->Clear();
  TString message = Experiment->GetMixType();
  fTbMsg->AddText(0,message); 
  gClient->NeedRedraw(fMsg);
}

void MyMainFrame::InitExp(){
  if(rnd) {
    delete rnd;
    rnd = 0;
  }
  rnd = new TRandom3();
  
  if(Experiment){
    delete Experiment;
    Experiment = 0;
  }
  Experiment = new TExperiment(rnd);
  Experiment->SetMixID(MixId);
  Experiment->SetPressure(Pressure);
  Experiment->SetThickness(Thickness);
  Experiment->Generalsetup();  
}


void MyMainFrame::ReadInput(){
  Events = atoi(((TGTextBuffer*)entry[0])->GetString());
  EnI = atof(((TGTextBuffer*)entry[1])->GetString());
  EnF = atof(((TGTextBuffer*)entry[2])->GetString());
  MixId = atoi(((TGTextBuffer*)entry[3])->GetString());
  Pressure = atof(((TGTextBuffer*)entry[4])->GetString());
  Thickness = atof(((TGTextBuffer*)entry[5])->GetString());
  Source = ((TGTextBuffer*)entry[6])->GetString();
  polAngle = atof(((TGTextBuffer*)entry[7])->GetString());
  polDeg  = atof(((TGTextBuffer*)entry[8])->GetString());		  
}

void MyMainFrame::ReadConfigFile(){
  infofileIn.open("info.dat",ios::in);
  for (int j = 0; j < 9; j++) {
    infofileIn >> varray[j] >> varray[j] >> varray[j];
  }
  infofileIn.close();
};

void MyMainFrame::WriteConfigFile(){
  
  infofileOut.open("info.dat",ios::out);
  infofileOut << "Events num: " << ((TGTextBuffer*)entry[0])->GetString() << endl;
  infofileOut << "EnergyFrom (keV): " << ((TGTextBuffer*)entry[1])->GetString() << endl;
  infofileOut << "EnergyTO (kev): " << ((TGTextBuffer*)entry[2])->GetString() << endl;
  infofileOut << "Gas ID: " << ((TGTextBuffer*)entry[3])->GetString() << endl;
  infofileOut << "Pressure (atm): " << ((TGTextBuffer*)entry[4])->GetString() << endl;
  infofileOut << "Thickness (cm): " << ((TGTextBuffer*)entry[5])->GetString() << endl;
  infofileOut << "Source type: " << ((TGTextBuffer*)entry[6])->GetString() << endl;
  infofileOut << "Polarization angle: " << ((TGTextBuffer*)entry[7])->GetString() << endl;
  infofileOut << "Polarization degree: " << ((TGTextBuffer*)entry[8])->GetString() << endl;
  infofileOut.close();
};

void MyMainFrame::HandleMenu(Int_t id)
{
  switch (id) {
  case M_FILE_OPEN:
    {
      static TString dir(".");
      fi.fIniDir    = StrDup(dir);
      printf("fIniDir = %s\n", fi.fIniDir);
      new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
      printf("Open file: %s (dir: %s)\n", fi.fFilename, fi.fIniDir);
      dir = fi.fIniDir;
    }
    break;
    
  case M_FILE_SAVE:
    printf("M_FILE_SAVE\n");
    break;
    
  case M_FILE_EXIT:
    WriteConfigFile();
    gApplication->Terminate(0);
    break;
   
  case M_HELP_ABOUT:
    ed = new Editor(this, 500, 250);
    ed->LoadFile("note.txt",0,-1);
    ed->Popup();
    break;
  }
}

MyMainFrame::~MyMainFrame()
{
  gApplication->Terminate(0);
}

ClassImp(MyMainFrame)
