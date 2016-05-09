#ifndef __MyMainFrame__
#define __MyMainFrame__

#include "MC.h"
#include "TTreeAnalysis.h"
#include "TExperiment.h"
#include "TEditor.h"
#include "TImage.h"
#include "TRegexp.h"

enum ETestCommandIdentifiers {
 M_FILE_OPEN,
 M_FILE_SAVE,
 M_FILE_SAVEAS,
 M_FILE_PRINT,
 M_FILE_PRINTSETUP,
 M_FILE_EXIT,

 M__NUMBERENTRY,

 M_HELP_CONTENTS,
 M_HELP_ABOUT,
};
 

class TGWindow;
class TGMainFrame;
class TRootEmbeddedCanvas;
class MyMainFrame : public TGMainFrame {
private:

  ifstream infofileIn;
  ofstream infofileOut;

  double polAngle, polDeg, eSource;
  double En, Tr;
  int nlines;
  int Events, MixId;
  float EnI, EnF, Pressure, Thickness;

  string varray[9];
  vector <TGTextBuffer*> entry;
  vector <TGTextBuffer*>::iterator pos;
  vector <double> Energy;
  vector <double> Trasp;

  TExperiment *Experiment;
  Editor *ed;

  TCanvas *fCanvas;
  TRootEmbeddedCanvas *fEcanvas;

  TGFileInfo fi;
  TString Source;
  TRandom3 *rnd;

  TGLayoutHints *fMenuBarLayout, *fMenuBarItemLayout, *fMenuBarHelpLayout;
  TGCompositeFrame *fMainFrameHor,*fRightFrame,*fLeftVertFrame;
  TGCompositeFrame *fRightVertFrame;
  TGCompositeFrame *fLeftFrame, *hframe;
  TGGroupFrame *fF1,*fF2,*fF3;

  TGTableLayout* tlo;
  TGTableLayoutHints* tloh;

  TGPopupMenu *fMenuFile, *fMenuHelp;
  TGMenuBar *fMenuBar;
  TGHProgressBar *fHProg;
  TGStatusBar *status;

  TGTextEntry *tent; 
  TGTextEntry *fMsg;

  TGTextButton *drawEff, *drawRange, *drawAbsL, *drawPheX, *drawAbsPr;
  TGTextButton *runMC, *runScan;
  TGTextBuffer *fTbMsg;
  TGTextBuffer *tbuf;

  TGComboBox *fCombo;

  TGraph *trasparency;

 public:
  MyMainFrame(const TGWindow *p,UInt_t w,UInt_t h);
  virtual ~MyMainFrame();
  void HandleMenu(Int_t id);
  void HandlePopup() { printf("menu popped up\n"); }
  void HandlePopdown() { printf("menu popped down\n"); }
  void DrawEff();
  void DrawElectronRange();
  void DrawPhElX();
  void DrawAbsLenght();
  void DrawAbsProb();
  void ReadConfigFile();
  void WriteConfigFile();
  void ReadInput();
  void InitExp();
  void WriteInfoBar();
  void ExecEvent(int,int,int,TObject*);
  void SelectCombo(int);
  void runSimulation();
  void runEnergyScan();
  TGraph* GasEff(double, double);
  ClassDef(MyMainFrame,1);
};

#endif
