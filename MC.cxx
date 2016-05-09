//MC simulation for Pixi detector - Author: G.Spandre (INFN Pisa)
//==============================================================

#include "MC.h"
#include "MyMainFrame.h"

MyMainFrame *theGUI;
//TRint *theApp;

int main(int argc, char**argv){
  TRint theApp("MonteCarlo", &argc, argv);
 //theApp = new TRint("Pixy", 0, 0);
  if (gROOT->IsBatch()) {
    fprintf(stderr, "%s: cannot run in batch mode\n", argv[0]);
    return 1;
  }
  theGUI = new MyMainFrame(gClient->GetRoot(),200,200);
  theApp.Run();
  gSystem->ProcessEvents();
}
