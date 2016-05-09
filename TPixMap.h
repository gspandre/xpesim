#ifndef TPIXMAP_HH
#define TPIXMAP_HH

#include "MC.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "TCanvas.h"
#include "THexagon3d.h"
#include "TMarker.h"

class TPixMap
{
 public:
  TPixMap();
  ~TPixMap(){;}
  int Sampling(double x, double y);
  void Draw3D(int Color, double Size, int LineWidth);
  void Draw3DPixel(int channel,int Color, double Size);
  //  void myDraw3DPixel(double x, double y,int Color, double Size);
  inline double GetX(int channel) {return X[channel];}
  inline double GetY(int channel) {return Y[channel];}
  
  
 private:
  TCanvas *PixMapCanvas;
  ifstream PixelMapFile;
  int N_X,N_Y,N_Ch;

  double MaxX,MinX,MaxY,MinY,Pitch;

  double *X;
  double *Y;
  int *M;
  int *B;
  
};

#endif
