#ifndef TREADOUT_HH
#define TREADOUT_HH

#include <vector>
#include <pair.h>
#include "TH2D.h"
#include "TGaxis.h"

#include "TRandom.h"
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include "TDimension.h"
#include "THexagon3d.h"
#include "THexagon2d.h"
#include "TDimension.h"
class TReadout 
{ 
 public:
  
  TReadout(){;}
  TReadout(TDimension *Dimension);
  ~TReadout(){;}
  
    double              Modulo(double x, double y);
    void                Setdimension(int xmin , int xmax, int ymin, int ymax);
    
    std::pair<int,int>  PositiontoMatrix(double Xa, double Ya);
    
    void                CenterHexagontoXY (std::pair<int,int> ij ); 
    
    std::pair<double,double>  GetCenterHexagontoXY(); 
    void                IJtoClusterChannel (int rig,int col);
    
    std::pair<int,int>  GetClusterChannelfromIJ();
    void                ClusterChanneltoIJ (int clu,int chan);
    
    std::pair<int,int>  GetIJfromClusterChannel();
    std::pair<int,int>  GetISqJSqfromIJ(int i, int j);
    void                ChannelfromIJ(int i, int j);
    inline int          GetChannel(){return GChannel;}
    std::pair<int,int>  IJfromChannel(int ch);
    void                myDraw3DPixel(int chan, double Size);
    void                myDraw2DPixel(int chan, double Size);

    void                myDraw3D(int Color, double Size, int LineWidth);
    void                myDraw2D(int Color, double Size, int LineWidth);

    void                SetPitch(double p); 
    std::pair<double,double>  XYorigin(std::pair<double,double> xymod);    
    inline double       GetPitchR(){return Pitch;}       
 private:
    
    TCanvas *PixMapCanvas;
    TDimension *Dimension;
    int         IMin,IMax,JMin,JMax;
    double      Pitch;
    double      L,X,Y,PassoX,PassoY;
    double      Xx,I,J;
    int         j;
    double      YHexagonCenter,XHexagonCenter;
    int         Cluster,Channel,Iint,Jint,GChannel;
    
};

#endif

  
