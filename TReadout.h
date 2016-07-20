#ifndef TREADOUT_HH
#define TREADOUT_HH

#include "MC.h"
#include "TDimension.h"
#include "THexagon3d.h"
#include "THexagon2d.h"


class TReadout 
{ 
 public:
  
  TReadout(){;}
  TReadout(TDimension *Dimension);
  ~TReadout(){;}
  
    double              Modulo(double x, double y);
    void                Setdimension(int xmin , int xmax, int ymin, int ymax);
    
    std::pair<int,int>  PositiontoMatrix(double Xa, double Ya);
    std::pair<int,int>  Position2Xpol(double Xa, double Ya);

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
    void                myDraw2D(int Color, double Size, double LineWidth);

    void                SetPitch(double p); 
    std::pair<double,double>  XYorigin(std::pair<double,double> xymod);    
    inline double       GetPitchR(){return Pitch;}       
 private:
    
    TCanvas *PixMapCanvas;
    TDimension *Dimension;
    int         IMin,IMax,JMin,JMax;
    double      Pitch, Size, Pitch_r;
    double      L,X,Y,PassoX,PassoY;
    double      Xx,I,J;
    double      Xxpol, Yxpol;
    int         j;
    double      YHexagonCenter,XHexagonCenter;
    int         Cluster,Channel,Iint,Jint,GChannel;
    
};

#endif

  
