#ifndef TCLUSTER_HH
#define TCLUSTER_HH
#include <vector>
#include "TReadout.h"
#include "TDetector.h"


#define kPI TMath::Pi()


struct Pixel
{
  double X;
  double Y;
  double Charge;
};
struct Square
{
  int I;
  int J;
  double Charge;
};

class TCluster

{
 public:
  TCluster(double a =1.5,double b=2.5,double c=50000.);
  ~TCluster(){;}
  
  std::pair<double,double>                Baricenter(std::vector<Pixel> VectorIn );
  double                                 *PhiMajorAxis(std::vector<Pixel> VectorIn,double Xb,double Yb);
  double                                  ThirdMomentum(std::vector<Pixel> VectorIn,double Xb,double Yb, double T);
  std::pair<double,double>                ImpactPoint(std::vector<Pixel> VectorIn, double T,double Xb, double Yb, double ThirdM, double SecondMom);
  std::vector<Pixel>                     ADCtoPixel (std::vector<ADC> VectorIn, TReadout Readout);
  
  std::vector<Pixel>                     ADCtoPixelinRange (std::vector<ADC> VectorIn,std::vector<double> RangeV, TReadout Readout);
   std::vector<double>                  ClusterRecord(std::vector<ADC> VectorIn, TReadout Readout,int  TriggerLevel=1200, int  Width = 10);
  
  std::vector<Pixel>                     RotateVector(std::vector<Pixel> VectorIn,double cosROT,double sinROT , double XOrigin, double YOrigin);
  inline double                          GetBaricenterX() {return X;}  
  inline double                          GetBaricenterY() {return Y;}
  inline double                          GetImpactX() {return ImpactX;}
  inline double                          GetImpactY() {return ImpactY;}
  std::vector<Pixel>                     PartialWeight(std::vector<Pixel> VectorIn, double ImpX, double ImpY);
  void                                   Analysis( std::vector <Pixel>  V);
  inline double                          Getbaricenter2X() {return X2;}
  inline double                          Getbaricenter2Y() {return Y2;}
  inline double                          GetTheta2() {return Theta2;}
  inline double                          GetTheta() {return Theta;}
  inline double                          GetThirdMomentum() {return Third;}
  inline double                          GetMaxMom() {return  MaxM;} 
  inline double                          GetMinM() {return MinM;}
  inline int                             GetNStatistic() {return  NStatistic;} 
  inline double                          GetMaxSw() {return MaxSw;}
  inline double                          GetMinSw() {return MinSw;}
  std::pair<double,double>               Froc(double Xzero,double Yzero, double Angle);
  // void                                   AnalysisCentered( std::vector <Pixel> V,TReadout Readout);
  void                                   ClusterProfile(std::vector<Pixel> VectorIn,TH1D *XProfile, TH1D *YProfile);
  inline  int                            GetStatistic() {return NStatistic;}
  inline std::vector<std::pair<double,double> >   GetTriggeredPositionXY(){return TriggerPixel;}
  inline std::vector<double>     GetSquareCharge() {return ChargeSquare;}
 private:
  int     NStatistic, Pixelinrange;
  double Y,Y2;
  double X,X2;
  double Theta,Theta2;
  double Third;
  double ImpactX;
  double ImpactY;
  double MaxM,MinM;
  double SmallCircle;
  double LargeCircle;
  double WeightScale;
  double  MaxSw,MinSw;
  double Pitch;
  
  std::vector<std::pair<double,double> > TriggerPixel;
  std::vector<double> ChargeSquare;
};
#endif
