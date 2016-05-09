#ifndef TDIMENSION_HH
#define TDIMENSION_HH

#include "MC.h"
class TDimension

{ 
 public:
  
  TDimension(TString DimensionFile = "Dimensions.txt");
  void ReadDimensionFromFile(TString FileName);
  void WriteDimensionFile(TString FileName);
  inline double                    GetZ_Drift() {return Z_Drift;}
  inline double                    GetZ_Gem() {return Z_Gem;};
  inline double                    GetZ_Readout(){return Z_ReadOut;}
  inline double                    GetGem_Pitch(){return  Gem_Pitch;}
  inline double                    GetGem_Radius(){return Gem_Radius;}
  inline double                    GetPitch(){return Pitch;}
  inline std::pair<int,int>  GetPixel(){return Pixel;}
  inline double                    GetPressure(){return Pressure;}
  inline void                      SetZ_Drift(double Zd) {Z_Drift=Zd;}
  inline void                      SetPressure(double PRESSURE) {Pressure = PRESSURE;}
 
  inline void                      SetZ_Gem(double Zg) {Z_Gem=Zg;}
  inline void                      SetZ_ReadOut(double Zr) {Z_ReadOut=Zr;}
  void                             SetGem_Pitch(double Gp); 
  inline void                      SetGem_Gain(double Gg) {Gem_Gain=Gg;}
  inline void                      SetGem_Radius(double Gr) {Gem_Radius=Gr;}
  inline void                      SetPitch(double P) {Pitch=P;}
  void                             SetPixel(int ,int);
  inline double                    GetGain() {return Gain;}
  inline double                    GetPac() {return Pac;}
  inline double                    GetElectronsADC() {return ElectronADC;}
  inline double                    GetHole() {return Hole;}
 
 
 private:
  double                    Pac;
  double                    Z_Drift;
  double                    Z_Gem;
  double                    Z_ReadOut;
  double                    Gem_Pitch;
  double                    Gem_Gain;
  double                    Gem_Radius;
  double                    Pitch;
  double                    Hole;
  std::pair<int,int>  Pixel;       
  double                    Gain;
  double                    ElectronADC;
  double                    Pressure;
};

#endif
