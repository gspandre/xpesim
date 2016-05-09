#include "TDimension.h"

TDimension::TDimension(TString DimensionFile) 
{
  Z_Drift     = 1.6; // cm
  Z_Gem       = 0.6;  //cm
  Z_ReadOut   = 0.0;  //cm 
  Gem_Radius  = 1.0; 
  Gem_Pitch   = 50.; //mu
  Pitch       = 50.; // mu  
  Hole        = 10;//mu
  Pixel= make_pair(352,300);//(160,138);//(160,138)(352,300)=XMax,YMax=JMax,IMax buono per Pixmap da 80
  Gain = 360;//1000
  Pac = 20;//40;
  ElectronADC=60.;

 
  Pressure = 1.0;
  ReadDimensionFromFile(DimensionFile);
}  

void TDimension::ReadDimensionFromFile(TString FileName)
{
  std::ifstream FileIn;
  FileIn.open(FileName,std::ios::in);
  if(!FileIn.is_open())
    {
      std::cout<<"Error opening "<<FileName<<std::endl;
    }
  std::cout<<" Reading parameters from file: "<<FileName<<std::endl;
  char dummy[20];
  for (int i=0;i<6;i++)
    FileIn>>dummy;
  double Thickness;
  FileIn >> Thickness;
  Z_Drift     =  Thickness+Z_Gem; // cm1.6 valori precedenti
  FileIn >> Gem_Radius;
  FileIn >> Gem_Pitch;
  FileIn >> Pitch;
  
  if(Pitch==80)
    Pixel= make_pair(160,138);
  else 
    {
      Pixel= make_pair(352,300);//Pixel= make_pair(160,138);
    }
  if (Gem_Pitch == 90)
    {
      Hole = 30;
    }
  else Hole = 10;
  
  FileIn >> Gain;
  FileIn >> ElectronADC;
  FileIn >> Pressure;
  FileIn.close();
}
void TDimension::SetGem_Pitch(double Gp)
{
  Gem_Pitch=Gp;
  if (Gem_Pitch == 90)
    {
      Hole = 30;
    }
  else 
    {
      Hole = 10;
    }
}
void TDimension::WriteDimensionFile(TString FileName)
{
  std::ofstream FileOut(FileName);
  if(!FileOut.is_open())
    {
      std::cout<<"Error opening "<<FileName<<std::endl;
    }
  std::cout<<" Writing parameters to file: "<<FileName<<std::endl;
  
  FileOut<<"Thickness Gem_Radius Gem_Pitch Pitch Gain ElectronADC"<<std::endl;
  double Thickness;
  Thickness=Z_Drift-Z_Gem;
  FileOut << Thickness<<" "<<Gem_Radius<<" "<<Gem_Pitch<<" "<<Pitch<<" "<<Gain<<" "<<ElectronADC<<" "<<Pressure<<std::endl;
  FileOut.close();
}

void TDimension::SetPixel(int Px,int Py)

{
  Pixel= make_pair(Px,Py);
}
