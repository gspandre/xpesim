#include "TCluster.h"

TCluster::TCluster(double a,double b,double c)
{
  SmallCircle =a;//1.5
  LargeCircle =b;//2.5
  WeightScale =c;//50000.
  Pixelinrange = 0;
}
//trasforma un vettore di ADC(canale, carica, carica digitalizzata) in un vettore Pixel (posizioneX, posizioneY,carica)
std::vector<Pixel>   TCluster::ADCtoPixel (std::vector<ADC> VectorIn, TReadout Readout)
{
  
  std::vector<Pixel>  Cluster;
  std::vector<ADC>::iterator pos;
  for(pos = VectorIn.begin() ; pos!=VectorIn.end();++pos)
    {
      
      Pixel PixelIn;
      std::pair<int,int> IJCluster=Readout.IJfromChannel((*pos).Channel);
      Readout.CenterHexagontoXY(IJCluster);
      std::pair<double,double> XYCluster= Readout.GetCenterHexagontoXY();
      PixelIn.X=XYCluster.first;
      PixelIn.Y=XYCluster.second;
      PixelIn.Charge = TMath::Max(1400-(*pos).DigiCharge,0);//MonteCarlo
      // PixelIn.Charge = (*pos).Charge;//RealDATA
      Cluster.push_back(PixelIn);
	
    }
  Pitch = Readout.GetPitchR();
  return Cluster;
}

std::vector<Pixel>   TCluster::ADCtoPixelinRange (std::vector<ADC> VectorIn,std::vector<double> RangeV,TReadout Readout)
{
  double Xmin = RangeV[0];
  double Ymin = RangeV[1];
  double Xmax = RangeV[2];
  double Ymax = RangeV[3];
  
  std::vector<Pixel>  Cluster;
  std::vector<ADC>::iterator pos;
  for(pos = VectorIn.begin() ; pos!=VectorIn.end();++pos)
    {      
      Pixel PixelIn;
      std::pair<int,int> IJCluster=Readout.IJfromChannel((*pos).Channel);
      Readout.CenterHexagontoXY(IJCluster);
      std::pair<double,double> XYCluster= Readout.GetCenterHexagontoXY();
      PixelIn.X=XYCluster.first;
      PixelIn.Y=XYCluster.second;
      PixelIn.Charge = TMath::Max(1400-(*pos).DigiCharge,0);//(*pos).Charge
      if (PixelIn.X>Xmin && PixelIn.X<Xmax && PixelIn.Y<Ymax && PixelIn.Y>Ymin )
	{
	  Cluster.push_back(PixelIn);
	}
     }
  Pitch = Readout.GetPitchR();
  return Cluster;
}

class ConfrontaSquare {
public:
  std::pair<int,int> IJsq;
  
  ConfrontaSquare(std::pair<int,int> CC) 
    : IJsq(CC){;}  
  
  bool operator() (Square Sq)
  { 
    //std::cout<<"Square.first    "<<Square.first <<std::endl;
    if (Sq.I == IJsq.first && Sq.J == IJsq.second)  
      return true;
    else return false;
  }
};

std::vector<double>   TCluster::ClusterRecord (std::vector<ADC> VectorIn, TReadout Readout, int  TriggerLevel,int  Width)
{ 
  std::vector<Square>  SquareV;
  std::vector<ADC>::iterator pos;
  std::vector<Square>::iterator itr;
  for(pos = VectorIn.begin() ; pos!=VectorIn.end();++pos)
    {
      std::pair<int,int> IJCluster= Readout.IJfromChannel((*pos).Channel);
      std::pair<int,int> IJSquare = Readout.GetISqJSqfromIJ(IJCluster.first, IJCluster.second);
      double Charge = (*pos).Charge;
      Square NewSquare;
      NewSquare.I = IJSquare.first;
      NewSquare.J = IJSquare.second;
      NewSquare.Charge = Charge;

      if (SquareV.size()==0)SquareV.push_back(NewSquare);
      else{
	itr = find_if(SquareV.begin(),SquareV.end(),ConfrontaSquare(IJSquare));
	if(itr == SquareV.end())
	  {
	    SquareV.push_back(NewSquare);
	  }
	else
	  {
	    (*itr).Charge += Charge;
	  }
      }	
    }
  std::vector<Square>::iterator pos2;
  int IMin=10000;
  int IMax=0;
  int JMin=10000;
  int JMax=0;
  for(pos2 = SquareV.begin() ; pos2!=SquareV.end();++pos2)
    {
      ChargeSquare.push_back((*pos2).Charge);
      if ((*pos2).Charge > TriggerLevel)//soglia trigger
	{
	  if (IMin>(*pos2).I) IMin =(*pos2).I; 
	  if (IMax<(*pos2).I) IMax =(*pos2).I; 
	  if (JMin>(*pos2).J) JMin =(*pos2).J; 
	  if (JMax<(*pos2).J) JMax =(*pos2).J; 
	  
	  2 * (*pos2).I - 1;//Posizione esagono basso sin del quadrato che triggera
	  2 * (*pos2).J - 1;
	  std::pair<int,int> PosTriggerIJ = make_pair(2 * (*pos2).I - 1, 2 * (*pos2).J - 1);
	  Readout.CenterHexagontoXY(PosTriggerIJ);
	  TriggerPixel.push_back(Readout.GetCenterHexagontoXY());//vettore con posizione quadrati che hanno triggerato
	}
    }
  if (IMin==10000 )
    {
      std::vector<double> XYRecord;
      XYRecord.push_back(0);
      XYRecord.push_back(0);
      XYRecord.push_back(0);
      XYRecord.push_back(0);
      XYRecord.push_back(0);
      return XYRecord;
    }  
  else
    {
      std::pair<int,int> MinIJ = make_pair(2 * IMin -1 -Width ,2 * JMin -1 -Width);
      Readout.CenterHexagontoXY(MinIJ);
      std::pair<double,double> XYMin = Readout.GetCenterHexagontoXY();
      
      std::pair<int,int> MaxIJ = make_pair(2 * IMax -1  +Width ,2 * JMax -1 + Width);
      Readout.CenterHexagontoXY(MaxIJ);
      std::pair<double,double> XYMax = Readout.GetCenterHexagontoXY();
      std::vector<double> XYRecord;
      //XYRecord primi 4 valori XMin,Ymin,Xmax, Ymax, poi ha numero totale pixel registrati
      XYRecord.push_back(XYMin.first);
      XYRecord.push_back(XYMin.second);
      XYRecord.push_back(XYMax.first);
      XYRecord.push_back(XYMax.second);
      
      XYRecord.push_back( (2 * (JMax-JMin)+ 2 * Width) * (2 * (IMax-IMin)+ 2 * Width)  );
      return XYRecord;
    }
}

std::pair<double,double>   TCluster:: Baricenter(std::vector<Pixel> VectorIn)
{
  double x=0.;
  double y=0.;
  double TotalCharge=0.;
  std::vector<Pixel> ::iterator pos;
  for(pos =   VectorIn.begin();pos!=  VectorIn.end();++pos)
    {
      x += ((*pos).X)* ((*pos).Charge);
      y += ((*pos).Y)* ((*pos).Charge);
      TotalCharge +=1.*((*pos).Charge);
      
    }
  y/=TotalCharge;
  x/=TotalCharge;
  std::pair<double,double> a =  std::make_pair (x,y);
  return a;
}

//restituisce theta (angolo momento II maggiore) e i valori dei Momenti secondi max e min     
double    *TCluster::PhiMajorAxis(std::vector<Pixel> VectorIn,double Xb,double Yb)
{
  double Norm = 0.;
  double *b = new double[3];
  double A=0;
  double B=0;
  double Th=0.;
  double MaxMom=0;
  std::vector<Pixel>::iterator pos;
  for (pos = VectorIn.begin(); pos != VectorIn.end() ; ++pos)
    {
      double ResX = (*pos).X - Xb;
      double ResY = (*pos).Y - Yb;
      A+=((*pos).Charge) * (ResX) * (ResY);
      B+=((*pos).Charge)*( pow(ResY,2)-pow(ResX,2));
    }
  if ((B < 0.000000001) && (B > -0.000000001)) { Th = kPI / 4.0;}
  else
    {Th = -0.5*atan (2.0*A/B);}//theta in radianti
  double CosT = cos(Th);
  double SinT = sin(Th);
  
  double MomentumTheta = 0;
  double MomentumThetaCompl = 0;
  std::vector<Pixel>::iterator pos2;
  for (pos2 = VectorIn.begin(); pos2 != VectorIn.end() ; ++pos2)
    {
      double ResX = ((*pos2).X - Xb);
      double ResY = ((*pos2).Y - Yb);
      MomentumThetaCompl += ((*pos2).Charge) *pow(CosT*ResX + SinT*ResY,2);
      MomentumTheta += ((*pos2).Charge) *pow(-SinT*ResX + CosT*ResY,2); 
      Norm+=(*pos2).Charge;
     
    }
  MomentumThetaCompl/= Norm;//X
  MomentumTheta/=Norm;//Y  
  /* voglio individuare l'asse maggiore, il piu' lungo, se X>Y e' theta, altrimenti devo invertire*/
  if (MomentumTheta>MomentumThetaCompl)
    {
      if (Th > 0.0)   Th -= 0.5*kPI;
      else Th += 0.5*kPI;
    }
  MaxMom=TMath::Max(MomentumTheta,MomentumThetaCompl);
  double MinMom = TMath::Min(MomentumTheta,MomentumThetaCompl);
   //std::cout << "Theta in gradi= "<<360*Th/(2*kPI)<<std::endl;
  b[0] = Th;
  b[1] = MaxMom;
  b[2] = MinMom;
  return b;
}

double     TCluster::ThirdMomentum(std::vector<Pixel> VectorIn,double Xb,double Yb, double T )
{ 
  double  ThirdMom = 0.; 
  double Norm = 0.;
  double CosT  = cos(T);
  double SinT  = sin(T);
  std::vector<Pixel>::iterator pos;
  for (pos = VectorIn.begin(); pos != VectorIn.end() ; ++pos)
    {
      double ResX = ((*pos).X - Xb);
      double ResY = ((*pos).Y - Yb);
      ThirdMom += ((*pos).Charge )* pow((CosT*ResX + SinT*ResY),3);
      //std::cout<<"Third mom = "<< ThirdMom/Norm<<std::endl;
      Norm+=(*pos).Charge;
    }

  ThirdMom /=Norm;
  // std::cout<<"Third mom = "<< ThirdMom<<std::endl;
  return ThirdMom;
}

std::pair<double,double> TCluster::ImpactPoint(std::vector<Pixel> VectorIn, double T,double Xb, double Yb, double ThirdM, double SecondMom)
{
  Pixelinrange = 0;
  double CosT = cos(T);
  double SinT = sin(T);
  double ImX = 0;
  double ImY = 0;
  double NormImpact = 0.;
  NStatistic = 0;
  std::vector<Pixel>::iterator pos;
  for (pos = VectorIn.begin(); pos != VectorIn.end() ; ++pos)
    {
      double ResX = ((*pos).X - Xb);
      double ResY = ((*pos).Y - Yb);
      if (SecondMom!= 0.0 && ThirdM != 0.0)
	{
	  double DistanceFromBaricenter = sqrt(ResX*ResX+ResY*ResY)/sqrt(SecondMom);
	 
	  if (DistanceFromBaricenter>SmallCircle && DistanceFromBaricenter<LargeCircle &&
	      ((CosT*ResX + SinT*ResY)/ThirdM) > 0.0 )
	    { 
	      Pixelinrange += 1;
	      ImX +=((*pos).Charge ) * ((*pos).X );
	      ImY +=((*pos).Charge ) * ((*pos).Y );
	      NormImpact +=(*pos).Charge;
	      NStatistic = 1;
	    }
	}
      else 
	{
	  NormImpact=1.;
	  ImX=Xb;
	  ImY=Yb;
	}
    }
  
  if (NormImpact==0)  
    {
     return  std::make_pair (Xb,Yb);
    } 
  ImX/=NormImpact;
  ImY/=NormImpact;
  //std::cout<<"ImpactX = "<< ImpactX<<" ImpactY ="<< ImpactY<<std::endl;
  std::pair<double,double> b = std::make_pair (ImX,ImY);

  return b;
}

// pesa diversamente la carica dei pixel in funzione della distanza da un punto, che sara' il punto d'impatto ricostruito

std::vector<Pixel> TCluster::PartialWeight (std::vector<Pixel> VectorIn, double ImpX, double ImpY)
{
  
  std::vector<Pixel>::iterator pos;
  double Cut = 10;
 
  for(pos =  VectorIn.begin();pos!= VectorIn.end();++pos)
    {
      double Distance=sqrt(pow((*pos).X - ImpX,2.) + pow((*pos).Y - ImpY,2.0));//sqrt(SecondMom);  TOLTA NORMALIZZAZIONE
  
      // Tentativo tagli (26/5/05)
      if (VectorIn.size() > 200) Cut = 6.;
      else if(VectorIn.size() > 100) Cut = 5.;
      
      if (Distance< Cut*100)
      {
	
	(*pos).Charge *= exp(-Distance/WeightScale);
      }
      else
	{  
	  (*pos).Charge = 0;
	}
    }
  
  return VectorIn;
  
}
void            TCluster::Analysis( std::vector <Pixel> V)
{ 
  // Calcolo Baricentro
  std::pair <double,double> Bar = Baricenter(V);  //<<<<<===================  BARICENTER
  X = Bar.first;
  Y = Bar.second;
  double *ThMMmM;
  //Calcolo Momenti secondi e Angolo Theta
  ThMMmM= PhiMajorAxis(V , X , Y);               //<<<<<===================    M2
  Theta =  ThMMmM[0];
  MaxM = ThMMmM[1];
  MinM =  ThMMmM[2];
  delete[] ThMMmM;
  //Momento Terzo
  Third = ThirdMomentum( V , X , Y , Theta);     //<<<<<===================    M3
  std::pair<double,double> Impact;
  //Controllo che Small Circle non sia troppo piccolo, che ci siano almeno 4 pixel
   while (Pixelinrange <= 4 && SmallCircle > 0.2)
    {    
      Impact = ImpactPoint( V ,Theta , X , Y , Third , MaxM);
      SmallCircle -= 0.2;
      // LargeCircle -=0.2;//REAL     
    }
  ImpactX = Impact.first;
  ImpactY = Impact.second;
  //ripeso tutto il vettore di pixel in funzione della distanza dal punto d'impatto
  std::vector<Pixel> Wvector;
   
   int  Mode = 2;
  if (Theta>=0. && Third<=0.)//UP
    {
      Mode=1;
    }
  if (Theta>0. && Third>0.)//down
    {
      Mode=0;
    }
  if (Theta<=0. && Third<=0.)//down
    {
      Mode=0;
    }
  if (Theta<0. && Third>0.)//UP
    {
      Mode=1;
    }
  switch(Mode)
    {
    case 1:
      if (Theta<0)
	{Theta+=kPI;}
      break;
      /////////
    case 0:
       if (Theta<=0)
	{Theta=2*kPI + Theta;}
      else {Theta+=kPI;}
       break;
    }
  double TypeAnCut =10;//20;
  if (Pitch == 80) TypeAnCut =8;
  
  if ( (V.size()<TypeAnCut && MaxM/MinM<3) || (ImpactX == X && ImpactY== Y ))//(V.size()<TypeAnCut || (ImpactX == X) || MaxM/MinM<1.5)
    {
      //std::cout<<" FIRST STEP ANALYSIS "<<std::endl;
      ImpactX = X;
      ImpactY = Y;
      Theta2 = Theta;
    }
  else 
    {
      // std::cout<<" SECOND "<<std::endl;
      Wvector = PartialWeight (V , ImpactX , ImpactY);
      
      //calcolo nuovo baricentro e nuovo momento secondo, con Theta2 che mi da' l'angolo ricostruito al secondo passo
      std::pair <double,double> Barb = Baricenter (Wvector);
      X2 = Barb.first ;
      Y2 = Barb.second ;
      double *ThMM2;
      ThMM2 = PhiMajorAxis(Wvector , X2, Y2 );
      Theta2 =  ThMM2[0];
      
      delete[] ThMM2;
      //scelgo l'orientamento di Theta primo passo
      
      //scelgo l'orientamento di Theta2, ricalcolando il momento terzo, rispetto al baricentro X2, Y2
      double Thirdbis = ThirdMomentum( V , X2 , Y2 , Theta2);
      int  Modebis = 2;
      if (Theta2>=0. && Thirdbis<=0.)//down
	{
	  Modebis=0;
	}
      if (Theta2>0. && Thirdbis>0.)//up
	{
	  Modebis=1;
	}
      if (Theta2<=0. && Thirdbis<=0.)//up
	{
	  Modebis=1;
	}
      if (Theta2<0. && Thirdbis>0.)//down
	{
	  Modebis=0;
	}
      
      switch(Modebis)
	{
	case 1:
	  if (Theta2<0)
	{Theta2+=kPI;}
	  break;
	case 0:
	  if (Theta2<=0)
	    {Theta2=2*kPI + Theta2;}
	  else {Theta2+=kPI;}
	  break;
	}
    }
}
//Angle in radianti
std::pair<double,double>    TCluster::Froc(double Xzero,double Yzero, double Angle)
  
{
 
  double Xf=400*cos(Angle)+Xzero;
  double Yf=400*sin(Angle)+Yzero;
  std::pair<double,double> Ex = make_pair(Xf,Yf);
  
  return Ex;

}


//Angle in radianti
//Ruota attorno al punto Origin tutto il cluster
std::vector<Pixel>          TCluster:: RotateVector(std::vector<Pixel>  VectorIn,double cosROT,double sinROT , double XOrigin, double YOrigin)
      
{
  std::vector<Pixel> VectorOut;
  std::vector<Pixel>::iterator pos;
  for(pos = VectorIn.begin();pos!=VectorIn.end();++pos)
    { 
      //
      if ((*pos).Charge>4)
      //
	{
	  Pixel Rot;
	  // Rot.X=((*pos).X - XOrigin)*cosROT + ((*pos).Y - YOrigin)*sinROT;
	  //Rot.Y= -((*pos).X - XOrigin)*cosROT + ((*pos).Y - YOrigin)*sinROT;
	  Rot.X=((*pos).X - XOrigin)*cosROT - ((*pos).Y - YOrigin)*sinROT;
	  Rot.Y= -((*pos).X - XOrigin)*sinROT + ((*pos).Y - YOrigin)*cosROT;
	  Rot.Charge = (*pos).Charge;
	  VectorOut.push_back(Rot);
	}
    }
  return VectorOut;
}

void     TCluster::ClusterProfile(std::vector<Pixel> VectorIn, TH1D *XProfile, TH1D *YProfile)

{
   std::vector<Pixel>::iterator pos;
   for (pos = VectorIn.begin();pos!=VectorIn.end();++pos)
     {
       XProfile->Fill((*pos).X,(*pos).Charge/10);
       YProfile->Fill((*pos).Y,(*pos).Charge/10);
  
     }
  
}
  

  
