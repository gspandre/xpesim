#include "TGem.h"

TGem::TGem(TRandom *RNG,TDimension *Dimension)
{
  rng=RNG;
  Pitch=Dimension->GetGem_Pitch();
  Gain=Dimension->GetGain();
  Pac=Dimension->GetPac();
  Hole=Dimension->GetHole();
  XStep = Pitch;
  YStep = Pitch*sqrt(3.0);
  TraslX =floor(10000/XStep);
  TraslY =floor(10000/YStep);
}

int TGem::GetSecondaryElectrons()
{
  return (int) (rng->Exp(Gain));
}

std::pair<double,double> TGem::Sampling(double x, double y)
{  
  float  px  = Pitch;
  double py  = px*sqrt(3.)/2.;
  int nx,ny;
  double cx, cy;
  std::pair<double,double> XY;
  if(x==0) nx = 0;
  else nx = (int) (2.*x/px + 0.5*x/fabs(x));
  if(y==0) ny = 0;
  else ny = (int) (y/py + 0.5*y/fabs(y));
  cx = nx*px/2.;
  cy = ny*py;
  
  if((nx+ny)%2==0) 
    {
      XY = std::make_pair(cx,cy);
    }
  else 
    {
      double rx = (x-cx);
      double ry = (y-cy);
      double MinimumDist=2*px/sqrt(3.);
      double Dist;
      double Xm;
      double Ym;
      for(int i = -1; i<=1 ; i++)
	{
	  for(int j = -(1-abs(i)); j<=(1-abs(i)) ;j+=2)
	    {
	      double xm = rx-i*px/2.;
	      double ym = ry-j*py;
	      double Dist = sqrt(pow(xm,2.)+pow(ym,2.));
	      if(Dist < MinimumDist)
		{
		  Xm          = i*px/2.;
		  Ym          = j*py;
		  MinimumDist = Dist;
		} 
	    }
	}
      XY = std::make_pair(cx + Xm,cy + Ym);
      if(MinimumDist>px/sqrt(3.)) std::cout<<"puppa"<<std::endl;
    }
  
  return XY;
} 
/*creazione elettroni secondari con distribuzione gaussiana campionati espressi
  in cm (posizione xy)*/

std::vector<std::pair<double,double> >  TGem::DiffusionofSecondaryElectrons(std::vector<std::pair<double,double> > VectorIn ) 
{ 
  AvalangeSigma=1./10000.* Pitch/2.5 ;// da /3 a /2 per i dati reali(25/11/05), nessun dato simulato con questa diffusione

  std::vector<std::pair<double,double> >   XYSecondary;
  double x;
  double y;
  int NumberPac;
  std::vector<std::pair<double,double> >::iterator pos;
  for (pos=VectorIn.begin();pos !=VectorIn.end(); ++pos )
    {     
      x=(*pos).first;
      y=(*pos).second;
      //
      std::pair<double,double> xyGem =  GemSampling(x,y);
      x= xyGem.first;
      y= xyGem.second;
	//
	
	NumberPac = GetSecondaryElectrons()/Pac;//1+( (double)GetSecondaryElectrons()/((double )Pac));
      //std::cout<<"number secondary electron"<<NumberPac*Pac<<std::endl; 
      for (int i=0; i<NumberPac; i++)
	{
	  double xp=x + rng->Gaus(0,AvalangeSigma);
	  double yp=y + rng->Gaus(0,AvalangeSigma);
	  std::pair<double,double> xy = std::make_pair( xp, yp);	  
	  XYSecondary.push_back(xy);
	}
    }
  //  std::cout<<"Dimensio of Avalance"<<XYSecondary.size()<<std::endl;
  return XYSecondary;
}

std::pair<double,double> TGem::GemSampling(double xin, double yin)
{  
  xin = xin*10000 +  XStep*TraslX;
  yin = yin*10000 +  YStep*TraslY;
  double NumberStepX = floor(xin/XStep) * XStep;
  double NumberStepY = floor(yin/YStep) * YStep;
  double XRed = xin - NumberStepX;
  double YRed = yin - NumberStepY;
  std::set<double> c;  // contenitore ordinato in ordine crescente di oggetti tutti diversi (il primo elemento- c.begin()- ha il valore piu' piccolo)
  double DistanceI = sqrt(XRed*XRed + YRed*YRed);
  double DistanceII = sqrt((XRed-XStep)*(XRed-XStep) + YRed*YRed);
  double DistanceIII= sqrt((XRed-XStep/2.)*(XRed-XStep/2.) + (YRed-YStep/2.)*(YRed-YStep/2.));
  double DistanceIV = sqrt(XRed*XRed + (YRed - YStep)*(YRed - YStep));
  double DistanceV = sqrt((XRed-XStep)*(XRed-XStep) + (YRed - YStep)*(YRed - YStep));
  c.insert(DistanceI);
  c.insert(DistanceII);
  c.insert(DistanceIII);
  c.insert(DistanceIV);
  c.insert(DistanceV);
  std::set<double>::iterator pos;
  pos = c.begin();
  std::pair<double,double> exit;
  if (*pos == DistanceI) 
    {
      exit = std::make_pair(NumberStepX,NumberStepY);
    }
  else if (*pos == DistanceII)
     {
      exit = std::make_pair(NumberStepX + XStep,NumberStepY);
     }
   else if (*pos == DistanceIII)
     {
      exit = std::make_pair(NumberStepX + XStep/2,NumberStepY + YStep/2);
     }
   else if (*pos == DistanceIV)
     {
      exit = std::make_pair(NumberStepX ,NumberStepY + YStep);
     }
   else if (*pos == DistanceV)
     {
      exit = std::make_pair(NumberStepX + XStep,NumberStepY + YStep);
     }
  //la valanga si produce ai bordi del foro della gem,non nel centro
  double PhiRandom = rng->Uniform(0.,2.*TMath::Pi());
  //std::cout<<" PhiRandom = "<< PhiRandom<<std::endl;
  exit.first =(exit.first- XStep*TraslX+Hole*cos(PhiRandom))/10000;
  exit.second =(exit.second - YStep*TraslY+Hole*sin(PhiRandom))/10000;
  return exit;
}
