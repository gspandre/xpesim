#include "TReadout.h"

TReadout::TReadout(TDimension *Dimension)
{  
  std::pair<int,int>  Pi=Dimension->GetPixel();
  Setdimension(1, Pi.first,1, Pi.second);
  SetPitch(Dimension->GetPitch());
}
void TReadout::SetPitch(double p)
{
  Pitch=p;
}
double TReadout::Modulo(double x, double y)
{
  double iy = floor(x/y);
  return (double) x-iy*y;
}

void TReadout::Setdimension(int xmin ,int xmax,int ymin,int ymax)
  //(1->160,1->138)cluster verticali//
{
  JMin=xmin;
  JMax=xmax;//160
  IMin=ymin;
  IMax=ymax;//138 
}

//I riga contata dal basso, dipende dalla posizione Y,J colonna contata da sinistra, dipende dalla X, X e Y in cm 
//con origine degli assi al centro del rivelatore
//origine riferimento meta lato in basso dell'esagono, riga e colonna 0

std::pair<int,int>  TReadout::PositiontoMatrix(double Xa, double Ya)
{//AGGIUNGO 1 A IMAX E JMAX
  Y=Ya*10000.+Pitch*(IMax+2)/2.+25;//aggiungo 25 per centrare// 
  L=Pitch/sqrt(3.);
  X=(Xa*10000.+(JMax+3)*Pitch*sqrt(3.)/4.)-L+7; //7 per centrare 
  //std::cout<<"Xmod "<<X<<" Y mod "<<Y<<std::endl;
  PassoX=Pitch*sqrt(3.)/2;
  PassoY=Pitch/2;
  Xx=Modulo(X,PassoX);
  int  j = int (X/PassoX);
  
  if (Xx <= L/2 && j%2 == 0)// per j pari colonna parte da pitch/2
    { // std::cout<<"caso1"<<std::endl;
      I=(Y-Pitch/2)/Pitch;
      J=j;
    }
  else
    {
      if (Xx<=L/2 && j%2 == 1)
	{// std::cout<<"caso2"<<std::endl;
	  I=Y/Pitch;
	  J=j;
	}   
      else
	{
	  if (Xx >= L && j%2== 0)
	    {// std::cout<<"caso3"<<std::endl;
	      J=1+(X/PassoX);
	      I=Y/Pitch;//(Y-Pitch/2)/Pitch;
	    }
	  else
	    {
	      if (Xx >= L && j%2 ==1)
		{// std::cout<<"caso4"<<std::endl;
		  J=1+(X/PassoX);//(Y/Pitch);
		  I=(Y-Pitch/2)/Pitch;//(X/PassoX+1);
		}
	      else
		{//
		  
		  double ysub=Modulo(Y-Pitch/2,Pitch);
		  // std::cout<<ysub<<std::endl;
		  // std::cout<<Pitch-sqrt(3)*(Xx-L/2)<<std::endl;
		  //int Yy=int(ysub);
		  //if (Yy%2== 0 && j== 0 && Yy<Pitch-sqrt(3)*Xx)
		  if (ysub>=Pitch/2  && j%2== 0 && ysub<Pitch-sqrt(3.)*(Xx-L/2))  
		    {// std::cout<<"caso5"<<std::endl;
		      J=(X/PassoX);
		      I=((Y-Pitch/2)/Pitch);
		    }
		  else
		    {
		      if (ysub>=Pitch/2  && j%2== 0 && ysub>=Pitch-sqrt(3.)*(Xx-L/2))
			//  if (Yy%2== 0 && j== 0 && Yy>=Pitch-sqrt(3)*Xx)
			{ // std::cout<<"caso6"<<std::endl;
			  J=1+X/PassoX;
			  I=Y/Pitch;
			}
		      else  
			{
			 if (ysub<=Pitch/2  && j%2== 0 && ysub<=sqrt(3.)*(Xx-L/2))
			   // if (Yy%2==1 && j== 0 && Yy<sqrt(3)*Xx-Pitch/2)
			   {// std::cout<<"caso7"<<std::endl;
			     J=1+X/PassoX;
			     I=Y/Pitch;
			   }
			 else
			   { 
			     if (ysub<=Pitch/2 && j%2== 0 && ysub>sqrt(3.)*(Xx-L/2)) 
			       // if (Yy%2== 1 && j== 0 && Yy>=sqrt(3)*Xx-Pitch/2) 
			       {// std::cout<<"caso8"<<std::endl; 
				 J=X/PassoX;
				 I=(Y-Pitch/2)/Pitch;
			       }
			     else
			       {
				 double ys=Modulo(Y,Pitch);
				 if (ys>=Pitch/2 && j%2== 1 && ys<=Pitch-sqrt(3.)*(Xx-L))
				   // if (Yy%2== 1 && j== 1 && Yy<Pitch-sqrt(3.)*Xx)
				   {// std::cout<<"caso9"<<std::endl;
				     J=X/PassoX;
				     I=Y/Pitch;
				   }
				 else 
				   {
				     if (ys>=Pitch/2 && j%2== 1 && ys>Pitch-sqrt(3.)*(Xx-L))
				       //if (Yy%2== 1 && j== 1 && Yy>=Pitch-sqrt(3.)*Xx)
				       {// std::cout<<"caso10"<<std::endl;
					 J=1+X/PassoX;
					 I=(Y-Pitch/2)/Pitch;
				       }
				     else
				       {
					 if (ys<=sqrt(3.)*(Xx-L))
					   {// std::cout<<"caso11"<<std::endl;
					     J=1+X/PassoX;
					     I=(Y-Pitch/2)/Pitch;
					   }
					 else
					   {// std::cout<<"caso12"<<std::endl;
					     J=X/PassoX;
					     I=Y/Pitch;
					   }
				       }
				   }
			       }
			   }
			}
		    }
		}
	    }
	}
    }
  
  if ((int)I>=IMin && 
      (int)I<IMax &&
      (int)J>=JMin &&
      (int)J<JMax )
    { 
      std::pair<int,int> a = make_pair((int) I,(int)J);
      return a;
    }
  else 
    {
      std::cout <<"ATTENZIONE coordinate fuori dal detector!!!  "<< " X= "<< Xa<< " Y= "<< Ya <<" I ="<<I<<" J = "<<J<<std::endl;
      std::pair<int,int> a = make_pair(-9999,-9999); // aggiunto il 13/03/08
      return  a;
    }
}

void  TReadout::CenterHexagontoXY (std::pair<int,int> ij )//i,j
{
   XHexagonCenter=(ij.second)*Pitch*sqrt(3.)/2;
   //XHexagonCenter=(ij.second)*Pitch*sqrt(3.)/2 -3507;//MODIFICATO
  if (ij.second%2 ==0)
    {
      YHexagonCenter=(ij.first+1)*Pitch;
      //YHexagonCenter=(ij.first+1)*Pitch-3525;//MODIFICATO
    }
  else
    {
      YHexagonCenter= (ij.first+1/2.)*Pitch;
      //YHexagonCenter= (ij.first+1/2.)*Pitch-3525;//MODIFICATO
    }
}

std::pair<double,double>  TReadout::GetCenterHexagontoXY()
{ 
  std::pair<double,double> XYHexagonCenter = std::make_pair( XHexagonCenter,YHexagonCenter);
  return XYHexagonCenter;
}

//8 clusters messi in verticale con 20 colonne e 138 righe ciascuno
//da canale 0 a canale 137,non generalizzato

void TReadout::IJtoClusterChannel (int rig,int col)
{
  Cluster=(col-1)/20;
  if (col%2==0)
    {
      Channel=IMax*((col-1)%20)+IMax-rig;//-1per partire dal canale 0
    }
  else
    {
      Channel=rig+IMax*((col-1)%20)-1;
    }
}

std::pair<int,int>  TReadout::GetClusterChannelfromIJ()
{
  std::pair<int,int> ClusterChannelfromIJ = make_pair( Cluster,Channel );
  return ClusterChannelfromIJ;
}
 

void TReadout::ClusterChanneltoIJ (int clu,int chan)
{
  J=1+clu*20.+chan/IMax;
  Jint=(int)J;
  if (Jint%2==0)
    {
      Iint=IMax-chan%IMax;
     
    }
  else
    {
      Iint=chan%IMax+1;
      
    }
}

std::pair<int,int> TReadout::GetIJfromClusterChannel()
{
  std::pair<int,int> IJfromClusterChannel = make_pair( Iint, Jint );
  return  IJfromClusterChannel;
}

void  TReadout::ChannelfromIJ(int i, int j)
{
  //Per chip 22k
  /*if (j%2==1)
    {
      GChannel=i-1+(j-1)*IMax;//138;
    }
  else
    {
      GChannel=(j-1)*IMax+IMax-i;//138+138-i;
    } 
 */

  //Per chip 105k
  GChannel=i+(j)*IMax;//i-1+(j-1)*IMax;
}
std::pair<int,int>   TReadout::IJfromChannel(int ch)
{
  // Per chip 22k
  /*
    j=(ch)/IMax+1;//138+1;tolto +1 a ch numeratore
  if (j%2==0)
    {
      I=IMax-(ch)%IMax;//138.-ch%138;
    }
  else
    { 
      I=(ch)%IMax+1.;
    }
  J=j*1.;
  std::pair<int,int> IJfromC = make_pair((int) I, j); 
  return IJfromC;
  */
 
  //Per chip 105k
  j=(ch)/IMax;//+1;
    I=(ch)%IMax;//+1.;
  J=j*1.;
  std::pair<int,int> IJfromC = make_pair((int) I, j); 
  return IJfromC;
}
std::pair<int,int>   TReadout::GetISqJSqfromIJ(int i, int j)
{
  int ISq = (int) (i+1)/2;
  int JSq = (int) (j+1)/2;
  std::pair<int,int> ISqJSq = make_pair(ISq,JSq);
  return ISqJSq;
}

void TReadout::myDraw3DPixel(int chan, double Size)
{ 
  PixMapCanvas->cd();
  std::pair<int,int> ij =  IJfromChannel(chan);
  CenterHexagontoXY (ij);
  
  if (Size > 1.0) Size = 1.0;
  if (Size > 0.0)
    {
      THexagon3d *pixel = new THexagon3d(GetCenterHexagontoXY().first,GetCenterHexagontoXY().second,0, Size * Pitch);
      pixel->Draw(2,2);
      PixMapCanvas->Update();
    }
  
}

void TReadout::myDraw2DPixel(int chan, double Size)
{ 
  PixMapCanvas->cd();
  std::pair<int,int> ij =  IJfromChannel(chan);
  CenterHexagontoXY (ij);
  
  if (Size > 1.0) Size = 1.0;
  if (Size > 0.0)
    {
      THexagon2d *pixel = new THexagon2d(GetCenterHexagontoXY().first,GetCenterHexagontoXY().second, Size * Pitch/2);//Size
      pixel->Draw(1,1);
    }
  
}

/*
void TReadout::myDraw3D(int Color, double Size, int LineWidth)
{
  PixMapCanvas = new TCanvas("PixMapCanvas","PixMapCanvas",800,800);
  PixMapCanvas->SetFillColor(10);
  TView *view = new  TView(1); 
 /// Gloria 20/07 commented because class TView in root5.16 is an abstract class and cannot be instatiated
  double SF    = 1.;
  double MaxX = 10000.;
  double MaxY = 10000.;
  double MinX = 0.0;
  double MinY = 0.0;
  view->SetRange(SF*MinX, SF*MinY, 0, SF*MaxX, SF*MaxY,1);
  view->Top();
  PixMapCanvas->Update();
}
*/

void TReadout::myDraw2D(int Color, double Size, double LineWidth)
{
  PixMapCanvas = new TCanvas("PixMapCanvas","PixMapCanvas",1200,800);
  PixMapCanvas->SetFillColor(10);
  PixMapCanvas->Draw();
  double SF    = 1.;
  double MaxX = 10000.;
  double MaxY = 10000.;
  double MinX = 0.0;
  double MinY = 0.0;
  PixMapCanvas->Range(-2000,-2000,1.25*SF*MaxX,1.25*SF*MaxY);
  TGaxis *xaxis = new TGaxis(SF*MinX,SF*MinY,SF*MaxX,SF*MinY,SF*MinX,SF*MaxX);
  TGaxis *yaxis = new TGaxis(SF*MinX,SF*MinY,SF*MinX,SF*MaxY,SF*MinY,SF*MaxY);

  xaxis->Draw();
  yaxis->Draw();
  PixMapCanvas->Update();
}

std::pair<double,double>  TReadout:: XYorigin(std::pair<double,double> xymod)

{
  double Xorigin = (( xymod.first+Pitch/sqrt(3.))-(JMax+1)*Pitch*sqrt(3.)/4.)/10000;
  double Yorigin = (xymod.second-Pitch*(IMax+1)/2.)/10000.;
  std::pair<double,double> XYout = make_pair( Xorigin, Yorigin);
  return XYout;

}
