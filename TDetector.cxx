#include "TDetector.h"

TDetector::TDetector(TGasMixture *MIXTURE, TRandom *RND, Bool_t VERBOSE,TDimension *Dimension)
{
  rng = RND;

  Gem_Radius  =   Dimension->GetGem_Radius();
  Gem_Pitch   =   Dimension->GetGem_Pitch();    //mu
  Pitch       =   Dimension->GetPitch();    //mu
  Pi          =   Dimension->GetPixel();
  PixelX      =   Pi.first;   //J  
  PixelY      =   Pi.second;   //I
  
  NoiseRMS =  2.0;
  Pedestal = 1400;
  ElectronsPerADC = Dimension->GetElectronsADC();///60

  Gem    =    new TGem(rng,Dimension);

}



void TDetector::TransferToGEM(TTrack Track)
{
  Track.Drift();//Drift nuovo,in ttrack ritorna una coppia, non modifica il vettore primaryionozation
}

int TDetector::Digitization(int Charge) /// ritorna in conteggi di ADC il segnale + piedistallo + noise 
{
  int DigiSignal = (Int_t) (Charge/ElectronsPerADC);
  int DigiNoise  = (Int_t) rng->Gaus(0.0,NoiseRMS);
  return TMath::Max(0,Pedestal - DigiSignal + DigiNoise);
}


std::vector<ADC>   TDetector::mySampling(std::vector<std::pair<double,double> >  XYSec, TReadout Readout)
{
  int nelectr=Gem->GetPac();
  std::vector<int> chvector;
  std::list<int> chlist;
  std::pair<int,int> a; 
  std::vector<std::pair<double,double> >::iterator pos;
  for (pos=XYSec.begin();pos !=XYSec.end(); ++pos )
    {
      //a  = Readout.PositiontoMatrix((*pos).first,(*pos).second);
      a  = Readout.Position2Xpol((*pos).first,(*pos).second);
      //std::cout<<"x = "<<(*pos).first<<" Y= "<<(*pos).second<<std::endl;
      //std::cout<<"I = "<<a.first<<" J= "<<a.second<<std::endl;
      Readout.ChannelfromIJ(a.first,a.second);
      int c=Readout.GetChannel();
      chvector.push_back(c);
      chlist.push_back(c);
      //if ((c==52050) || (c==52350) || (c==52650)){
      //std::cout<<"c= " << c <<" x= "<<(*pos).first<<" y= "<<(*pos).second <<" I= "<<a.first<<" J= "<<a.second<<std::endl;
      //}
    }
  
  chlist.sort();
  chlist.unique();
  std::vector<ADC> Output;
  ADC Signal; 
  std::list<int>::iterator poss;
  int ncount=0;
  for (poss = chlist.begin() ;poss !=chlist.end(); ++poss )
    {
      Signal.Channel =(*poss);
      int counter=count(chvector.begin(),chvector.end(),Signal.Channel);
      Signal.Charge=nelectr*counter;
      Signal.DigiCharge=Digitization(Signal.Charge);
      Output.push_back(Signal);
      ncount+=Signal.Charge;
     }
  //int Os=Output.size();
  return Output;
}


int TDetector::CountChannels(std::vector<int> chvect, int ch)
{
  std::vector<int>::iterator pos;
  int count=0;
  sort(chvect.begin(),chvect.end());
  for (pos = chvect.begin() ;pos !=chvect.end(); ++pos )
    {
      if((*pos)==ch) count++;
      else if((*pos)>ch) return count;
    }
  return count;
}


void TDetector::myDrawSignal(std::vector<ADC> Digi,TReadout Readout)
{
  Readout.myDraw2D(4,1,.5);
  int N = Digi.size();
  std::list<double> x;
  std::list<double> y;
  std::list<int> Sizes; 
  int effClusize = 0;
  double sigmax = 0.;
  double signal = 0.; 
  int isigmax,jsigmax,chann;
  for(int i = 0; i< N; i++)
    {
      // if(Digi[i].DigiCharge<1400)
   	{
	  //double maxCharge = Pedestal;
	  Readout.CenterHexagontoXY(Readout.IJfromChannel(Digi[i].Channel));
	  x.push_back((Readout.GetCenterHexagontoXY()).first);
	  y.push_back((Readout.GetCenterHexagontoXY()).second);
	  double Sig = (Pedestal-Digi[i].DigiCharge);
	  signal = 1400-Digi[i].DigiCharge;
	  if(signal) {
	    effClusize++;
	    if (signal>sigmax) {
	      sigmax = signal;
	      isigmax = Readout.IJfromChannel(Digi[i].Channel).first;
	      jsigmax = Readout.IJfromChannel(Digi[i].Channel).second;
	      chann = Digi[i].Channel;
	    }
	  }
	  Sizes.push_back((int)Sig);

	  //std::cout<<"Clusize  " << N << " i = " << Readout.IJfromChannel(Digi[i].Channel).first << " j = " << Readout.IJfromChannel(Digi[i].Channel).second
	  //	   <<" charge <- " <<Digi[i].DigiCharge << " X = "<<(Readout.GetCenterHexagontoXY()).first << " Y =" << (Readout.GetCenterHexagontoXY()).second
	  //	   << std::endl;
	}
    }
  std::cout << " ===>> Effective Clusize (Signal > Pedestal (= 1400 ADC counts)): " << effClusize << std::endl;
  //std::cout << " coord of highest pixel ==> I: " << isigmax << "  J: " << jsigmax << " -- channel: " << chann << std::endl;
  x.sort();
  y.sort();
  Sizes.sort();
  int MaxSize = Sizes.back();
  for(int i = 0; i< N; i++)
    {
      double Size = (Pedestal-Digi[i].DigiCharge);
      
      Size/=MaxSize;
      if(Size>0)
	{
	  Readout.myDraw2DPixel(Digi[i].Channel,Size);
	}
    }

  double MaxX=x.back();
  double MaxY=y.back();
  double MinX=x.front();
  double MinY=y.front();
  
  TGaxis *Xax = new TGaxis(MinX-100,MinY-140,MaxX+170,MinY-140,MinX-100-7664.11,MaxX+170-7664.11);
  TGaxis *Yax = new TGaxis(MinX-100,MinY-140,MinX-100,MaxY+170,MinY-140-7575,MaxY+170-7575);
  
  Xax->SetTitle("#mum ");
  Yax->SetTitle("#mum ");
  Xax->SetLabelSize(0.027);
  Xax->SetTitleOffset(0.7);
  Yax->SetLabelSize(0.027);
  Yax->SetTextAlign(31);
  Yax->SetTitleOffset(1.1);
  Xax->Draw();
  Yax->Draw();
  
  gPad->Range( MinX-500,MinY-500,MaxX+500,MaxY+500);
  gPad->Update();
  
}
