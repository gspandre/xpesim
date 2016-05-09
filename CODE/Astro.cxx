#include <iostream>
#include <vector>
#include "TGraph.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "TExperiment.h"
#include "FastSimulator.h"
#include "TMDP.h"


 void RunFastSimulator(int Nevents,double Energy, double &Mu,double &MuE, double &AngleErr)
{
  std::cout<<"... processing "<<Nevents<<" events of energy  "<<Energy<<" keV"<<std::endl;
  FastSimulator FS(false);
  FS.SetOption(5,1,1);
  FS.SetEnergy(Energy);
  FS.SetPolarizationFraction(0.1);
  FS.SetPolarizationAngle(0.0);
  FS.ProcessEvents(Nevents);
  FS.GetResults(Mu,MuE,AngleErr);
  std::cout<<"... result: Mu "<<Mu<<" +- "<<MuE<<", deltaPhi= "<<AngleErr<<std::endl;
 
}

void PlotSpectrum()
{    
  TRandom3 *rnd = new TRandom3();
  double EMIN=2.0;
  double EMAX=10.0;
  int    MixId = 5;
  double Pressure = 1.0;
  double Thickness = 1.0;
  
  TExperiment *Experiment = new TExperiment(rnd,"pippo");

  Experiment->SetMixID(MixId);
  Experiment->SetPressure(Pressure);
  Experiment->SetThickness(Thickness);
  Experiment->Generalsetup();  
  TMDP *mdp = new TMDP();
  
  TF1* blazar = mdp->Spectrum2(2.47,9.9e-3);
  mdp->SetMirror("XEUSNEW.txt");      
  TGraph *myMirror = mdp->MirrorGraph(1.0);
  TGraph *myWindow = mdp->GetWindowEfficiency();
  TGraph *myGasEff = Experiment->GetEfficiencyGraph();
  
  myGasEff->SetLineColor(2);
  //////////////////////////////////////////////////

  TH1D *SourceSpectrum = mdp->Spectrum(blazar,EMIN,EMAX);
  SourceSpectrum->Scale(5e4);
    
  TH1D *SourceSpectrumAeff = mdp->Spectrum(blazar,myMirror,EMIN,EMAX);
  TH1D *SourceSpectrumAeffWind = mdp->Spectrum(blazar,myMirror,myWindow,EMIN,EMAX);
  TH1D *SourceSpectrumAeffWindGas = mdp->Spectrum(blazar,myMirror,myWindow,myGasEff,EMIN,EMAX);
  //SourceSpectrum->SetLineColor();
  SourceSpectrumAeff->SetLineColor(2);
  SourceSpectrumAeffWind->SetLineColor(3);
  SourceSpectrumAeffWindGas->SetLineColor(4);
  //////////////////////////////////////////////////  
  double TOBS=1e5;
  double R1 = mdp->Rate(SourceSpectrum,EMIN,EMAX);
  double R2 = mdp->Rate(SourceSpectrumAeff,EMIN,EMAX);
  double R3 = mdp->Rate(SourceSpectrumAeffWind,EMIN,EMAX);
  double R4 = mdp->Rate(SourceSpectrumAeffWindGas,EMIN,EMAX);
  std::cout<<"Source Rate *50000           : "<<R1<<" Tobs : "<<TOBS<<" Nevents: "<<R1*TOBS<<std::endl;
  std::cout<<"Source Rate Mirror           : "<<R2<<" Tobs : "<<TOBS<<" Nevents: "<<R2*TOBS<<std::endl;
  std::cout<<"Source Rate Mirror Window    : "<<R3<<" Tobs : "<<TOBS<<" Nevents: "<<R3*TOBS<<std::endl;
  std::cout<<"Source Rate Mirror Window Gas: "<<R4<<" Tobs : "<<TOBS<<" Nevents: "<<R4*TOBS<<std::endl;
  
  TCanvas *ac = new TCanvas("ac","ac");

  ac->Divide(2,2);
  ac->cd(1);
  myWindow->Draw("alp");
  myGasEff->Draw("lp");
  ac->cd(2);
  gPad->SetLogx();
  gPad->SetLogy();
  myMirror->Draw("alp");

  ac->cd(3);
  gPad->SetLogx();
  gPad->SetLogy();
  blazar->Draw();

  ac->cd(4);
  gPad->SetLogx();
  gPad->SetLogy();
  //  obsSpectrum->SetLineStyle(2);
  SourceSpectrum->SetLineStyle(2);
  SourceSpectrum->Draw();
  SourceSpectrumAeff->Draw("same");
  SourceSpectrumAeffWind->Draw("same");
  SourceSpectrumAeffWindGas->Draw("same");
  //////////////////////////////////////////////////

  static const int Nenergies = 8;
  static const int Ntimes    = 10;
  
  double modulation_factors[Nenergies][Ntimes];
  double modulation_factors_err[Nenergies][Ntimes];
  double reconstruction_angle_err[Nenergies][Ntimes];

  int NumberOfEvents[Nenergies][Ntimes];
  double *TOB = new double [Ntimes];

  for (int t = 0; t < Ntimes; t++)
    {
      TOB[t]=pow(10,5+0.4*t);
      
      for (int i = 0; i<Nenergies;i++)
	{
	  NumberOfEvents[i][t] =(int) (mdp->Rate(SourceSpectrumAeffWindGas,2.+i,3.+i) * TOB[t]);
	  std::cout<<" ****************************** The Number of events between "<<2.+i<<" and "<<3.+i<<" keV in "<<TOB[t]<<" seconds, is "<< NumberOfEvents[i][t]<<std::endl;
	  if(NumberOfEvents[i][t] < 10000)
	    {
	      
	      modulation_factors[i][t]=-1.0;
	      modulation_factors_err[i][t]=0.0;
	      reconstruction_angle_err[i][t]=0.0;
	    }
	  else
	    {
	      RunFastSimulator(NumberOfEvents[i][t],2+i,
			       modulation_factors[i][t],
			       modulation_factors_err[i][t],
			       reconstruction_angle_err[i][t]);
	    }
	  //	  modulation_factors[i][t]=1.0;
	  //modulation_factors_err[i][t]=0.1;
	  //reconstruction_angle_err[i][t]=0.1;
	}
    }
  
  
  /*  
      double *MU  = new double [Ntimes];
      double *MUE = new double [Ntimes];
      double *ANE = new double [Ntimes];
  */
  
  
  
  TCanvas *cacca = new TCanvas("cacca","cacca",700,900);  
  cacca->Divide(1,2);
  TLegend *leg_1 = new TLegend(.8,.7,.99,.99);
  TLegend *leg_2 = new TLegend(.8,.7,.99,.99);
  for (int i = 0; i<Nenergies;i++)
    {
      std::vector<double> MU;
      std::vector<double> MUE;
      std::vector<double> ANE;
      std::vector<double> TIMES;
      
      for (int t = 0; t < Ntimes; t++)
	{
	  std::cout<<TOB[t]<<" "<<2+i<<" "<<NumberOfEvents[i][t]<<" "<<modulation_factors[i][t]<<" "<<
	    modulation_factors_err[i][t]<<" "<<reconstruction_angle_err[i][t]<<std::endl;
	  if(modulation_factors[i][t]>0)
	    {
	      TIMES.push_back(TOB[t]);
	      MU.push_back(modulation_factors[i][t]);
	      MUE.push_back(modulation_factors_err[i][t]);
	      ANE.push_back(reconstruction_angle_err[i][t]);
	    }
	}
      std::cout<<"Size: "<<TIMES.size()<<std::endl;
      
      TGraphErrors *gr_mu = new TGraphErrors(TIMES.size(), &TIMES[0], &MU[0], 0, &MUE[0]);
      gr_mu->SetMarkerStyle(20+i);
      gr_mu->SetLineColor(1+i);

      char label1[100];
      char label2[100];
      sprintf(label1, "[%d-%d] keV",i+2,i+3);
      sprintf(label2, "[%d-%d] keV",i+2,i+3);
      
      leg_1->AddEntry(gr_mu,label1,"lp");
      
      TGraph *phi = new TGraph(TIMES.size(), &TIMES[0],&ANE[0]);
      phi->SetMarkerStyle(20+i);
      //      phi->SetLineStyle(2);
      phi->SetLineColor(1+i);
      leg_2->AddEntry(phi,label2,"lp");
      
      if(i==0) 
	{
	  gr_mu->SetTitle("Measured Polarization (\%)");
	  phi->SetTitle("Polarization Angle Estimtated Error (deg)");
	  cacca->cd(1); gr_mu->Draw("alp"); 
	  cacca->cd(2); phi->Draw("alp"); 
	  
	}
      else
      {

	cacca->cd(1); gr_mu->Draw("lp");
	cacca->cd(2); phi->Draw("lp");
      }
    }
  cacca->cd(1); leg_1->Draw();
  cacca->cd(2); leg_2->Draw();
}


