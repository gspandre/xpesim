#include "TSource.h"

//*******************************************************************************
/// \brief Basic constructor.                                                   *
/// Polarization angle to be given in degrees.                                  *
//*******************************************************************************
TSource::TSource(Double_t ENERGY, Double_t POLARIZATION_ANGLE,
		 Double_t POLARIZATION_DEGREE, TRandom *RND, Bool_t VERBOSE)
{
  // Make constructor parameters class members.
  Energy             = ENERGY;
  PolarizationAngle  = POLARIZATION_ANGLE * kPI/180.0;
  PolarizationDegree = POLARIZATION_DEGREE;
  Verbose            = VERBOSE;
  rnd = RND;
  ThetaDistribution =
    new TF1("ThetaDistribution","[0]*pow(sin(x), 3.0)/pow(1.0-[1]*cos(x), 4.0)", 0.0, kPI);
  ThetaDistribution->SetParameter(0, 1.0);
  PhiDistribution =
    new TF1("PhiDistribution","[0] + [1]*pow((cos(x-[2])), 2.0)", 0.0, 2.0*kPI);  /// gradodiPolariz=[1]/(2*[0]+[1])
  PhiDistribution->SetParameter(0, (1-PolarizationDegree));
  PhiDistribution->SetParameter(1, 2*PolarizationDegree); 
  PhiDistribution->SetParameter(2, PolarizationAngle);
}


//*******************************************************************************
/// \brief Basic destructor.                                                    *
//*******************************************************************************
TSource::~TSource()
{
  delete ThetaDistribution;
  delete PhiDistribution;
}


//*******************************************************************************
/// \brief Returns the energy of the photon emitted by the source.              *
//*******************************************************************************
Double_t TSource::GetEnergy()
{
  return Energy;
}


//*******************************************************************************
/// \brief Sets the polarization angle for the source.                          *
/// Angle to be provided in degrees.                                            *
//*******************************************************************************
void TSource::SetPolarizationAngle(Double_t POLARIZATION_ANGLE)
{
  PhiDistribution->SetParameter(2, POLARIZATION_ANGLE*kPI/180.0);
}


//*******************************************************************************
/// \brief Returns the energy of the photon emitted by the source.              *
//*******************************************************************************
void TSource::SetPolarizationDegree(Double_t POLARIZATION_DEGREE)
{
  PhiDistribution->SetParameter(0, (1-POLARIZATION_DEGREE));
  PhiDistribution->SetParameter(1, 2*POLARIZATION_DEGREE);
}

void TSource::SetPattern(TString Set)
{
  PatternType = Set;
}

std::pair<double,double> TSource::GetConversionXY()
{
  std::pair<double,double> XY = make_pair(0.,0.);
  if (PatternType=="GAUS")
    {
      double X=rnd->Gaus(0.,0.05);
      double Y=rnd->Gaus(0.,0.05);
	
      XY = make_pair(X,Y);
    }
  else if(PatternType=="UNIFORM")
    {
      double X=rnd->Uniform(-0.007,0.007);
      double Y=rnd->Uniform(-0.007,0.007);
      XY = make_pair(X,Y);
    }
  else XY = make_pair(0.,0.);
  return XY;
}
