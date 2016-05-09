#include "TPhoton.h"
#include "TMath.h"

//*******************************************************************************
/// \brief Basic constructor.                                                   *
//*******************************************************************************
TPhoton::TPhoton(TSource *SOURCE, TXYZ CONVERSION_POINT,
		 TRandom *RND, Bool_t VERBOSE)
{
  // Make constructor parameters class members.
  Source             = SOURCE;
  ConversionPoint    = CONVERSION_POINT;
  Verbose            = VERBOSE;
  rnd = RND;
  Energy = Source->GetEnergy();  
}


//*******************************************************************************
/// \brief Basic destructor.                                                    *
//*******************************************************************************
TPhoton::~TPhoton()
{

}


//*******************************************************************************
/// \brief Returns the Theta angle of the photoelectron.                        *
//*******************************************************************************
Double_t TPhoton::GetPhotoelectronTheta()
{
  Double_t Beta = sqrt(1.0 - pow((Energy/ELECTRON_MASS + 1), -2.0));//sqrt 23/4!!!!
  Source->GetThetaDistribution()->SetParameter(1, Beta);
  return TMath::Pi()-Source->GetThetaDistribution()->GetRandom();//P-()per correggere sistema riferimento23/4
}


//*******************************************************************************
/// \brief Returns the Phi angle of the photoelectron.                          *
//*******************************************************************************
Double_t TPhoton::GetPhotoelectronPhi()
{
  return Source->GetPhiDistribution()->GetRandom();
}


//*******************************************************************************
/// \brief Returns the Theta angle of the Auger electron.                       *
//*******************************************************************************
Double_t TPhoton::GetAugerElectronTheta()
{
  return acos(rnd->Uniform(-1,1));
}


//*******************************************************************************
/// \brief Returns the Phi angle of the Auger electron.                         *
//*******************************************************************************
Double_t TPhoton::GetAugerElectronPhi()
{
  return rnd->Uniform(0, 2*kPI);
}


