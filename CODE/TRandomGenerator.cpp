#include "TRandomGenerator.h"


//*******************************************************************************
/// \brief Basic constructor.                                                   *
//*******************************************************************************
TRandomGenerator::TRandomGenerator()
{
  RandomGenerator = new TRandom();
}


//*******************************************************************************
/// \brief Basic destructor.                                                    *
//*******************************************************************************
TRandomGenerator::~TRandomGenerator()
{

}


//*******************************************************************************
/// \brief Chooses a random seed for the random number generator.               *
/// The routine uses time() as to set the seed.                                 *
//*******************************************************************************
void TRandomGenerator::SetRandomSeed()
{
  UInt_t Seed = (UInt_t)time('\0');
  SetSeed(Seed);
}


//*******************************************************************************
/// \brief Sets the seed of the random number generator.                        *
//*******************************************************************************
void TRandomGenerator::SetSeed(UInt_t SEED)
{
  RandomGenerator->SetSeed(SEED);
}
