#ifndef TRANDOMGENERATOR_HH
#define TRANDOMGENERATOR_HH

#include "MC.h"
class TRandomGenerator{
 public:
  TRandomGenerator();
  ~TRandomGenerator();
  inline TRandom* GetRandomGenerator() {return RandomGenerator;}
  void            SetRandomSeed();
  void            SetSeed(UInt_t SEED);

 private:
  TRandom* RandomGenerator;
};

#endif
