#ifndef TXYZ_HH
#define TXYZ_HH
#include "MC.h"

class TXYZ{
 public:
  TXYZ(){;}
  TXYZ(double x, double y, double z);
  inline  double X(){return m_X;}
  inline double Y(){return m_Y;}
  inline double Z(){return m_Z;}
  inline void SetXYZ(double x, double y, double z) {m_X=x; m_Y=y; m_Z=z;}
 private:
  double m_X,m_Y,m_Z;
};

#endif
