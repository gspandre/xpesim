#ifndef THEXAGON2D_HH
#define THEXAGON2D_HH

#include "MC.h"
class THexagon2d {
  
 public:
  
  THexagon2d(double X, double Y, double Height);
  ~THexagon2d(){;}
  void Draw(int color, int LineWidth);
  Float_t HitX, HitY, Size;
  Float_t x[6];
  Float_t y[6];
};

#endif
