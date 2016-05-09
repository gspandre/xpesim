#ifndef THEXAGON3D_HH
#define THEXAGON3D_HH

class THexagon3d {
  
 public:
  
  THexagon3d(double X, double Y, double Z, double Height);
  ~THexagon3d(){;}
  void Draw(int color, int LineWidth);
  Float_t HitX, HitY, Size;
  Float_t x[ 6 ];
  Float_t y[ 6 ];
  Float_t z;
};

#endif
