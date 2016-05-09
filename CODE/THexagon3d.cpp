#include "THexagon3d.h"
#include "TPolyLine3D.h"

//THexagon3d::THexagon3d(Int_t NumHits, HitPixel *fHits ) {

THexagon3d::THexagon3d(double X, double Y, double Z, double Height) {
  Size  = Height;
  // Find Hexagon vertex points
  x[0] = ( X + Size* 0.5 );
  y[0] = ( Y + Size* 0.8660254 );
  x[1] = ( X + Size );
  y[1] = ( Y );
  x[2] = ( X + Size* 0.5 );
  y[2] = ( Y - Size* 0.8660254 );
  x[3] = ( X - Size* 0.5 );
  y[3] = ( Y - Size* 0.8660254 );
  x[4] = ( X - Size );
  y[4] = ( Y );
  x[5] = ( X - Size* 0.5 );
  y[5] = ( Y + Size* 0.8660254 );
  z    =   Z;
}

void THexagon3d::Draw(int color, int LineWidth) {
  TPolyLine3D *hexagon = new TPolyLine3D (7);
  TPolyLine3D *face1 = new TPolyLine3D (3);
  TPolyLine3D *face2 = new TPolyLine3D (3);
  TPolyLine3D *face3 = new TPolyLine3D (3);
  TPolyLine3D *face4 = new TPolyLine3D (3);
  TPolyLine3D *face5 = new TPolyLine3D (3);
  TPolyLine3D *face6 = new TPolyLine3D (3);
  
  hexagon->SetPoint(0,x[ 0 ], y[ 0 ], z);
  hexagon->SetPoint(1,x[ 1 ], y[ 1 ], z);
  hexagon->SetPoint(2,x[ 2 ], y[ 2 ], z);
  hexagon->SetPoint(3,x[ 3 ], y[ 3 ], z);
  hexagon->SetPoint(4,x[ 4 ], y[ 4 ], z);
  hexagon->SetPoint(5,x[ 5 ], y[ 5 ], z);
  hexagon->SetPoint(6,x[ 0 ], y[ 0 ], z);
  hexagon->SetLineColor(color);
  hexagon->SetLineWidth(LineWidth);
  
  hexagon-> Draw();
  
  if(z>0)
    {
      face1->SetPoint(0,x[ 0 ], y[ 0 ], z);
      face1->SetPoint(1,x[ 0 ], y[ 0 ], 0);
      face1->SetPoint(2,x[ 1 ], y[ 1 ], 0);
      face2->SetPoint(0,x[ 1 ], y[ 1 ], z);
      face2->SetPoint(1,x[ 1 ], y[ 1 ], 0);
      face2->SetPoint(2,x[ 2 ], y[ 2 ], 0);
      face3->SetPoint(0,x[ 2 ], y[ 2 ], z);
      face3->SetPoint(1,x[ 2 ], y[ 2 ], 0);
      face3->SetPoint(2,x[ 3 ], y[ 3 ], 0);
      face4->SetPoint(0,x[ 3 ], y[ 3 ], z);
      face4->SetPoint(1,x[ 3 ], y[ 3 ], 0);
      face4->SetPoint(2,x[ 4 ], y[ 4 ], 0);
      face5->SetPoint(0,x[ 4 ], y[ 4 ], z);
      face5->SetPoint(1,x[ 4 ], y[ 4 ], 0);
      face5->SetPoint(2,x[ 5 ], y[ 5 ], 0);
      face6->SetPoint(0,x[ 5 ], y[ 5 ], z);
      face6->SetPoint(1,x[ 5 ], y[ 5 ], 0);
      face6->SetPoint(2,x[ 0 ], y[ 0 ], 0);
      face1->SetLineColor(color);
      face2->SetLineColor(color);
      face3->SetLineColor(color);
      face4->SetLineColor(color);
      face5->SetLineColor(color);
      face6->SetLineColor(color);
      face1->SetLineWidth(LineWidth);
      face2->SetLineWidth(LineWidth);
      face3->SetLineWidth(LineWidth);
      face4->SetLineWidth(LineWidth);
      face5->SetLineWidth(LineWidth);
      face6->SetLineWidth(LineWidth);
      face1-> Draw();
      face2-> Draw();
      face3-> Draw();
      face4-> Draw();
      face5-> Draw();
      face6-> Draw();
    }
}
