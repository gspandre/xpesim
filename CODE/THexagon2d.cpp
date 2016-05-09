#include "THexagon2d.h"
#include "TPolyLine.h"

THexagon2d::THexagon2d(double X, double Y, double Height) {
  Size  = Height;
  // Find THexagon vertex points
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
}

void THexagon2d::Draw(int color, int LineWidth) {
  TPolyLine *hexagon = new TPolyLine (7);
  hexagon->SetPoint(0,x[ 0 ], y[ 0 ]);
  hexagon->SetPoint(1,x[ 1 ], y[ 1 ]);
  hexagon->SetPoint(2,x[ 2 ], y[ 2 ]);
  hexagon->SetPoint(3,x[ 3 ], y[ 3 ]);
  hexagon->SetPoint(4,x[ 4 ], y[ 4 ]);
  hexagon->SetPoint(5,x[ 5 ], y[ 5 ]);
  hexagon->SetPoint(6,x[ 0 ], y[ 0 ]);
  hexagon->SetLineColor(color);
  hexagon->SetLineWidth(LineWidth);
  
  hexagon-> Draw();
}
