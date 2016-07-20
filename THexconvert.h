#ifndef THEXCONV_HH
#define THEXCONV_HH

/* 
define struct for coord types 
Copied from 
http://www.redblobgames.com/grids/hexagons/codegen/output/lib.cc
Documented in 
http://www.redblobgames.com/grids/hexagons/implementation.htm
*/

#define SQRT3 1.732

// a generic point in space
//struct Point
//{
//    const double x;
//    const double y;
//    Point(double x_, double y_): x(x_), y(y_) {}
//};

// Cube/Hex coordinate
// Since cube is just and extension of hex we use only this one
struct Hex
{
    const int q;
    const int r;
    const int s;
    Hex(int q_, int r_, int s_): q(q_), r(r_), s(s_) {}
};

// As above, but allowing float index
struct FractionalHex
{
    const double q;
    const double r;
    const double s;
    FractionalHex(double q_, double r_, double s_): q(q_), r(r_), s(s_) {}
};

// Offset coordinate system, the one used in the ASIC
struct OffsetCoord
{
    const int col;
    const int row;
    OffsetCoord(int col_, int row_): col(col_), row(row_) {}
};

// Function to round Cube/Hex to int
Hex hex_round(FractionalHex h)
{
    int q = int(round(h.q));
    int r = int(round(h.r));
    int s = int(round(h.s));
    double q_diff = abs(q - h.q);
    double r_diff = abs(r - h.r);
    double s_diff = abs(s - h.s);
    if (q_diff > r_diff and q_diff > s_diff)
    {
        q = -r - s;
    }
    else
        if (r_diff > s_diff)
        {
            r = -q - s;
        }
        else
        {
            s = -q - r;
        }
    return Hex(q, r, s);
};

// Conversion functions cube <-> offset
// Let's use only "even-r" horizontal layout

Hex offset_to_hex(OffsetCoord h)
{
    int q = h.col - int((h.row + (h.row & 1)) / 2);
    int r = h.row;
    int s = -q - r;
    return Hex(q, r, s);
};

OffsetCoord hex_to_offset(Hex h)
{
    int col = h.q + int((h.r + (h.r & 1)) / 2);
    int row = h.r;
    return OffsetCoord(col, row);
};


// Convert a location in space it into a hex grid coordinate
FractionalHex location_to_hex(double x, double y, double size)
{
  double q = (x * SQRT3/3 - y / 3.)/size;
  double r = y * 2./3 / size;
  return FractionalHex(q, r, -q - r);
};
  

OffsetCoord location_to_offset(double x, double y, double size)
{
  FractionalHex h = location_to_hex(x, y, size);
  Hex hh = hex_round(h);
  return hex_to_offset(hh);
};

#endif
