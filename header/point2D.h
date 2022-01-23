#ifndef POINT2D_H_INCLUDED
#define POINT2D_H_INCLUDED

typedef struct
{
  double x;
  double y;
}TPoint2D;

double point2D_get_x(TPoint2D p2d);
double point2D_get_y(TPoint2D p2d);

#endif // POINT2D_H_INCLUDED
