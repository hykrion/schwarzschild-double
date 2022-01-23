#ifndef TORTOISE_H_INCLUDED
#define TORTOISE_H_INCLUDED

#include "point2D.h"
#include <glib.h>

struct tortoise_xyParams
{
  double rS;
};

void tortoise_init(void);
void tortoise_destroy(void);

void tortoise_calculate_yx_xy(void);
int
tortoise_xy_integration(double a,
                        double b,
                        int nodes,
                        double ci,
                        void *param);
int
tortoise_yx_integration(double a,
                        double b,
                        int nodes,
                        double ci,
                        void *param);
void tortoise_sort_asc(void);
GArray* tortoise_get_xy(void);
TPoint2D* tortoise_get_yx(void);
double tortoise_get_xy_y(double x);
double tortoise_get_yx_x(double y);
TPoint2D tortoise_get_yx_i(int i);
void tortoise_xy_analytical(void);
double* tortoise_get_xy_analytical(void);
void tortoise_yx_analytical(void);
TPoint2D* tortoise_get_yx_analytical(void);

#endif // TORTOISE_H_INCLUDED
