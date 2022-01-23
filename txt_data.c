#include "header/txt_data.h"

#include "header/ui.h"
#include "header/phi.h"
#include "header/tortoise.h"

#include <stdio.h>

void
txt_data_wave(void)
{
  FILE *fp;
  fp = fopen("wave.txt", "w");

  double x = ui_get_a();
  double h = ui_get_h_forward();
  double* phiRea = phi_get_rea();
  double* phiImg = phi_get_img();
  double* phiV = phi_get_v();

  int i;
  for (i = 0; i < NODES; i++)
  {
    fprintf(fp, "%.32f\t%.32f\t%.32f\t%.32f\n", x, phiRea[i], phiImg[i], phiV[i]);
    x += h;
  }
  fclose(fp);
}

/**
  @brief  Compare the calculated value to the analytical one
*/
void
txt_data_turtle()
{
  FILE *fp;
  fp = fopen("turtle.txt", "w");

  /* Calculated value */
  GArray *xy = tortoise_get_xy();
  /* Calculated value. We must invert it to be xy in gnuplot */
  TPoint2D *yx = tortoise_get_yx();
  /* Analytical value */
  tortoise_xy_analytical();
  double *xyAnalytical = tortoise_get_xy_analytical();

  int i;
  for (i = 0; i < NODES; i++)
  {
    fprintf(fp, "%.32f\t%.32f\t%.32f\t%.32f\n", point2D_get_x(g_array_index(xy, TPoint2D, i)), point2D_get_y(g_array_index(xy, TPoint2D, i)), xyAnalytical[i], point2D_get_y(yx[i]));
  }
  fclose(fp);
}

void
txt_data_rt()
{
  FILE *fp;
  fp = fopen("coefficients.txt", "w");

  int nW = ui_get_nW();
  double wMin = ui_get_wMin();
  double hW = ui_get_hW();
  double *R = phi_get_R();
  double *T = phi_get_T();
  double w;

  int i;
  for (i = 0; i < nW; i++)
  {
    w = wMin + hW*i;
    fprintf(fp, "%.32f\t%.32f\t%.32f\t%.32f\n", w, R[i], T[i], R[i] + T[i]);
  }
  fclose(fp);
}

void
txt_data_sigma()
{

  FILE *fp;
  fp = fopen("sigma-l.txt", "w");

  double wMin = ui_get_wMin();
  int nW = ui_get_nW();
  double hW = ui_get_hW();
  double *deltaL = phi_get_deltaL();
  double w;

  int i;
  for (i = 0; i < nW; i++)
  {
    w = wMin + hW*i;
    fprintf(fp, "%.32f\t%.32f\n", w, deltaL[i]);
  }
  fclose(fp);
}
