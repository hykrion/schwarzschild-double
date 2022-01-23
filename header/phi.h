#ifndef PHI_H_INCLUDED
#define PHI_H_INCLUDED

struct TPhiParams
{
  double w;
  double l;
  double u; /* For the system of differential equations */
  double v; /* idem */
};

void phi_wave_yx_xy(double l);
void phi_rea_fwd(void);
void phi_img_fwd(void);
void phi_calculate_RT(double w, int i);
void phi_calculate_v(void);
double* phi_get_rea(void);
double* phi_get_img(void);
double* phi_get_v(void);
double* phi_get_R(void);
double* phi_get_T(void);
double* phi_get_deltaL(void);

int
phi_integration_rea(double a,
                     double b,
                     int nodes,
                     void *intParams);
int
phi_integration_img(double a,
                     double b,
                     int nodes,
                     void *intParams);
void
phi_calculate_sigma_l(double l);

#endif // PHI_H_INCLUDED
