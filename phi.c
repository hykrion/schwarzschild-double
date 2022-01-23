#include "header/phi.h"

#include "header/ui.h"
#include "header/tortoise.h"

#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex_math.h>

#define SYS_DIM 2
/* Arbitrary indexes to solve the system of equations */
#define PHI_SYS_IND_0 15
#define PHI_SYS_IND_1 20

static double m_phiRea[NODES];
static double m_phiImg[NODES];
static double m_T[NODES];
static double m_R[NODES];
static double m_v[NODES];
static double m_deltaL[NODES];

/*
  @brief  Potential well / barrier. You need to calculate xy[] before
*/
static double
phi_v(double x)
{
  double l = ui_get_l();
  double rS = ui_get_rS();

  double y = tortoise_get_xy_y(x);

  return (l*(1 + l)/pow(y, 2) + rS/pow(y, 3))*(1 - rS/y);
}

/**
  @brief  ODE system

  @param  r: independent variable
  @param  y: system variables for the 1st order system
  @param  dydx: system derivatives
  @param  params: equation parameters
*/
static int
funcPhi (double r,
         const double y[],
         double dydx[],
         void* params)
{
  struct TPhiParams *par = (struct TPhiParams *)params;
  double w = (par->w);

  double u = y[0];
  double v = y[1];

  dydx[0] = v;
  dydx[1] = (phi_v(r) - w*w)*u;

  return GSL_SUCCESS;
}

/**
  @brief  Integrate phi

  @param  IN, a: interval inferior limit
  @param  IN, b: interval superior limit
  @param  IN, nodes: #nodes
  @param  IN, ci: initial condition
  @param  IN, intParams: params needed for the integration
  @param  IN, OUT phi: array with values of phi

  @return integration status
*/
static int
phi_integration(double a,
                double b,
                int nodes,
                void *intParams,
                double *phi)
{
  int status = GSL_SUCCESS;

  struct TPhiParams *intPar = (struct TPhiParams*)intParams;
  double w = (intPar->w);
  double l = (intPar->l);
  double u = (intPar->u);
  double v = (intPar->v);

  double h = (b - a)/(nodes - 1);

  double x0 = a;
  double x1 = x0 + h;
  double epsAbs = 0;
  double epsRel = 1e-6;

  struct TPhiParams par = {w, l, 0, 0};

  gsl_odeiv2_system sys = {funcPhi, NULL, SYS_DIM, &par};
  const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rk8pd; /* rkf45 rk8pd , rk4imp bsimp msadams msbdf*/
  gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new (&sys, T, h, epsAbs, epsRel);

  double y[2] = {u, v};
  int i;
  for (i = 0; i < nodes; i++)
  {
    status = gsl_odeiv2_driver_apply(d, &x0, x1, y);
    x0 = x1;
    x1 = x0 + h;

    if (status != GSL_SUCCESS)
    {
      printf ("error, return value = %d\n", status);
      printf("x = %f, y = %f\n", x0, y[0]);
      break;
    }
    else
      phi[i] = y[0];
  }
  gsl_odeiv2_driver_free(d);

  return status;
}

/* -----------------------------------------------------------------------------
  PUBLIC
----------------------------------------------------------------------------- */

/**
  @brief  Calculate coeficients R, T within an interval [a, b], for an angular
          frequency range [wMin, wMax] for a particular angular momentum solving
          the ODE x'(y)

  @param  IN, ci: condición inicial
*/
void
phi_wave_yx_xy(double l)
{
  double a = ui_get_a();
  double b = ui_get_b();
  double wMin = ui_get_wMin();
  double hW = ui_get_hW();
  int nW = ui_get_nW();
  double w;

  /*
  Far away points from the barrier where we solve the system of equations to get
  the R & T values
  */
  int i;
  for(i = 0; i < nW; i++)
  {
    /*
    Calculate 'q' from the Sturm-Liouville equation and seek indexes at the
    barrier
    */
    w = wMin + hW*i;

    /* Backward integration */

    /* Real part */
    double u = cos(w*b);
    double v = -w*sin(w*b);
    struct TPhiParams intParams = {w, l, u, v};
    phi_integration_rea(b, a, NODES, &intParams);

    /* Imaginary part */
    u = sin(w*b);
    v = w*cos(w*b);
    intParams.u = u;
    intParams.v = v;
    phi_integration_img(b, a, NODES, &intParams);

    /* Sort the arrays as if they were a forward integration */
    phi_rea_fwd();
    phi_img_fwd();

    phi_calculate_RT(w, i);
  }
}

/**
  @brief  phi integration, real part

  @param  IN, a: interval inferior limit
  @param  IN, b: interval superior limit
  @param  IN, nodes: #nodes
  @param  IN, ci: initial condition
  @param  IN, intParams: params for the integration
  @param  IN, OUT phi: array with phi values

  @return integration state
*/
int
phi_integration_rea(double a,
                    double b,
                    int nodes,
                    void *intParams)
{
  return phi_integration(a, b, nodes, intParams, m_phiRea);
}

/**
  @brief  phi integration, imaginary part

  @param  IN, a: interval inferior limit
  @param  IN, b: interval superior limit
  @param  IN, nodes: #nodes
  @param  IN, ci: initial condition
  @param  IN, intParams: params for the integration
  @param  IN, OUT phi: array with phi values

  @return integration state
*/
int
phi_integration_img(double a,
                    double b,
                    int nodes,
                    void *intParams)
{
  return phi_integration(a, b, nodes, intParams, m_phiImg);
}

/**
  @brief  Reverse phi vector
*/
void
phi_reverse(double *phi)
{
  double tmp[NODES];
  int i;
  for(i = 0; i < NODES; i++)
    tmp[i] = phi[NODES - 1 - i];
  for(i = 0; i < NODES; i++)
    phi[i] = tmp[i];
}
/**
  @brief  Short phiRea as if it were a forward integration
*/
void
phi_rea_fwd(void)
{
  phi_reverse(m_phiRea);
}

/**
  @brief  Short phiImg as if it were a forward integration
*/
void
phi_img_fwd(void)
{
  phi_reverse(m_phiImg);
}

/**
  @brief  Calculate R & T solving the system equations
*/
void
phi_calculate_RT(double w,
                 int i)
{
  double a = ui_get_a();
  double h = ui_get_h_forward();
  double B0 = m_phiRea[PHI_SYS_IND_0];
  double B1 = m_phiImg[PHI_SYS_IND_0];
  double B2 = m_phiRea[PHI_SYS_IND_1];
  double B3 = m_phiImg[PHI_SYS_IND_1];

  /* Leave constants with positive values... */
  double k1 = -1*(a + PHI_SYS_IND_0*h);
  double k2 = -1*(a + PHI_SYS_IND_1*h);
  gsl_complex z1 = gsl_complex_polar(1, k1*w);
  gsl_complex z2 = gsl_complex_polar(1, k2*w);

  double e0r = GSL_REAL(z1);
  double e0i = GSL_IMAG(z1);
  double e1r = GSL_REAL(z2);
  double e1i = GSL_IMAG(z2);

  double Ar;
  double Ai;
  double Br;
  double Bi;
  gsl_complex AA;
  gsl_complex BB;

  double A[] = { e0r, e0i, e0r, -e0i,
                -e0i, e0r, e0i, e0r,
                 e1r, e1i, e1r, -e1i,
                -e1i, e1r, e1i, e1r };

  double B[] = { B0,
                 B1,
                 B2,
                 B3 };
  gsl_matrix_view m = gsl_matrix_view_array (A, 4, 4);
  gsl_vector_view b = gsl_vector_view_array (B, 4);
  gsl_vector *x = gsl_vector_alloc (4);
  int signum;
  gsl_permutation *permutation = gsl_permutation_alloc (4);
  gsl_linalg_LU_decomp (&m.matrix, permutation, &signum);
  gsl_linalg_LU_solve (&m.matrix, permutation, &b.vector, x);

  Ar = gsl_vector_get(x, 0);
  Ai = gsl_vector_get(x, 1);
  Br = gsl_vector_get(x, 2);
  Bi = gsl_vector_get(x, 3);
  AA = gsl_complex_rect(Ar, Ai);
  BB = gsl_complex_rect(Br, Bi);

  m_T[i] = 1.0 / gsl_complex_abs2(AA);
  m_R[i] = gsl_complex_abs2(BB) / gsl_complex_abs2(AA);

  gsl_permutation_free(permutation);
  gsl_vector_free(x);
}

/**
  @brief  Potential values (for plotting)
*/
void
phi_calculate_v(void)
{
  double a = ui_get_a();
  double h = ui_get_h_forward();

  double x;
  int i;
  for (i = 0; i < NODES; i++)
  {
    x = a + i*h;
    m_v[i] = phi_v(x);
  }
}

double*
phi_get_rea(void)
{
  return m_phiRea;
}

double*
phi_get_img(void)
{
  return m_phiImg;
}

double*
phi_get_v(void)
{
  return m_v;
}

double*
phi_get_R(void)
{
  return m_R;
}

double*
phi_get_T(void)
{
  return m_T;
}

/**
  @brief  Calculate sigma_l for different angular momentums
*/
void
phi_calculate_sigma_l(double l)
{
  double wMin = ui_get_wMin();
  double nW = ui_get_nW();
  double hW = ui_get_hW();
  double *R = phi_get_R();
  double w;

  int i;
  for(i = 0; i < nW; i++)
  {
    w = wMin + hW*i;
    m_deltaL[i] += (M_PI*(2*l + 1)*(1 - R[i]))/(w*w);
  }
}

double*
phi_get_deltaL(void)
{
  return m_deltaL;
}
