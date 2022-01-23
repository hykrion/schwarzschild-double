#include "header/tortoise.h"

#include "header/ui.h"

#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>

#define TOR_SYS_DIM 1
#define TOR_MIN_INF -1000000

static GArray *m_xy;
static TPoint2D m_yx[NODES];
static double m_xyAnalytical[NODES];

/**
  @brief  Jacobian xy[]: rS/r(y)^2
*/
static int
jacoXY(double x,
      const double y[],
      double *dfdy,
      double dfdt[],
      void *params)
{
  (void)(x);
  struct tortoise_xyParams *par = (struct tortoise_xyParams*)params;
  double rS = (par->rS);

  gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 1, 1);
  gsl_matrix *m = &dfdy_mat.matrix;
  gsl_matrix_set (m, 0, 0, rS/(y[0]*y[0]));
  dfdt[0] = 0.0;

  return GSL_SUCCESS;
}

/**
  @brief  IV ODE to solve x[y]' that appear in the calculation of V
  @param  t     Independent variable
  @param  y[]   Left side of the system of equations
  @param  sys[] Right side of the system of equations
*/
static int
funcXY(double x,
       const double y[],
       double sys[],
       void* params)
{
  (void)(x);
  struct tortoise_xyParams* par = (struct tortoise_xyParams*)params;
  double rS = (par->rS);
  sys[0] = 1 - rS/y[0];

  return GSL_SUCCESS;
}

/**
  @brief  Jacobian yx[]: -rS/(r^2*(1 - rS/r)^2)
*/
static int
jacoYX(double x,
      const double y[],
      double *dfdy,
      double dfdt[],
      void *params)
{
  (void)(y);
  struct tortoise_xyParams *par = (struct tortoise_xyParams*)params;
  double rS = (par->rS);

  gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 1, 1);
  gsl_matrix *m = &dfdy_mat.matrix;
  gsl_matrix_set (m, 0, 0, -rS/(x*x*(1 - rS/x)*(1 - rS/x)));
  dfdt[0] = 0.0;

  return GSL_SUCCESS;
}

/**
  @brief  IV ODE to solve y[x]' that appear in the calculation of V
  @param  t     Independent variable
  @param  y[]   Left side of the system of equations
  @param  sys[] Right side of the system of equations
*/
static int
funcYX(double x,
       const double y[],
       double sys[],
       void* params)
{
  (void)(y);
  struct tortoise_xyParams* par = (struct tortoise_xyParams*)params;
  double rS = (par->rS);
  sys[0] = 1.0/(1.0 - rS/x);

  return GSL_SUCCESS;
}

/**
  @brief  Compare function to order doubles
*/
static gint
compare_double_gth(gconstpointer a,
                   gconstpointer b)
{
  int result = 0;

  double *aa = (double*)a;
  double *bb = (double*)b;

  if(*aa - *bb > 0)
    result = 1;
  else if(*aa - *bb < 0)
    result = -1;

  return result;
}

/**
  @brief  Compare function to find 'x' in 'yx'
*/
static gint
compare_yx_gth(gconstpointer a,
               gconstpointer b)
{
  int result = 0;

  const TPoint2D *aa = a;
  const TPoint2D *bb = b;

  if(point2D_get_y(*aa) - point2D_get_y(*bb) < 0)
    result = -1;

  return result;
}

/* -----------------------------------------------------------------------------
  PUBLIC
----------------------------------------------------------------------------- */
void
tortoise_init(void)
{
  m_xy = g_array_new(FALSE, FALSE, sizeof(TPoint2D));
}

void
tortoise_destroy(void)
{
  g_array_free(m_xy, FALSE);
}

/**
  @brief  Calculate coordintes 'yx' and 'xy'
*/
void
tortoise_calculate_yx_xy(void)
{
  double a = ui_get_a();
  double b = ui_get_b();
  double rS = ui_get_rS();

  struct tortoise_xyParams torParam = {rS};
  /* Calculate yx from analytic expression */
  tortoise_yx_analytical();
  /* Now that we have y(x) we can use it to get the x'(y) IC */
  double ci = b - rS*log(b/rS - 1);
  ci = tortoise_get_yx_x(b);
  tortoise_xy_integration(b, a, NODES, ci, &torParam);
  tortoise_sort_asc();
}

int
tortoise_xy_integration(double a,
                        double b,
                        int nodes,
                        double ci,
                        void *param)
{
  int status = GSL_SUCCESS;

  struct tortoise_xyParams *par = (struct tortoise_xyParams*)(param);
  double rS = par->rS;
  int n = nodes - 1;
  double h = (b - a)/ n;
  double x0 = a;
  double x1 = x0 + h;
  double epsAbs = 0;
  double epsRel = 1e-6;

  struct tortoise_xyParams params = {rS};

  gsl_odeiv2_system sys = {funcXY, jacoXY, TOR_SYS_DIM, &params};
  const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rk8pd; /* rk2 rk4 rkf45 rk8pd, msbdf msadams rk4imp */
  gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new (&sys, T, h, epsAbs, epsRel);

  double y[1] = {ci};
  int i;
  for (i = 0; i < nodes; i++)
  {
    status = gsl_odeiv2_driver_apply(d, &x0, x1, y);
    x0 = x1;
    x1 = x0 + h;

    if (status != GSL_SUCCESS)
    {
      printf ("error, return value = %d\n", status);
      break;
    }
    /*
    y[0] function value
    y[1] derivative function value
    */
    TPoint2D p2d = {x0, y[0]};
    g_array_append_val(m_xy, p2d);
  }
  gsl_odeiv2_driver_free(d);

  return status;
}

int
tortoise_yx_integration(double a,
                        double b,
                        int nodes,
                        double ci,
                        void *param)
{
  int status = GSL_SUCCESS;

  struct tortoise_xyParams *par = (struct tortoise_xyParams*)(param);
  double rS = par->rS;
  int n = nodes - 1;
  double h = (b - a)/ n;
  double x0 = a;
  double x1 = x0 + h;
  double epsAbs = 0;
  double epsRel = 1e-6;

  struct tortoise_xyParams params = {rS};

  gsl_odeiv2_system sys = {funcYX, jacoYX, TOR_SYS_DIM, &params};
  const gsl_odeiv2_step_type* T = gsl_odeiv2_step_msbdf; /* rk2 rk4 rkf45 rk8pd, msbdf msadams rk4imp */
  gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new (&sys, T, h, epsAbs, epsRel);

  double y[1] = {ci};
  int i;
  for (i = 0; i < nodes; i++)
  {
    status = gsl_odeiv2_driver_apply(d, &x0, x1, y);
    x0 = x1;
    x1 = x0 + h;

    if (status != GSL_SUCCESS)
    {
      printf ("error, return value = %d\n", status);
      break;
    }
    /*
    y[0] function value
    y[1] derivative function value
    */
    TPoint2D p2d = {x0, y[0]};
    m_yx[i] = p2d;
  }
  gsl_odeiv2_driver_free(d);

  return status;
}

/**
  @brief  Reverse the array
*/
void
tortoise_sort_asc(void)
{
  g_array_sort(m_xy, (GCompareFunc)compare_double_gth);
}

GArray*
tortoise_get_xy(void)
{
  return m_xy;
}

TPoint2D*
tortoise_get_yx(void)
{
  return m_yx;
}

/**
  @brief  Search the 'x' in 'yx' that has the more approximate value of 'y'

  NOTE    In that case we have built 'x' so we've got the formula
*/
double
tortoise_get_xy_y(double x)
{
  double a = ui_get_a();
  double h = ui_get_h_forward();
  gint i = (x - a)/h;
  TPoint2D p2d = g_array_index(m_xy, TPoint2D, i);

  return point2D_get_y(p2d);
}

/**
  @brief  Search the 'y' in 'yx' that has the more approximate value of 'x'
*/
double
tortoise_get_yx_x(double y)
{
  double result;
  GArray *tmp = g_array_new(FALSE, FALSE, sizeof(TPoint2D));

  /* Short the array */
  g_array_append_vals(tmp, &m_yx, NODES);
  g_array_sort(tmp, (GCompareFunc)compare_double_gth);
  /* Search for the more approximate value of 'y' */
  guint index;
  TPoint2D target = {0, y};
  gboolean find = g_array_binary_search(tmp, &target, (GCompareFunc)compare_yx_gth, &index);
  /*
  We know that the value is >= but not optimal, so we search for it now that
  we're close
  */
  if(find)
  {
    gint i = index;
    TPoint2D p2d = g_array_index(tmp, TPoint2D, i);
    double yTarget = point2D_get_y(p2d);
    while(yTarget - y > 0 && i >= 0)
    {
      i--;
      p2d = g_array_index(tmp, TPoint2D, i);
      yTarget = point2D_get_y(p2d);
    }
    p2d = g_array_index(tmp, TPoint2D, i);
    result = point2D_get_x(p2d);
  }
  else
  {
    index = NODES - 1;
    TPoint2D p2d = g_array_index(tmp, TPoint2D, index);
    result = point2D_get_x(p2d);
  }

  g_array_free(tmp, FALSE);

  return result;
}

TPoint2D
tortoise_get_yx_i(int i)
{
  return m_yx[i];
}

/**
  @brief  Analytical result (forward)
*/
void
tortoise_xy_analytical(void)
{
  double x = ui_get_a();
  double h = ui_get_h_forward();
  double rS = ui_get_rS();
  int i;
  for(i = 0; i < NODES; i++)
  {
    /*
    NOTE  I'm using yx because if I use xy the graphic is no so good (it should
          be difficult to plot the exponential... I guest)
    */
    double yx = x + rS*log(x/rS - 1);

    /* log(x <= 0) = -inf */
    if(isnan(yx))
      yx = TOR_MIN_INF;

    m_xyAnalytical[i] = yx;
    x += h;
  }
}

double*
tortoise_get_xy_analytical(void)
{
  return m_xyAnalytical;
}

void
tortoise_yx_analytical(void)
{
  double x = ui_get_a();
  double h = ui_get_h_forward();
  double rS = ui_get_rS();
  int i;
  for(i = 0; i < NODES; i++)
  {
    double yx = x + rS*log(x/rS - 1);

    /* log(x <= 0) = -inf */
    if(isnan(yx))
      yx = TOR_MIN_INF;

    TPoint2D p2d = {x, yx};
    m_yx[i] = p2d;
    x += h;
  }

}

TPoint2D*
tortoise_get_yx_analytical(void)
{
  return m_yx;
}
