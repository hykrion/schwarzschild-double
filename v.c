#include "header/v.h"

#include "header/globals.h"
#include "header/ui_parameters.h"

#include <math.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#define SYS_DIM 1

/**
  @brief  ODE IV de 1er grado para resolver x[y]' que aparece en el cálculo del
          potencial V.
  @param  t     Variable independiente
  @param  y[]   Parte izq del sistema de ecuaciones de primer grado
  @param  sys[] Parte dcha del sistema de ecuaciones de primer grado
*/
int
funcX (double t,
       const double y[],
       double sys[],
       void* params)
{
  double rS = getRS();
  /* Evitar aviso de variable no usada */
  (void)(t);
  (void)(params);

  double x = y[0];

  sys[0] = 1 - rS/x;

  return GSL_SUCCESS;
}

int
jacoX(double t,
      const double y[],
      double *dfdy,
      double dfdt[],
      void *params)
{
  double x = y[0];

  (void)(t); /* avoid unused parameter warning */
  struct xyParams *par = (struct xyParams*)params;
  double rS = (par->rS);

  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 1, 1);
  gsl_matrix *m = &dfdy_mat.matrix;
  gsl_matrix_set (m, 0, 0, rS/(x*x));
  dfdt[0] = 0.0;

  return GSL_SUCCESS;
}

int
xyIntegration()
{
  int status = GSL_SUCCESS;

  double a = getA();
  double rS = getRS();
  double h = getH();
  int n = getN();
  /*
  h = 1e-6;
  n = 10000;
  */
  /*
  DEBUG
    Con h = 1e-6 y n = 2000 ya tengo suficientes puntos para mostrar el potencial
    en el rango [-20, 20].
    Extrañamente... con h = 1e-6, n = 4000 y [-40, 40] saca resultados, pero si
    pongo [-50, 50] no ¿¿por qué si el paso es fijo??

    Con estos valores lo grafica bien en [-150, 150] y paso fijo
  h = 1e-6;
  n = 10000;
  */
  a = -60;
  double x0 = a;
  double x1 = x0 + h;
  double epsAbs = 0;
  double epsRel = 1e-6;

  /*
  Dato del problema: CI en -inf
  NOTE  Hemos de sumar 1e-9 si queremos usar infinitos que vayan más allá de unas
        40 unidades desde el origen porque de otro modo, no se integra.
  */
  double y0 = rS*(1 + exp(x0/rS - 1));
  struct xyParams params = {rS};

  gsl_odeiv2_system sys = {funcX, NULL, SYS_DIM, &params};
  const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rkf45; /* rk8pd gsl_odeiv2_step_rkf45 */
  gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new (&sys, T, h, epsAbs, epsRel);

  double yy[1] = {y0};
  int i;
  for (i = 0; i < NODES; i++)
  {
    /*
    NOTE  Si epsAbs es muy bajo... da error

    Fijo:

    status = gsl_odeiv2_driver_apply_fixed_step(d, &x0, h, n, yy);

    Adaptativo:

    status = gsl_odeiv2_driver_apply(d, &x0, x1, yy);
    x0 = x1;
    x1 = x0 + h;
    */
    status = gsl_odeiv2_driver_apply(d, &x0, x1, yy);
    x0 = x1;
    x1 = x0 + h;

    if (status != GSL_SUCCESS)
    {
      printf ("error, return value = %d\n", status);
      break;
    }
    /*
    yy[0] valor de la función
    yy[1] valor de la derivada de la función (en este caso no existe)
    */
    double foo = yy[0];
    xy[i] = yy[0];
  }
  gsl_odeiv2_driver_free(d);

  return status;
}

/*
  @brief  Ecuación de la barrera/pozo potencial. Es necesario haber calculado
          antes xy[]
*/
double
V(double y)
{
  double a = getA();
  double ll = getLL();
  double rS = getRS();
  double h = getH();

  /* Nos pasan un double... pero necesitamos un int */
  int i = (y - a) / h;

  return (ll*(1 + ll)/(xy[i]*xy[i]) + rS/(xy[i]*xy[i]*xy[i]))*(1 - rS/xy[i]);
  /*return (ll*(1 + ll)*(rS - xy[i])/(xy[i]*xy[i]*xy[i]));*/
}
