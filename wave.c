#include "header/wave.h"

#include "header/globals.h"
#include "header/ui_parameters.h"
#include "header/v.h"
#include "header/numeric_integration.h"

#include <string.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex_math.h>

static int m_systemIndex0 = 6;
static int m_systemIndex1 = 10;

/* -----------------------------------------------------------------------------
  Funciones de ayuda
----------------------------------------------------------------------------- */
/**
  @brief  Encontrar los valores R y T con nuestro sistema de ecuaciones
*/
void
findRT(double B0,
       double B1,
       double B2,
       double B3,
       double k,
       double* R2,
       double* T2)
{
  /*
  Ángulos -19.76 y -19.6 con rango [-20, 20] y 1000 subintervalos (h = 0.04)
  Los puntos fijos escogidos son x6 y x10, es decir -20 + 0.04*6 y -20 + 0.04*10
  gsl_complex z1 = gsl_complex_polar(1, 19.76*k);
  gsl_complex z2 = gsl_complex_polar(1, 19.6*k);
  Con 10000 subintervalos tenemos -19.976 y -19.96

  Si cogemos un rango de [-150, 150] y 10000 subintervales (h = 0.03), con
  nuestros puntos fijos (x6 y x10) tendríamos unos ángulos de -149.82 y -149.7
  gsl_complex z1 = gsl_complex_polar(1, 149.82*k);
  gsl_complex z2 = gsl_complex_polar(1, 149.7*k);
  */
  /*
  TODO Esto ya no me gusta... debería tener una clase que controle los
  parámetros del GUI
  */
  double aa = getA();
  double h = getH();

  /* Dejamos las constantes en valor positivo */
  double k1 = -1*(aa + m_systemIndex0*h);
  double k2 = -1*(aa + m_systemIndex1*h);
  gsl_complex z1 = gsl_complex_polar(1, k1*k);
  gsl_complex z2 = gsl_complex_polar(1, k2*k);

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

  /* Ángulos -19.76 y -19.6 *
  double A[] = { 0.613395, 0.789776, 0.613395, -0.789776,
                 -0.789776, 0.613395, 0.789776, 0.613395,
                 0.731386, 0.681964, 0.731386, -0.681964,
                 -0.681964, 0.731386, 0.681964, 0.731386 };
*/
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

  *T2 = 1.0 / gsl_complex_abs2(AA);
  *R2 = gsl_complex_abs2(BB) / gsl_complex_abs2(AA);

  gsl_permutation_free(permutation);
  gsl_vector_free(x);
}

/* -----------------------------------------------------------------------------
  Funciones públicas
----------------------------------------------------------------------------- */

/**
  @brief  Valores del potencial. Así podemos graficarlo
*/
void
potentialValues(double vx[])
{
  double a = getA();
  double h = getH();
  int n = getN();

  double x;
  int i;

  for (i = 0; i <= n; i++)
  {
    x = a + i*h;
    vx[i] = V(x);
  }
}

/**
  @brief  delta_l No sé qué es, pero Gonzalo lo grafica
*/
double
deltaL(double Rw,
       double w)
{
  double ll = getLL();

  return (M_PI*(2*ll + 1)*(1 - Rw))/(w*w);
}

/**
  TODO
*/
void wave()
{
  int status = GSL_SUCCESS;

  double a = getA();
  double b = getB();
  double h = getH();
  int n = getN();
  double wMin = getWMin();
  double hW = getHW();
  int nW = getNW();
  double ll = getLL();
  double rS = getRS();

  double ur[NODES];
  double ui[NODES];
  double qf[NODES];
  double qb[NODES];
  double s[NODES];

  double R2;
  double T2;
  double R2T2;
  /*
  Puntos alejados de la barrera donde resolveremos el sistema de ecuaciones
  para obtener los valores de R y T
  */
  double w;

  int i;
  for(i = 0; i < nW; i++)
  {
    /* Limpiar buffers */
    memset(uRe, 0, NODES*sizeof(double));
    memset(uIm, 0, NODES*sizeof(double));
    memset(ur, 0, NODES*sizeof(double));
    memset(ui, 0, NODES*sizeof(double));
    memset(qf, 0, NODES*sizeof(double));
    memset(qb, 0, NODES*sizeof(double));
    memset(s, 0, NODES*sizeof(double));
    memset(xy, 0, NODES*sizeof(double));

    /* Calcular 'q' de la ecuación y buscar índices de la barrera */
    w = wMin + hW*i;

    /* Tenemos que integrar porque solo conocemos su derivada y una CI en -inf */
    status = xyIntegration();

    if(status != GSL_SUCCESS)
    {
      int foo = 0;
    }
    double x;
    int j;
    for (j = 0; j <= n; j++)
    {
      x = a + j*h;
      qf[j] = w*w - (ll*(1 + ll)*(rS - xy[j]))/(xy[j]*xy[j]*xy[j]);
      qb[n - j] = qf[j];
    }

    /* ---------------------------------------------------------------------------
    Integración hacia atrás para, a partir del valor de B calculado obtener el
    valor calculado de A.
    --------------------------------------------------------------------------- */
    /*
    backguardIntegration(b, m, planck, k, ur, ui, qb, s, h, n);
    int status = backguardIntegrationRK(b, ur, ui, qb, s, h, n);
    */

    /*
      Parte real

      Convertimos la ecuación de 2o grado en un sistema: u, v
    printf("Parte real\n");
    */
    double u = cos(w*b);
    double v = -w*sin(w*b);
    int status = backguardIntegrationRKSchwarzschild(b, u, v, ur, -h, w);

    if(status != GSL_SUCCESS)
    {
      int foo = 0;
    }
    /* TODO Manejar status */
    /*
      Parte imaginaria

      Convertimos la ecuación de 2o grado en un sistema: u, v
    */
    u = sin(w*b);
    v = w*cos(w*b);
    status = backguardIntegrationRKSchwarzschild(b, u, v, ui, -h, w);

    if(status != GSL_SUCCESS)
    {
      int foo = 0;
    }

    for(j = 0; j <= n; j++)
    {
      uRe[j] = ur[n - j];
      uIm[j] = ui[n - j];
    }

    findRT(uRe[m_systemIndex0], uIm[m_systemIndex0], uRe[m_systemIndex1], uIm[m_systemIndex1], w, &R2, &T2);
    RR[i] = R2;
    TT[i] = T2;
    /* Por definición debe ser 1 */
    R2T2 = R2 + T2;
    RRTT[i] = R2 + T2;
  }
}
