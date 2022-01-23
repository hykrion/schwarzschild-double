#ifndef V_H_INCLUDED
#define V_H_INCLUDED

struct xyParams
{
  double rS;
};

/* Barrera o pozo potencial */
double V(double y);
/* Vector x[y] */
int xyIntegration();


#endif
