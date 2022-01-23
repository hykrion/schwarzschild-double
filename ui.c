#include "header/ui.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
NOTE  The number and the order of the parameters have to be the same as the
      declared in 'cmd_visualizar' of the file 'gui.tcl'
*/
/* Number of parameters + 1 */
#define UI_PARAM_SIZE 8
/* Max length of the name of parameters */
#define UI_PARAM_LEN 16
/* This order have to be the same as the gui.tcl file */
#define UI_PARAM_A 0
#define UI_PARAM_B 1
#define UI_PARAM_L 2
#define UI_PARAM_RS 3
#define UI_PARAM_WMIN 4
#define UI_PARAM_WMAX 5
#define UI_PARAM_NW 6
#define UI_PARAM_NL 7

static double m_uiParams[UI_PARAM_SIZE];

/**
  @brief  Read the different parameters that we can change from the GUI.
          For simplicity these parameters are in text mode in a file, parameters.txt,
          that will have one parameter for each line.
  */
void
ui_init_only(void)
{
  char str[UI_PARAM_SIZE][UI_PARAM_LEN];

  int i;
  for(i = 0; i < UI_PARAM_SIZE; i++)
    memset(str[i], '\0', UI_PARAM_LEN);

  /* Read the parameters from the file */
  FILE *fp;
  fp = fopen("parameters.txt", "r");

  i = 0;
  while(!feof(fp))
  {
    fgets(str[i], UI_PARAM_LEN - 1, fp);
    i++;
  }

  fclose(fp);

  i = 0;
  while(i < UI_PARAM_SIZE)
  {
    m_uiParams[i] = atof(str[i]);
    i++;
  }
}

double
ui_get_a(void)
{
  return m_uiParams[UI_PARAM_A];
}

double
ui_get_b(void)
{
  return m_uiParams[UI_PARAM_B];
}

double
ui_get_h_forward(void)
{
  double a = ui_get_a();
  double b = ui_get_b();
  double n = ui_get_n();

  return (b - a) / n;
}

double
ui_get_h_backguard(void)
{
  double a = ui_get_a();
  double b = ui_get_b();
  double n = ui_get_n();

  return (a - b) / n;
}

int
ui_get_nodes(void)
{
  return NODES;
}

int
ui_get_n(void)
{
  return NODES - 1;
}

double
ui_get_wMin(void)
{
  return m_uiParams[UI_PARAM_WMIN];
}

double
ui_get_wMax(void)
{
  return m_uiParams[UI_PARAM_WMAX];
}

double
ui_get_hW(void)
{
  double wm = ui_get_wMin();
  double wM = ui_get_wMax();
  double nW = ui_get_nW();

  return (wM - wm) / nW;
}

int
ui_get_nW(void)
{
  return m_uiParams[UI_PARAM_NW];
}

double
ui_get_rS(void)
{
  return m_uiParams[UI_PARAM_RS];
}

double
ui_get_l(void)
{
  return m_uiParams[UI_PARAM_L];
}

double
ui_get_nL(void)
{
  return m_uiParams[UI_PARAM_NL];
}

void ui_set_l(double l)
{
  m_uiParams[UI_PARAM_L] = l;
}
