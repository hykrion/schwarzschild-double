#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "header/ui.h"
#include "header/phi.h"
#include "header/txt_data.h"
#include "header/tortoise.h"

/**
  @brief  Calculate phi solving the ODE y'(x)
*/
void
calculate_phi_yx_xy(void)
{
  double l = ui_get_l();
  tortoise_calculate_yx_xy();
  /* Do the same calculation for a serie of different angular momentums */
  int nL = ui_get_nL();
  int i;
  for(i = 0; i < nL; i++)
  {
    phi_wave_yx_xy(l);
    phi_calculate_sigma_l(l);
    l += 1.0;
    /*
    The angular momentum comes from the GUI... so we have to save the modification
    */
    ui_set_l(l);
  }
}
/* -----------------------------------------------------------------------------
----------------------------------------------------------------------------- */
int
main()
{
  /* Read parameters from GUI */
  ui_init_only();

  tortoise_init();

  /* Calculate phi solving y'(x) so we can set an IC to solve x'(y) */
  calculate_phi_yx_xy();

  /* Get the potential values to be able to plot it */
  phi_calculate_v();

  /* phi Re, Im & potential */
  txt_data_wave();

  /* Tortoise coordinates */
  txt_data_turtle();

  /* Coefficient values R, T y R + T */
  txt_data_rt();

  /* sigma_l */
  txt_data_sigma();

  tortoise_destroy();

  return 0;
}
