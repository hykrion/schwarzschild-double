#ifndef UI_H_INCLUDED
#define UI_H_INCLUDED

/* Number of nodes. We cannot go beyond 32768 (but this's not an issue at all) */
#define NODES 1001

void ui_init_only(void);

double ui_get_a(void);
double ui_get_b(void);
double ui_get_h_forward(void);
double ui_get_h_backguard(void);
int ui_get_nodes(void);
int ui_get_n(void);
double ui_get_wMin(void);
double ui_get_wMax(void);
double ui_get_hW(void);
int ui_get_nW(void);
double ui_get_rS(void);
double ui_get_l(void);
void ui_set_l(double l);
double ui_get_nL(void);

#endif // UI_H_INCLUDED
