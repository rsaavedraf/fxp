/*
 * fxp_l.h
 *
 * Header file for the alternative (faster) *_l functions for
 * fxp's. Faster because of using longs instead of only
 * exclusively ints.
 *
 * This can be used only when the target system's
 * sizeof(long) >= 2 * syzeof(int)
 *
 * By Raul Saavedra, Bonn, Germany
 */

#include <limits.h>
#include "fxp.h"

extern int FXP_frac_bits;
extern int FXP_whole_min;
extern int FXP_BKM_PREC;
extern int FXP_BKM_PREC_L;
extern int FXP_BKM_L_CLZSHIFT;
extern unsigned long FXP_BKM_ONE_L;
extern unsigned long FXP_BKM_HALF_L;
extern int FXP_shifted_e;
extern int FXP_shifted_pi;
extern int FXP_shifted_ln_2;
extern int FXP_shifted_lg10_2;
extern int FXP_half;
extern int FXP_one;
extern int FXP_two;
extern unsigned long FXP_max_lshifted;
extern int FXP_lg2_l_maxloops;
extern int FXP_lg2_l_shift;

int fxp_mul_l(int fxp1, int fxp2);
int fxp_div_l(int fxp1, int fxp2);
int fxp_lg2_l(int fxp);
int fxp_lg2_mul_l(int fxp);
int fxp_ln_l(int fxp);
int fxp_lg10_l(int fxp);
