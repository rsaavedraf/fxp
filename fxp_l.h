/* SPDX-License-Identifier: MIT */
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
#include "fxp_extern.h"

int fxp_mul_l(int fxp1, int fxp2);
int fxp_div_l(int fxp1, int fxp2);

int fxp_lg2_l(int fxp1);
int fxp_lg2_mul_l(int fxp1);
int fxp_ln_l(int fxp1);
int fxp_lg10_l(int fxp1);
int fxp_pow2_l(int fxp1);
int fxp_exp_l(int fxp1);

int fxp_pow10_l(int fxp1);
int fxp_powxy_l(int fxp1, int fxp2);
