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
int fxp_sqrt_l(int fxp1);
int fxp_powxy_l(int fxp1, int fxp2);

int fxp_sin_l(int fxp1);
int fxp_cos_l(int fxp1);
int fxp_tan_l(int fxp1);
int fxp_asin_l(int fxp1);
int fxp_acos_l(int fxp1);
int fxp_atan_l(int fxp1);

unsigned long dmul_ulongs(unsigned long x, unsigned long y);
unsigned long dmul_ulong_x_uint(unsigned long x, unsigned int y);

typedef struct super_fxp_l {
        int sign;
        int nwbits;             // <- Number of whole bits
        unsigned long number;   // <- The super fxp number
} super_fxp_l;

super_fxp_l sfxp_l_from_fxp(int fxp1);
int sfxp_l_2_fxp(super_fxp_l x);
