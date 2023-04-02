/* SPDX-License-Identifier: MIT */
/*
 * fxp.h
 *
 * Implementation of binary Fixed-Point numbers
 * encoding them into integers.
 *
 * By Raul Saavedra, Bonn, Germany
 */

#include <limits.h>

// Number of bits to use for the whole vs. frac parts
int fxp_get_frac_bits();
int fxp_set_frac_bits(int nfb);

int fxp_get_whole_bits();
int fxp_get_frac_mask();
int fxp_get_frac_max();

// Desired max frac decimal value
int fxp_get_frac_max_dec();
int fxp_set_frac_max_dec(int nfmaxdec);
int fxp_set_auto_frac_max_dec();

// Max and min valid values for the whole part of the fxp's
int fxp_get_whole_max();
int fxp_get_whole_min();

// Max and min fxp values converted to floating-point types
float fxp_get_max_f();
float fxp_get_min_f();
double fxp_get_max_d();
double fxp_get_min_d();
long double fxp_get_max_ld();
long double fxp_get_min_ld();

// Create an fxp
int fxp(int whole);
int fxp_bin(int whole, int frac);
int fxp_dec(int whole, int dec_frac);

// Inspect an fxp
int fxp_get_whole_part(int fxp1);
int fxp_get_bin_frac(int fxp1);
int fxp_get_dec_frac(int fxp1);

// Aux function to get position of most signif bit
int fxp_nbits(unsigned int n);

// Safe implementations of the basic arithmetics
int fxp_add(int fxp1, int fxp2);
int fxp_sub(int fxp1, int fxp2);
int fxp_mul(int fxp1, int fxp2);
unsigned int mul_distrib(unsigned int x, unsigned int y);
int fxp_mul_wrefs(int fxp1, int fxp2, int * rw, int * rf);
int fxp_div(int fxp1, int fxp2);

// Unsafe arithmetic operations (beware)
int fxp_unsafe_add(int fxp1, int fxp2);
int fxp_unsafe_sub(int fxp1, int fxp2);
int fxp_unsafe_mul(int fxp1, int fxp2);
int fxp_unsafe_div(int fxp1, int fxp2);

// Important constants in fxp format
unsigned int fxp_get_e();
unsigned int fxp_get_pi();
unsigned int fxp_get_ln_2();
unsigned int fxp_get_lg10_2();

// Beyond the basic four ops
int fxp_lg2(int fxp1);
int fxp_ln(int fxp1);
int fxp_lg10(int fxp1);
int fxp_pow2(int fxp1);
int fxp_exp(int fxp1);
int fxp_pow10(int fxp_x);

int fxp_powxy(int fxp_x, int fxp_y);
int fxp_sqrt(int fxp1);

// Trigonometric functions

