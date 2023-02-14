/*
 * fxp.h
 * An implementation of binary Fixed-Point numbers
 * encoding them into integers.
 * By Raul Saavedra, Bonn, Germany
 *
 * v1.0: 2022-11-13
 * v2.0: 2023-01-07: runtime-modifiable whole vs. frac bits
 */

#include <limits.h>

#define FXP_INT_BITS ((int) (sizeof(int) * 8))
#define FXP_INT_BITS_M1 (FXP_INT_BITS - 1)
//#define FXP_INT_BITS_HALF (FXP_INT_BITS / 2)

// Infinity constants and Max and min fxp valid values
#define FXP_POS_INF INT_MAX
#define FXP_NEG_INF (-INT_MAX)
#define FXP_MAX (FXP_POS_INF - 1) // max whole - 1 (1 bit short of max frac)
#define FXP_MIN (FXP_NEG_INF + 1) // min whole + 1 (1 bit short of min frac)
#define FXP_MAX_L ((long) FXP_MAX)

// Saving the true most negative int for "undefined"
#define FXP_UNDEF INT_MIN

// Number of bits to use for the frac part
int fxp_get_frac_bits();
int fxp_set_frac_bits(int vfracbits);

// Number of bits for the whole (and sign) part
int fxp_get_whole_bits();

int fxp_get_frac_mask();
int fxp_get_frac_max();

// Desired max frac decimal value
int fxp_get_frac_max_dec();
int fxp_set_frac_max_dec(int vfracmaxdec);
int fxp_set_auto_frac_max_dec();

// Max and min valid values for the whole part of the fxp's
int fxp_get_whole_max();
int fxp_get_whole_min();

int fxp(int whole);
int fxp_bin(int whole, int frac);
int fxp_dec(int whole, int dec_frac);

int fxp_get_whole_part(int fxp);
int fxp_get_bin_frac(int fxp);
int fxp_get_dec_frac(int fxp);

int fxp_unsafe_add(int fxp1, int fxp2);
int fxp_unsafe_sub(int fxp1, int fxp2);
int fxp_unsafe_mul(int fxp1, int fxp2);
int fxp_unsafe_div(int fxp1, int fxp2);

int fxp_nbits(unsigned int n);

// Safe implementations using only ints
int fxp_add(int fxp1, int fxp2);
int fxp_sub(int fxp1, int fxp2);
int fxp_mul(int fxp1, int fxp2);
int fxp_div(int fxp1, int fxp2);

// Safe implementations using longs
//int fxp_add_l(int fxp1, int fxp2);
//int fxp_sub_l(int fxp1, int fxp2);
int fxp_mul_l(int fxp1, int fxp2);
int fxp_div_l(int fxp1, int fxp2);

// Beyond the basic four arithmetic ops
int fxp_get_e();
int fxp_get_pi();
int fxp_get_ln2();
int fxp_get_lg2e();
int fxp_lg2(int fxp);
int fxp_lg2_l(int fxp);
int fxp_lg2_mul(int fxp);
int fxp_lg2_mul_l(int fxp);
int fxp_ln(int fxp);
int fxp_ln_l(int fxp);

// Still to implement
//int fxp_sqrt(int fxp);
//int fxp_pow(int fxp_base, int fxp_exponent);
