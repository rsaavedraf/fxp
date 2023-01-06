/*
 * fxp.h
 * An implementation of binary Fixed-Point numbers
 * encoding them into integers
 * By Raul Saavedra, 2022-11-13
 */

#include <limits.h>

#define FXP_INT_BITS (sizeof(int) * 8)

// Number of bits to use for the frac part
#define FXP_FRAC_BITS 12

// Number of bits for the whole part
#define FXP_WHOLE_BITS (FXP_INT_BITS - FXP_FRAC_BITS)

// FXP_FRAC_MAX should correspond to 2^FXP_FRAC_BITS - 1
// Also used as mask for binary frac part of the int
#define FXP_FRAC_MAX ((1 << FXP_FRAC_BITS) - 1)

// Desired max frac decimal value
#define FXP_FRAC_MAX_DEC 999

// Infinity constants and Max and min fxp valid values
#define FXP_POS_INF INT_MAX
#define FXP_NEG_INF (-INT_MAX)
#define FXP_MAX (FXP_POS_INF - 1) // max whole + 1 bit short of max frac
#define FXP_MIN (FXP_NEG_INF + 1) // min whole + 1 bit short of min frac
#define FXP_MAX_L ((long) FXP_MAX)
#define FXP_MAX_LSHIFTED ((FXP_MAX_L) << FXP_FRAC_BITS)

// Saving the true most negative int for "undefined"
#define FXP_UNDEF INT_MIN

// Max and min valid values for the whole part of the fxp's
#define FXP_WHOLE_MAX (INT_MAX >> FXP_FRAC_BITS)
#define FXP_WHOLE_MIN (-FXP_WHOLE_MAX)

int fxp(int whole);
int fxp_bin(int whole, int frac);
int fxp_dec(int whole, int dec_frac);

int fxp_get_whole_part(int fxp);
int fxp_get_bin_frac(int fxp);
int fxp_get_dec_frac(int fxp);

int fxp_unsafe_sum(int fxp1, int fxp2);
int fxp_unsafe_sub(int fxp1, int fxp2);
int fxp_unsafe_mul(int fxp1, int fxp2);
int fxp_unsafe_div(int fxp1, int fxp2);

unsigned int nbits(unsigned int n, unsigned int maxbits);
unsigned int nbits_v0(unsigned int n);

int fxp_sum(int fxp1, int fxp2);
int fxp_sub(int fxp1, int fxp2);
int fxp_mul(int fxp1, int fxp2);
int fxp_div(int fxp1, int fxp2);
