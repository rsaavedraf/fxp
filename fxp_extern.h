/* SPDX-License-Identifier: MIT */
/*
 * fxp_extern.h
 *
 * All external variables defined and used by fxp.c,
 * which the auxiliary / satellite programs also need
 *
 * By Raul Saavedra, Bonn, Germany
 */

extern const int FXP_INT_BITS;
extern const int FXP_INT_BITS_M1;
extern const int FXP_INT_BITS_M1_NEG;
extern const int FXP_FRAC_BITS_DEF;
extern const int FXP_POS_INF;
extern const int FXP_NEG_INF;
extern const int FXP_MAX;
extern const long FXP_MAX_L;
extern const int FXP_MIN;
extern const int FXP_UNDEF;

extern const float FXP_PINF_F;
extern const float FXP_UNDEF_F;
extern const float FXP_NINF_F;
extern const double FXP_PINF_D;
extern const double FXP_UNDEF_D;
extern const double FXP_NINF_D;
extern const long double FXP_PINF_LD;
extern const long double FXP_UNDEF_LD;
extern const long double FXP_NINF_LD;

extern const unsigned int FXP_PI_I32;
extern const unsigned long FXP_PI_I64;
extern const unsigned long FXP_LG10_2_I64;
extern const unsigned long FXP_LN_2_I64;
extern const unsigned long FXP_LG10_2_I64;
extern const unsigned long FXP_LG2_E_I64;
extern const unsigned long FXP_LG2_10_I64;
extern const unsigned int FXP_LG2_E_WBITS;
extern const unsigned int FXP_LG2_10_WBITS;

extern int FXP_frac_bits;
extern int FXP_frac_bits_m1;
extern unsigned int FXP_frac_mask;
extern unsigned int FXP_frac_max;
extern unsigned int FXP_frac_max_p1;
extern int FXP_whole_bits;
extern int FXP_whole_bits_m1;
extern int FXP_whole_bits_m2;
extern int FXP_whole_max;
extern int FXP_whole_min;
extern int FXP_whole_min_m1;
extern long fxp_frac_max_dec;
extern long fxp_frac_max_dec_p1;

extern float FXP_min_fx;
extern float FXP_max_fx;
extern double FXP_min_dx;
extern double FXP_max_dx;
extern long double FXP_min_ldx;
extern long double FXP_min_ld;
extern long double FXP_max_ld;
extern long double FXP_max_ldx;

extern unsigned int FXP_shifted_e;
extern unsigned int FXP_shifted_pi;
extern unsigned int FXP_shifted_ln_2;
extern unsigned int FXP_shifted_lg10_2;
extern int FXP_half;
extern int FXP_one;
extern int FXP_two;
extern int FXP_almost1;
extern int FXP_lg2_maxloops;
extern int FXP_lg2_l_mshift;
extern int FXP_lg2_l_mshift_m1;
extern int FXP_lg2_l_mshift_p1;
extern unsigned long FXP_max_lshifted;

extern const int FXP_WORD_BITS;
extern const int FXP_WORD_BITS_M1;
extern const unsigned int FXP_LWORD_MASK;
extern const unsigned int FXP_RWORD_MASK;
extern const unsigned long FXP_LINT_MASK;
extern const unsigned long FXP_RINT_MASK;

extern const int FXP_LONG_BITS;
extern const int FXP_LONG_BITS_M1;
extern const unsigned long FXP_BKM_A_ONE_L;
extern const unsigned long FXP_BKM_A_POINT5_L;
extern const unsigned long FXP_BKM_X_ONE_L;

extern const unsigned int SFXP_MAX_WBITS;

extern const int FXP_LOGX_LOOPS;
extern const int FXP_POWX_LOOPS;
extern const int FXP_SQRT_LOOPS;
extern const int FXP_POWXY_LOOPS;

