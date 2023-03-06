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

extern int FXP_frac_bits;
extern int FXP_frac_bits_m1;
extern int FXP_frac_mask;
extern int FXP_frac_max;
extern int FXP_frac_max_p1;
extern int FXP_whole_bits;
extern int FXP_whole_bits_m1;
extern int FXP_whole_max;
extern int FXP_whole_min;
extern int FXP_whole_min_m1;
extern int fxp_frac_max_dec;
extern int fxp_frac_max_dec_p1;

extern float FXP_max_f;
extern float FXP_min_f;
extern double FXP_max_d;
extern double FXP_min_d;
extern long double FXP_max_ld;
extern long double FXP_min_ld;

//extern int FXP_BKM_L_CLZSHIFT;
extern unsigned long FXP_BKM_ONE_L;
extern unsigned long FXP_BKM_HALF_L;
extern int FXP_shifted_e;
extern int FXP_shifted_pi;
extern int FXP_shifted_ln_2;
extern int FXP_shifted_lg10_2;
extern int FXP_half;
extern int FXP_one;
extern int FXP_two;
extern int FXP_lg2_maxloops;
extern int FXP_lg2_l_mshift;
extern int FXP_lg2_l_mshiftm1;
extern int FXP_pow2_l_xshift;
extern int FXP_pow2_l_xshiftm1;
extern unsigned long FXP_max_lshifted;

extern const int FXP_WORD_BITS;
extern const int FXP_WORD_BITS_M1;
extern const unsigned int FXP_LWORD_MASK;
extern const unsigned int FXP_RWORD_MASK;