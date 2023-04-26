/* SPDX-License-Identifier: MIT */
/*
 * fxp_conv.c
 *
 * Provides functions to do conversions both ways
 * between fxp and float, or double, or long double
 */

#include <stdio.h>
#include "fxp_conv.h"

// Preserving these relationship as they stand for fxp's:
// UNDEF < -INF == -(+INF)
// FXP_UCOEF and FXP_NICOEF can always be adjusted
// (ie. FXP_UCOEF made slighty < 1.0,) if there is no
// +/- symmetry in the target floating-point
// implementation/platform
const long double FXP_UCOEF = 1.0;
const long double FXP_NICOEF = FXP_UCOEF - 0.01;
const float FXP_UNDEF_F = -(FLT_MAX * FXP_UCOEF);
const float FXP_NINF_F = FXP_UNDEF_F * FXP_NICOEF;
const float FXP_PINF_F = -FXP_NINF_F;
const double FXP_UNDEF_D = -DBL_MAX * FXP_UCOEF;
const double FXP_NINF_D = FXP_UNDEF_D * FXP_NICOEF;
const double FXP_PINF_D = -FXP_NINF_D;
const long double FXP_UNDEF_LD = -LDBL_MAX * FXP_UCOEF;
const long double FXP_NINF_LD = FXP_UNDEF_LD * FXP_NICOEF;
const long double FXP_PINF_LD = -FXP_NINF_LD;

/*
 * Returns a float value corresponding to an fxp
 */
float fxp2f(int fxp)
{
    if (fxp == FXP_UNDEF) return FXP_UNDEF_F;
    if (fxp == FXP_NEG_INF) return FXP_NINF_F;
    if (fxp == FXP_POS_INF) return FXP_PINF_F;
    // TODO: double-check how this rounding is done
    // An IEEE-754 float uses 24 bits for the
    // mantissa part (first implicit, 23 explicit)
    // So we will round up that last bit in it
    int pfxp = fxp < 0? -fxp: fxp;
    int nbits = fxp_nbits(pfxp);
    int pw = fxp_get_whole_part(fxp);
    int frac = fxp_get_bin_frac(fxp);
    if (frac < 0) frac = -frac;
    unsigned int twopower = 1 << FXP_frac_bits;
    float ffrac = 0.0;
    int rbit = 0;
    while (frac > 0) {
            int bit = frac & 1;
            if (nbits >= 25) {
                    rbit = bit;
            } else {
                    if ((frac & 1) || rbit)
                            ffrac += ((float) 1.0 / twopower);
                    rbit = 0;
            }
            frac = frac >> 1;
            nbits--;
            twopower = twopower >> 1;
    }
    float wf = (float) pw;
    return wf + ((fxp < 0)? -ffrac: ffrac);
}

/*
 * Returns the fxp corresponding to a float value
 */
int f2fxp(float x)
{
    if (isnan(x)) return FXP_UNDEF;
    if (isinf(x)) return (signbit(x)? FXP_NEG_INF: FXP_POS_INF);
    if (x <= FXP_UNDEF_F) return FXP_UNDEF;
    if ((x <= FXP_NINF_F) || (x <= FXP_min_fx)) return FXP_NEG_INF;
    if ((x >= FXP_PINF_F) || (x >= FXP_max_fx)) return FXP_POS_INF;
    float px = (x < 0)? -x: x;
    float pwf = truncf(px);
    float frac = px - pwf;
    unsigned int shift = 1 << FXP_frac_bits;
    int pfrac_i = (int) truncf(frac * shift);
    int pwhole_i = (int) pwf;
    if (pwhole_i > 0)
        if (x > 0.0)
            return fxp_bin(pwhole_i, pfrac_i);
        else
            return fxp_bin(-pwhole_i, pfrac_i);
    else    // pwhole == 0
        if (x >= 0.0)
            return fxp_bin(0, pfrac_i);
        else
            return fxp_bin(0, -pfrac_i);
}

/*
 * Returns the double value corresponding to an fxp
 */
double fxp2d(int fxp)
{
    if (fxp == FXP_UNDEF) return FXP_UNDEF_D;
    if (fxp == FXP_NEG_INF) return FXP_NINF_D;
    if (fxp == FXP_POS_INF) return FXP_PINF_D;
    unsigned int twopower = 1 << FXP_frac_bits;
    int frac = fxp_get_bin_frac(fxp);
    if (frac < 0) frac = -frac;
    double dfrac = 0.0;
    while (frac > 0) {
        if (frac & 1) dfrac += ((double) 1.0 / twopower);
        frac = frac >> 1;
        twopower = twopower >> 1;
    }
    double wd = (double) fxp_get_whole_part(fxp);
    return wd + ((fxp < 0)? -dfrac: dfrac);
}

/*
 * Returns the fxp corresponding to a double value
 */
int d2fxp(double x)
{
    if (isnan(x)) return FXP_UNDEF;
    if (isinf(x)) return (signbit(x)? FXP_NEG_INF: FXP_POS_INF);
    if (x <= FXP_UNDEF_D) return FXP_UNDEF;
    if (x <= FXP_min_dx) return FXP_NEG_INF;
    if (x >= FXP_max_dx) return FXP_POS_INF;
    double px = (x < 0)? -x: x;
    double pw = trunc(px);
    double frac = px - pw;
    unsigned int shift = 1 << FXP_frac_bits;
    int pfrac_i = (int) trunc(frac * shift);
    int pwhole_i = (int) pw;
    if (pwhole_i > 0)
        if (x > 0.0)
            return fxp_bin(pwhole_i, pfrac_i);
        else
            return fxp_bin(-pwhole_i, pfrac_i);
    else    // pwhole == 0
        if (x >= 0.0)
            return fxp_bin(0, pfrac_i);
        else
            return fxp_bin(0, -pfrac_i);
}

/*
 * Returns a long double value corresponding to a given fxp
 */
long double fxp2ld(int fxp)
{
    if (fxp == FXP_UNDEF) return FXP_UNDEF_LD;
    if (fxp == FXP_NEG_INF) return FXP_NINF_LD;
    if (fxp == FXP_POS_INF) return FXP_PINF_LD;
    //printf("\nfxp2ld: fxp is %X\n", fxp);
    unsigned int twopower = 1u << FXP_frac_bits;
    int frac = fxp_get_bin_frac(fxp);
    if (frac < 0) frac = -frac;
    long double ldfrac = 0.0L;
    while (frac > 0) {
        if (frac & 1u) ldfrac += ((long double) 1.0L) / twopower;
        frac = frac >> 1;
        twopower = twopower >> 1;
    }
    long double wld = (long double) fxp_get_whole_part(fxp);
    long double num = wld + ((fxp < 0)? -ldfrac: ldfrac);
    //printf("fxp2ld: ld num is %34.30Le\n", num);
    return num;
}

/*
 * Returns the fxp corresponding to a long double value
 */
int ld2fxp(long double x)
{
    if (isnan(x)) return FXP_UNDEF;
    if (isinf(x)) return (signbit(x)? FXP_NEG_INF: FXP_POS_INF);
    if (x <= FXP_UNDEF_LD) return FXP_UNDEF;
    if (x <= FXP_min_ldx) return FXP_NEG_INF;
    if (x >= FXP_max_ldx) return FXP_POS_INF;
    long double px = (x < 0)? -x: x;
    long double pwld = truncl(px);
    long double frac = px - pwld;
    unsigned int shift = 1 << FXP_frac_bits;
    int pfrac_i = (int) truncl(frac * shift);
    int pwhole_i = (int) pwld;
    if (pwhole_i > 0)
        if (x > 0.0)
            return fxp_bin(pwhole_i, pfrac_i);
        else
            return fxp_bin(-pwhole_i, pfrac_i);
    else    // pwhole == 0
        if (x >= 0.0)
            return fxp_bin(0, pfrac_i);
        else
            return fxp_bin(0, -pfrac_i);
}
