/*
 * fxp_conv.c
 *
 * Provides functions to do conversions both ways
 * between fxp and float, or double, or long double
 *
 */

#include "fxp_conv.h"

// The largest valid negative values of each type are used
// as the special values for 'undef' for the purposes of
// conversion into and from fxp's.
// Negating those largest values here assumes symmetry in the
// floating-point implementations. If that does not hold on
// the target platform, then FXP_UCOEF can simply be
// made slighty < 1.0
#define FXP_UCOEF 1.0
#define FXP_NICOEF (FXP_UCOEF - 0.01)
const float FXP_PINF_F = FXP_MAX;
const float FXP_UNDEF_F = -(FLT_MAX * FXP_UCOEF);
const float FXP_NINF_F = (FXP_UNDEF_F * FXP_NICOEF);
const double FXP_PINF_D = DBL_MAX;
const double FXP_UNDEF_D = -(DBL_MAX * FXP_UCOEF);
const double FXP_NINF_D = (FXP_UNDEF_D * FXP_NICOEF);
const long double FXP_PINF_LD = LDBL_MAX;
const long double FXP_UNDEF_LD = -(LDBL_MAX * FXP_UCOEF);
const long double FXP_NINF_LD = (FXP_UNDEF_LD * FXP_NICOEF);

/*
float fxp_get_undef_f()
{
        return FXP_UNDEF_F;
}

float fxp_get_neg_inf_f()
{
        return FXP_NINF_F;
}

float fxp_get_pos_inf_f()
{
        return FXP_PINF_F;
}

double fxp_get_undef_d()
{
        return FXP_UNDEF_D;
}

double fxp_get_neg_inf_d()
{
        return FXP_NINF_D;
}

double fxp_get_pos_inf_d()
{
        return FXP_PINF_D;
}

long double fxp_get_undef_ld()
{
        return FXP_UNDEF_LD;
}

long double fxp_get_neg_inf_ld()
{
        return FXP_NINF_LD;
}

long double fxp_get_pos_inf_ld()
{
        return FXP_PINF_LD;
}
*/

/*
 * Returns a float value corresponding to an fxp
 */
float fxp2f(int fxp)
{
    if (fxp == FXP_UNDEF) return FXP_UNDEF_F;
    if (fxp == FXP_NEG_INF) return FXP_NINF_F;
    if (fxp == FXP_POS_INF) return FXP_PINF_F;
    unsigned int twopower = 1 << FXP_frac_bits;
    int frac = fxp_get_bin_frac(fxp);
    if (frac < 0) frac = -frac;
    float ffrac = 0.0;
    while (frac > 0) {
        if (frac & 1) ffrac += ((float) 1.0 / twopower);
        frac = frac >> 1;
        twopower = twopower >> 1;
    }
    float wf = (float) fxp_get_whole_part(fxp);
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
    if (x < FXP_min_f) return FXP_NEG_INF;
    if (x > FXP_max_f) return FXP_POS_INF;
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
    if (x < FXP_min_d) return FXP_NEG_INF;
    if (x > FXP_max_d) return FXP_POS_INF;
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
    unsigned int twopower = 1 << FXP_frac_bits;
    int frac = fxp_get_bin_frac(fxp);
    if (frac < 0) frac = -frac;
    long double ldfrac = 0.0;
    while (frac > 0) {
        if (frac & 1) ldfrac += ((long double) 1.0 / twopower);
        frac = frac >> 1;
        twopower = twopower >> 1;
    }
    long double wld = (long double) fxp_get_whole_part(fxp);
    return wld + ((fxp < 0)? -ldfrac: ldfrac);
}

/*
 * Returns the fxp corresponding to a long double value
 */
int ld2fxp(long double x)
{
    if (isnan(x)) return FXP_UNDEF;
    if (isinf(x)) return (signbit(x)? FXP_NEG_INF: FXP_POS_INF);
    if (x <= FXP_UNDEF_LD) return FXP_UNDEF;
    if (x < FXP_min_ld) return FXP_NEG_INF;
    if (x > FXP_max_ld) return FXP_POS_INF;
    long double px = (x < 0)? -x: x;
    long double pwld = truncl(px);
    long double frac = px - pwld;
    unsigned int shift = 1 << FXP_frac_bits;
    int pfrac_i = (int) truncf(frac * shift);
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
