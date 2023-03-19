
/* SPDX-License-Identifier: MIT */
/*
 * fxp.c
 *
 * An implementation of Binary Fixed Point numbers
 * encoding them into integers, together with their
 * arithmetic operations +, -, *, and /, and more.
 *
 * By Raul Saavedra, Bonn, Germany
 *
 */

#include "fxp.h"
#include <stdio.h>
#include <stdlib.h>
//#include <assert.h>

//Used while testing and debugging trying to optimize division
//#include "fxp_aux.h"
//#define VERBOSE 1

const int FXP_INT_BITS = ((int) sizeof(int)) * 8;
const int FXP_INT_BITS_M1 = FXP_INT_BITS - 1;
const int FXP_INT_BITS_M1_NEG = -FXP_INT_BITS_M1;
const int FXP_INT_BITS_M2 = FXP_INT_BITS - 2;
const int FXP_POS_INF = INT_MAX;
const int FXP_NEG_INF = -INT_MAX;
const int FXP_MAX = FXP_POS_INF - 1;
const long FXP_MAX_L = (long) FXP_MAX;
const int FXP_MIN = FXP_NEG_INF + 1;
// Saving the true most negative int for "undefined"
const int FXP_UNDEF = INT_MIN;
const int FXP_FRAC_BITS_MIN = 4;
const int FXP_FRAC_BITS_DEF = 16;

const int FXP_WORD_BITS = FXP_INT_BITS >> 1;
const int FXP_WORD_BITS_M1 = FXP_WORD_BITS - 1;
const unsigned int FXP_RWORD_MASK = ((1 << FXP_WORD_BITS) - 1);
const unsigned int FXP_LWORD_MASK = FXP_RWORD_MASK \
                                    << FXP_WORD_BITS;
const unsigned long FXP_RINT_MASK = ((1l << FXP_INT_BITS) - 1);
const unsigned long FXP_LINT_MASK = FXP_RINT_MASK \
                                    << FXP_INT_BITS;
const int FXP_LONG_BITS = ((int) sizeof(long)) * 8;
const int FXP_LONG_BITS_M1 = FXP_LONG_BITS - 1;

// Allowing for no whole part, so other than the sign bit,
// all bits used for fraction part. This use case represents
// the range [-0.999..., 0.999...], or equivalently: (-1, 1)
// (e.g. actual whole values -1 and 1 outside the valid range)
const int FXP_FRAC_BITS_MAX = FXP_INT_BITS_M1;
const int FXP_FRAC_MAX_DEC = 9999999;

// Default number of bits to use for the frac part.
int FXP_frac_bits = FXP_FRAC_BITS_DEF;
// FXP_frac_bits can be changed dynamically calling
// fxp_set_frac_bits(), but it remains unchanged until calling
// that function again.
// All such static variables which will only change
// if modifying the number of frac bits are named here
// starting with "FXP_" (uppercase,) but the rest of the
// variable name in lowercase

int FXP_frac_bits_m1 = FXP_FRAC_BITS_DEF - 1;
static int FXP_frac_bits_mod2 = FXP_FRAC_BITS_DEF % 2;

// For improved-precision version of fxp_mul
static int FXP_frac_mshift = FXP_FRAC_BITS_DEF / 2;
static unsigned int FXP_frac_maskl = 4032;
static unsigned int FXP_frac_maskr = 63;

// Default number of bits for the whole part (includes sign bit)
int FXP_whole_bits = FXP_INT_BITS - FXP_FRAC_BITS_DEF;
int FXP_whole_bits_m1 = FXP_INT_BITS_M1 - FXP_FRAC_BITS_DEF;
int FXP_whole_bits_m2 = FXP_INT_BITS_M1 - FXP_FRAC_BITS_DEF - 1;

// FXP_FRAC_MASK should correspond to 2^FXP_FRAC_BITS - 1
unsigned int FXP_frac_mask = ((1u << FXP_FRAC_BITS_DEF) - 1);
int FXP_frac_max = ((1u << FXP_FRAC_BITS_DEF) - 1);
long FXP_frac_max_p1 = (1l << FXP_FRAC_BITS_DEF);

// Default max and min valid values for the whole part of the fxp's
int FXP_whole_max = FXP_MAX >> FXP_FRAC_BITS_DEF;
int FXP_whole_min = -(FXP_MAX >> FXP_FRAC_BITS_DEF);
int FXP_whole_min_m1 = -(FXP_MAX >> FXP_FRAC_BITS_DEF) - 1;

// Default max and min valid values for the conversion types
float FXP_max_f = ((float) (FXP_MAX >> FXP_FRAC_BITS_DEF)) \
                        + (((float) ((1 << FXP_FRAC_BITS_DEF) - 1)) \
                            / (1u << FXP_FRAC_BITS_DEF));
float FXP_min_f = -((float) (FXP_MAX >> FXP_FRAC_BITS_DEF)) \
                        - (((float) ((1 << FXP_FRAC_BITS_DEF) - 1)) \
                            / (1u << FXP_FRAC_BITS_DEF));

double FXP_max_d = ((double) (FXP_MAX >> FXP_FRAC_BITS_DEF)) \
                        + (((double) ((1 << FXP_FRAC_BITS_DEF) - 1)) \
                            / (1ul << FXP_FRAC_BITS_DEF));
double FXP_min_d = -((double) (FXP_MAX >> FXP_FRAC_BITS_DEF)) \
                        - (((double) ((1 << FXP_FRAC_BITS_DEF) - 1)) \
                            / (1ul << FXP_FRAC_BITS_DEF));

long double FXP_max_ld = ((long double) \
                    (FXP_MAX >> FXP_FRAC_BITS_DEF)) \
                    + (((long double) ((1 << FXP_FRAC_BITS_DEF) - 1)) \
                    / (1uL << FXP_FRAC_BITS_DEF));

long double FXP_min_ld = -((long double) \
                    (FXP_MAX >> FXP_FRAC_BITS_DEF)) \
                    - (((long double) ((1 << FXP_FRAC_BITS_DEF) - 1)) \
                        / (1uL << FXP_FRAC_BITS_DEF));

// transcendental constants (2 whole, 30 and 62 frac bits)
// e = 2.718281828...
const unsigned int FXP_E_I32 = 0xADF85459u;
const unsigned long FXP_E_I64 = 0xADF85458A2BB4A9Bul;

// pi = 3.14159265...
const unsigned int FXP_PI_I32 = 0xC90FDAA2u;
const unsigned long FXP_PI_I64 = 0xC90FDAA22168C235ul;

// transcendental constants (0 whole, 32 and 64 frac bits)
// ln(2) = 0.69314718...
const unsigned int FXP_LN_2_I32 = 0xB17217F8u;
const unsigned long FXP_LN_2_I64 = 0xB17217F7D1CF7C72ul;

// lg10(2) = 0.30102999...
const unsigned int FXP_LG10_2_I32 = 0x4D104D42u;
const unsigned long FXP_LG10_2_I64 = 0x4D104D427DE7FD01ul;

static const int FXP_DEF_SHIFT = FXP_INT_BITS - FXP_FRAC_BITS_DEF;
unsigned int FXP_shifted_e = (FXP_E_I32 >> (FXP_DEF_SHIFT - 2));
unsigned int FXP_shifted_pi = (FXP_PI_I32 >> (FXP_DEF_SHIFT - 2));
unsigned int FXP_shifted_ln_2 = (FXP_LN_2_I32 >> FXP_DEF_SHIFT);
unsigned int FXP_shifted_lg10_2 = (FXP_LG10_2_I32 >> FXP_DEF_SHIFT);

// Auxiliary variables used in the lg2 implementations
int FXP_half = 1 << (FXP_FRAC_BITS_DEF - 1);
int FXP_two = 1 << (FXP_FRAC_BITS_DEF + 1);
//int FXP_lg2_maxloops = FXP_FRAC_BITS_DEF + 1;
unsigned int FXP_one = 1u << FXP_FRAC_BITS_DEF;
unsigned long FXP_one_l =  1ul << FXP_FRAC_BITS_DEF;

// Default desired max frac decimal value
// (can be changed dynamically calling set_[auto_]frac_max_dec
int fxp_frac_max_dec = 9999;
int fxp_frac_max_dec_p1 = 10000;

// For the BKM lg2 and pow2 calculations.
// Values in this array are: a[k] = log2(1 + 1/2^k) represented as
// unsigned long fxp's (8 bytes,) with 31 frac bits
static const unsigned int FXP_BKM_LOGS[] = {
        0x80000000, 0x4AE00D1C, 0x2934F097, 0x15C01A39,
        0xB31FB7D, 0x5AEB4DD, 0x2DCF2D0, 0x16FE50B,
        0xB84E23, 0x5C3E0F, 0x2E24CA, 0x1713D6,
        0xB8A47, 0x5C53A, 0x2E2A3, 0x17153,
        0xB8A9, 0x5C55, 0x2E2A, 0x1715,
        0xB8A, 0x5C5, 0x2E2, 0x171,
        0xB8, 0x5C, 0x2E, 0x17,
        0xB, 0x5, 0x2, 0x1,
        0x0, 0x0, 0x0, 0x0
        //0x0, 0x0, 0x0, 0x0,
        //0x0, 0x0, 0x0, 0x0,
        //0x0, 0x0, 0x0, 0x0,
        //0x0, 0x0, 0x0, 0x0,
        //0x0, 0x0, 0x0, 0x0,
        //0x0, 0x0, 0x0, 0x0,
        //0x0, 0x0, 0x0, 0x0,
        //0x0
};

// This array contains the continuation of the frac values
// in the previous array
static const unsigned int FXP_BKM_LOGS_XTRA[] = {
        0x0, 0xFDEB4400, 0x9A371600, 0xFBD68800,
        0x64898B00, 0x63BF61C0, 0xB85A4540, 0x6EF08510,
        0x6BD563B8, 0xFC29D594, 0x6E87E8A8, 0x2F7957C3,
        0x6150DFE4, 0xC47E94D8, 0x2762FA6B, 0x5002E4A,
        0xDED47C11, 0x67F6E58, 0x89050622, 0x45F3D72B,
        0xA35640A7, 0x51C23599, 0xA8E6E01E, 0x5474E163,
        0xAA3ACD06, 0x551D7D98, 0x2A8EC491, 0x154763BA,
        0x8AA3B239, 0xC551D933, 0xE2A8EC9F, 0x71547651,
        0xB8AA3B28, 0x5C551D94, 0x2E2A8ECA, 0x17154765
        //0xB8AA3B2, 0x5C551D9, 0x2E2A8EC, 0x1715476,
        //0xB8AA3B, 0x5C551D, 0x2E2A8E, 0x171547,
        //0xB8AA3, 0x5C551, 0x2E2A8, 0x17154,
        //0xB8AA, 0x5C55, 0x2E2A, 0x1715,
        //0xB8A, 0x5C5, 0x2E2, 0x171,
        //0xB8, 0x5C, 0x2E, 0x17,
        //0xB, 0x5, 0x2, 0x1,
        //0x0
};

// BKM one is a 1 unsigned fxp with 2 whole bits
static const unsigned int FXP_BKM_ONE = 1u << (FXP_INT_BITS - 2);
static const unsigned int FXP_BKM_HALF = FXP_BKM_ONE >> 1;
static const unsigned int FXP_BKM_MAXU = ~(0u);

// Auxiliary variables for the implementations that use longs (fxp_l.c)
const unsigned long FXP_BKM_ONE_L = \
                    ((unsigned long) FXP_BKM_ONE) << FXP_INT_BITS;
const unsigned long FXP_BKM_HALF_L = (FXP_BKM_ONE_L >> 1);

unsigned long FXP_max_lshifted = (FXP_MAX_L << FXP_FRAC_BITS_DEF) \
                    | (((1 << FXP_FRAC_BITS_DEF) - 1) / 2);
int FXP_lg2_l_mshift = 2 * FXP_INT_BITS - 1 - FXP_FRAC_BITS_DEF;
int FXP_lg2_l_mshift_m1 = 2 * FXP_INT_BITS - 2 - FXP_FRAC_BITS_DEF;

/*
 * Given an fxp with x number of frac bits, returns
 * the rounded representation using y frac bits.
 * Used internally to adjust the unsigned
 * transcendental constants when changing the number
 * of frac bits to use.
 */
static unsigned int fxp_rshift_tconst(unsigned int fxp, int x, int y)
{
        int shift = x - y;
        if (shift <= 0) return (unsigned int) FXP_POS_INF;
        unsigned int rbit = (fxp >> (shift - 1)) & 1u;
        return (fxp >> shift) + rbit;
}

/*
 * Dynamic/runtime setting of the bits to use for the frac part.
 * Restricting the usable range of frac bits from FXP_FRAC_BITS_MIN
 * up to FXP_FRAC_BITS_MAX, so that there will always be at least
 * some bits for the fraction part, and at least 1 bit for the
 * sign. This eliminates the need to handle additional special
 * cases
 */
int fxp_set_frac_bits(int nfracbits)
{
        FXP_frac_bits = (nfracbits < FXP_FRAC_BITS_MIN? \
                            FXP_FRAC_BITS_MIN:
                            (nfracbits > FXP_FRAC_BITS_MAX?
                                FXP_FRAC_BITS_MAX: nfracbits));
        FXP_frac_bits_mod2 = FXP_frac_bits % 2;

        FXP_whole_bits = FXP_INT_BITS - FXP_frac_bits;
        FXP_whole_bits_m1 = FXP_whole_bits - 1;
        FXP_whole_bits_m2 = FXP_whole_bits - 2;

        // fxp_frac_mask should correspond to 2^FXP_FRAC_BITS - 1
        FXP_frac_mask = (1u << FXP_frac_bits) - 1;

        // Variables used in fxp_mul to process the
        // full multiplication of the frac parts avoiding
        // precision loss
        if (FXP_frac_bits == FXP_INT_BITS_M1) {
                FXP_frac_mshift = (FXP_INT_BITS_M1 - 1) / 2;
        } else {
                FXP_frac_mshift = (FXP_frac_bits / 2) \
                                    + (FXP_frac_bits % 2);
        }
        FXP_frac_maskr = (1u << FXP_frac_mshift) - 1;
        FXP_frac_maskl = FXP_frac_maskr << FXP_frac_mshift;

        // When using all bits for frac (except for the sign bit),
        // then our max valid frac cannot be equal to frac_mask
        // (because in that case that value is already the largest
        // positive integer == POS_INF), so we must substract one
        // from the frac mask to get the largest valid frac value
        FXP_frac_max = FXP_frac_mask - (FXP_whole_bits == 1? 1: 0);
        FXP_frac_max_p1 = FXP_frac_mask + 1;

        // Auxiliary variables for fxp_l.c implementations
        FXP_max_lshifted = ((FXP_MAX_L) << FXP_frac_bits) \
                                | (FXP_frac_mask / 2);

        // Max and min valid values for the whole part of the fxp's
        FXP_whole_max = FXP_MAX >> FXP_frac_bits;
        FXP_whole_min = (-FXP_whole_max);
        FXP_whole_min_m1 = FXP_whole_min - 1;

        // Max and min floating-point conversion values
        FXP_max_f = ((float) FXP_whole_max) \
                        + ((float) FXP_frac_max) \
                            / ((float) (((unsigned int) 1) \
                                << FXP_frac_bits));
        FXP_max_d = ((double) FXP_whole_max) \
                        + ((double) FXP_frac_max) \
                            / ((double) (((unsigned int) 1) \
                                << FXP_frac_bits));
        FXP_max_ld = (long double) FXP_max_d;
        FXP_min_f = -FXP_max_f;
        FXP_min_d = -FXP_max_d;
        FXP_min_ld = -FXP_max_ld;

        // Adjust precision of e, pi, etc. to the frac bits in use
        FXP_shifted_e = fxp_rshift_tconst(FXP_E_I32, \
                        FXP_INT_BITS_M2, FXP_frac_bits);
        FXP_shifted_pi = fxp_rshift_tconst(FXP_PI_I32, \
                        FXP_INT_BITS_M2, FXP_frac_bits);
        FXP_shifted_ln_2 = fxp_rshift_tconst(FXP_LN_2_I32, \
                        FXP_INT_BITS, FXP_frac_bits);
        FXP_shifted_lg10_2 = fxp_rshift_tconst(FXP_LG10_2_I32, \
                        FXP_INT_BITS, FXP_frac_bits);

        // Auxiliary variables used for lg2 calculations
        FXP_one = 1 << FXP_frac_bits;
        FXP_one_l = 1l << FXP_frac_bits;
        FXP_half = FXP_one >> 1;
        FXP_two = FXP_one << 1;
        //FXP_lg2_maxloops = FXP_frac_bits + 1;
        FXP_lg2_l_mshift = 2 * FXP_INT_BITS - 1 - FXP_frac_bits;
        FXP_lg2_l_mshift_m1 = FXP_lg2_l_mshift - 1;

        return FXP_frac_bits;
}

/*
 * Automatically set the fractional max decimal to use.
 * The number of nines in fxp_frac_max_dec will be
 * floor(fxp_frac_bits / 4), e.g. one nine for every 4 bits
 */
int fxp_set_auto_frac_max_dec()
{
        int nnines = FXP_frac_bits / 4;
        fxp_frac_max_dec = 9;
        for (int i=1; i < nnines; i++) {
            fxp_frac_max_dec = (fxp_frac_max_dec * 10) + 9;
        }
        fxp_frac_max_dec_p1 = fxp_frac_max_dec + 1;
        return fxp_frac_max_dec;
}

/*
 * Manually set the fractional max decimal to use.
 */
int fxp_set_frac_max_dec(int vfracmaxdec)
{
    fxp_frac_max_dec = (vfracmaxdec < 9? 9:
            (vfracmaxdec > FXP_FRAC_MAX_DEC?
                    FXP_FRAC_MAX_DEC: vfracmaxdec));
    fxp_frac_max_dec_p1 = fxp_frac_max_dec + 1;
    return fxp_frac_max_dec;
}

/*
 * Return the max decimal value 99...9 to be used as
 * max posible decimal fraction value
 */
inline int fxp_get_frac_max_dec()
{
    return fxp_frac_max_dec;
}

/*
 * Create an fxp number given its whole part only
 */
inline int fxp(int whole)
{
        return fxp_bin(whole, 0);
}

/*
 * Create an fxp number given its whole and (binary) frac parts.
 */
int fxp_bin(int whole, int bin_frac)
{
        if ((whole == FXP_UNDEF) || (bin_frac == FXP_UNDEF))
                return FXP_UNDEF;
        if (whole > FXP_whole_max)
                return FXP_POS_INF;
        if (whole < FXP_whole_min)
                return FXP_NEG_INF;
        if (bin_frac > FXP_frac_max)
                bin_frac = FXP_frac_max;
        else if (bin_frac < -FXP_frac_max)
                bin_frac = -FXP_frac_max;
        int sign = 1;
        if ((whole == 0) && (bin_frac < 0)) {
                // Special case for negative numbers when whole part is
                // zero, then the fxp gets sign from the frac
                sign = -1;
                bin_frac = -bin_frac;
        } else {
                // All other cases: fxp sign == sign of whole part
                if (whole < 0) {
                        sign = -1;
                        whole = -whole;
                }
                if (bin_frac < 0) {
                        bin_frac = -bin_frac;
                }
        }
        int positive_fxp = (whole << FXP_frac_bits) | bin_frac;
        return (sign == 1)? positive_fxp: -positive_fxp;
}

/*
 * Create an fxp number given its whole and (decimal) frac parts.
 * dec_frac is expected to be a decimal number between 0 and
 * fxp_frac_max_dec, e.g. between 0 and 999
 * Usage examples:
 *     For fxp=16.001, you would invoke: fxp_dec(16, 1)
 *     For 20.09: fxp_dec(24, 90)
 *     For 24.5:  fxp_dec(24, 500)
 * Note last example frac is not 5 but 500. (5 would correspond to
 * 24.005) This decimal value gets scaled into the binary range
 * available for frac. * For negative numbers with whole part=0, the
 * frac value must be negative. * If frac is too large it will get
 * truncated to its most * significant decimal digits until under the
 * value of FXP_FRAC_MAX_DEC, * e.g. for a max of 999, a frac=987654
 * will be truncated to 987
 * (Truncating and not rounding, because the latter would require
 * changing the whole part in some border cases)
 */
int fxp_dec(int whole, int dec_frac)
{
        if ((whole == FXP_UNDEF) || (dec_frac == FXP_UNDEF))
                return FXP_UNDEF;
        if (whole > FXP_whole_max)
                return FXP_POS_INF;
        if (whole < FXP_whole_min)
                return FXP_NEG_INF;
        int frac_sign = 1;
        if (dec_frac < 0) {
                frac_sign = -1;
                dec_frac = -dec_frac;
        }
        int trunc_frac = dec_frac;
        while (trunc_frac > fxp_frac_max_dec) {
                trunc_frac = trunc_frac / 10;
                //printf("   fxp_from_dec_frac: frac trimmed to: %d\n", \
                //    trunc_frac);
        }
        // Watch out this conversion itself can overflow when frac_bits
        // is large, and the frac_max_dec value is also large.
        // Using longs here because of this
        //int bin_frac = (trunc_frac * fxp_frac_max) / fxp_frac_max_dec ;
        int bin_frac = (int) (((long) trunc_frac * \
                                (long) FXP_frac_max) / \
                                    (long) (fxp_frac_max_dec));
        return fxp_bin(whole, ((frac_sign == 1)? bin_frac: -bin_frac));
}

inline int fxp_get_whole_part(int fxp)
{
        if (fxp < 0)
                if (fxp <= FXP_NEG_INF)
                        // Technically, the whole part of -INF would
                        // still be -INF, and both the whole and frac
                        // parts of UNDEF would still be UNDEF. But here
                        // returning just what our whole-part values for
                        // -INF (and also UNDEF) actually correspond to.
                        // However, see also the comment in
                        // fxp_get_bin_frac()
                        return -FXP_whole_max;
                else
                        return -((-fxp) >> FXP_frac_bits);
        else
                return (fxp >> FXP_frac_bits);
}

/*
 * Get the frac part directly (binary)
 */
inline int fxp_get_bin_frac(int fxp)
{
        if (fxp < 0)
                if (fxp <= FXP_NEG_INF)
                        // In the whole-part function above we are
                        // returning the same whole parts for both
                        // -INF and UNDEF. Here for -INF we return the
                        // value used, but for UNDEF something
                        // different and odd in order to enable
                        // differenciating -INF from UNDEF when
                        // checking their whole and frac constituents.
                        // It must be something that no other fxp
                        // could ever have, not even the +/-INF
                        // values. That can be some (alleged) frac
                        // bits that are actually invalid, e.g.
                        // outside the range of valid frac bits, like
                        // -frac_max - 1. Notice that this value
                        // technically "overflows" the frac bits, but
                        // in a way makes sense to mark and
                        // differenciate precisely only UNDEF this way.
                        // Conveniently, when the # of frac bits is
                        // 31, this scheme returns the full FXP_UNDEF
                        // as the frac part.
                        return -FXP_frac_mask \
                                - (fxp == FXP_UNDEF? 1: 0);
                else
                        return -((-fxp) & FXP_frac_mask);
        else
                return (fxp & FXP_frac_mask);
}

/*
 * Get the frac part as decimal between 0 and fxp_frac_max_dec,
 * e.g. 0 .. 9999
 */
int fxp_get_dec_frac(int fxp)
{
        // Watch out the bin to dec conversion itself can overflow when
        // the chosen frac_bits is large, and the frac_max_dec
        // value is also large. Using longs here because of this
        long num, ldiv;
        long positive_frac;
        if (fxp < 0) {
                positive_frac = (-fxp) & FXP_frac_mask;
                num = -(((long) positive_frac) \
                            * fxp_frac_max_dec_p1);
        } else {
                positive_frac = fxp & FXP_frac_mask;
                num = ((long) positive_frac) \
                            * fxp_frac_max_dec_p1;
        }
        ldiv = num / FXP_frac_max_p1;
        return (int) ldiv;
}

/*
 * Simpler unsafe implementations of the arithmetic operations
 * add, sub, mul, and div.
 */
int fxp_unsafe_add(int fxp1, int fxp2)
{
        return fxp1 + fxp2;
}

int fxp_unsafe_sub(int fxp1, int fxp2)
{
        return fxp1 - fxp2;
}

int fxp_unsafe_mul(int fxp1, int fxp2)
{
        long p = ((long) fxp1) * fxp2;
        return (int) (p >> FXP_frac_bits);
}

int fxp_unsafe_div(int fxp1, int fxp2)
{
        long n1 = ((long) fxp1) << FXP_frac_bits;
        long div = n1 / fxp2;
        return ((int) div);
}

/*
 * Default safe implementations of fxp1 + fxp2
 * using only ints.
 * Works for systems in which sizeof(long) is not
 * larger than sizeof(int).
 */
int fxp_add(int fxp1, int fxp2)
{
        // Check for undef or infinity arguments
        if (fxp1 == FXP_UNDEF || fxp2 == FXP_UNDEF)
                return FXP_UNDEF;
        if (fxp1 == FXP_POS_INF || fxp2 == FXP_POS_INF)
                return (fxp1 == FXP_NEG_INF || fxp2 == FXP_NEG_INF)?
                        FXP_UNDEF: FXP_POS_INF;
        if (fxp1 == FXP_NEG_INF || fxp2 == FXP_NEG_INF)
                return FXP_NEG_INF;

        // Rewrite closest to the INT32-C recommendation
        if (((fxp2 > 0) && (fxp1 > (FXP_MAX - fxp2))) || \
            ((fxp2 < 0) && (fxp1 < (FXP_MIN - fxp2))))
                // Return infinity of the appropriate sign
                return (fxp2 > 0? FXP_POS_INF: FXP_NEG_INF);
        // No overflow danger, do sum
        return fxp1 + fxp2;
}

/*
 * Safe fxp1 - fxp2 just using the safe add functions :)
 */
inline int fxp_sub(int fxp1, int fxp2)
{
        // Conservative check to never attempt to change sign
        // to the true most negative int (aka. FXP_UNDEF)
        if (fxp2 == FXP_UNDEF) return FXP_UNDEF;
        return fxp_add(fxp1, -fxp2);
}

/*
 * Return number of bits used by x, that is,
 * most significant bit in x that is a 1
 * (returned value is between 0 and FXP_INT_BITS)
 */
inline int fxp_nbits(unsigned int x)
{
    if (x == 0) return 0;
    // Replacing with gcc builtin function that counts leading zeros
    return FXP_INT_BITS - __builtin_clz(x);
}

/*
 * Multiplies two frac values represented as unsigned ints.
 * Notice that the multiplication of two frac values can never
 * overflow regardless of their representing magnitudes
 * (operands are < 1, so the product is always < 1.)
 * Assumes the arguments are already properly shifted, so that
 * they have the binary point right to the left of the most
 * significant bit in an unsigned int.
 * Uses the distributive multiplication approach:
 * Identifying xa, xb, ya and yb as the inner "words"
 * (e.g. half ints) of the unsigned int arguments, so
 * with WORD_BITS = (FXP_INT_BITS / 2) :
 *      x == (xa << WORD_BITS) | xb
 *      y == (ya << WORD_BITS) | yb
 * Calculates the product as:
 *      x * y = (xa * ya) + ql1 + ql2 + ql3, where
 *      ql1 = lword of xa * yb
 *      ql2 = lword of ya * xb
 *      ql3 = lword of qrsum
 * and qrsum = qr1 + qr2 + qr3, where
 *      qr1 = rword of xa * yb
 *      qr2 = rword of ya * xb
 *      qr3 = lword of xb * yb
 */
unsigned int mul_distrib( unsigned int x,
                          unsigned int y)
{
        unsigned int xa, xb, ya, yb, xaya, xayb, yaxb, xbyb;
        unsigned int qr1, qr2, qr3, qrsum, ql1, ql2, ql3, product;
        xa = x >> FXP_WORD_BITS;
        xb = (x & FXP_RWORD_MASK);
        ya = y >> FXP_WORD_BITS;
        yb = (y & FXP_RWORD_MASK);
        //if (VERBOSE) printf("\txa:%X, xb:%X, ya:%X, yb:%X\n", xa, xb, ya, yb);
        xayb = xa * yb;
        yaxb = ya * xb;
        qr1 = xayb & FXP_RWORD_MASK;
        qr2 = yaxb & FXP_RWORD_MASK;
        qr3 = xbyb >> FXP_WORD_BITS;
        //int rbit = (xbyb >> FXP_WORD_BITS_M1) & 1;
        qrsum = qr1 + qr2 + qr3; // + rbit;
        xaya = xa * ya;
        xbyb = xb * yb;
        ql1 = xayb >> FXP_WORD_BITS;
        ql2 = yaxb >> FXP_WORD_BITS;
        //rbit = ((qrsum & FXP_LWORD_MASK) \
        //              >> FXP_WORD_BITS_M1) & 1;
        ql3 = (qrsum >> FXP_WORD_BITS); // + rbit;
        product = xaya + ql1 + ql2 + ql3;
        return product;
}

/*
 * Safe implementation of fxp multiplication using only ints.
 *
 * Description of the distributive multiplication approach:
 * As in (a + x) * (b + y) = ab + ay + bx + xy,
 * for arguments fxp1 and fxp2 equal to the following:
 * fxp1 = (a << fxp_frac_bits) | x, and
 * fxp2 = (b << fxp_frac_bits) | y
 *
 * a, b their corresponding positive whole bits, and x, y their
 * positive frac bits, the product is calculated as:
 *
 * Frac part: frac_part(pfsum);
 *          where:
 *          pfsum = pf1 + pf2 + pf3
 *          pf1 = frac_part(a * y)
 *          pf2 = frac_part(b * x)
 *          pf3 = ((x' * y') >> frac_bits)
 *
 * Notice when all bits (except sign one) are used for the
 * frac part, the entire multiplication is effectively just pf3
 * (everything else will be zero.) Calculating pf3 avoiding
 * precision loss is particularly important for that case.
 *
 * Whole part: pwsum = pw1 + pw2 + pw3 + pw4
 *          where:
 *          pw1: (a * b)
 *          pw2: whole_part(a * y)
 *          pw3: whole_part(b * x)
 *          pw4: whole_part( pfsum )
 *
 * Notice that regardless of frac_bits, pw2 and pw3 can never
 * overflow, similarly to pf1 and pf2. However if pw1 > whole_max,
 * that already means overflow.
 * Also of course, if the sum pw1 + pw2 + pw3 + pw4 > whole_max,
 * then overflow as well.
 *
 * Works for systems in which sizeof(long) is not larger
 * than sizeof(int).
 * Does not use divisions.
 */
int fxp_mul(int fxp1, int fxp2)
{
        if (fxp1 == FXP_UNDEF || fxp2 == FXP_UNDEF)
                return FXP_UNDEF;
        if (fxp1 == FXP_POS_INF || fxp2 == FXP_POS_INF) {
                if (fxp1 == 0 || fxp2 == 0)
                        return FXP_UNDEF;
                else
                        return (fxp1 < 0 || fxp2 < 0)?
                                FXP_NEG_INF: FXP_POS_INF;
        }
        if (fxp1 == FXP_NEG_INF || fxp2 == FXP_NEG_INF) {
                if (fxp1 == 0 || fxp2 == 0)
                        return FXP_UNDEF;
                return (fxp1 > 0 || fxp2 > 0)?
                                FXP_NEG_INF: FXP_POS_INF;
        }
        int v1, v2, a, b, nba, nbb;
        int x, y, nbx, nby;
        int ab, ay, bx, xy;
        // Positive values of the arguments
        v1 = (fxp1 >= 0)? fxp1: -fxp1;
        v2 = (fxp2 >= 0)? fxp2: -fxp2;

        // Whole parts, and nbits in them
        a = fxp_get_whole_part(v1);
        b = fxp_get_whole_part(v2);
        nba = fxp_nbits(a);
        nbb = fxp_nbits(b);

        //printf("a:%d (%d bits),  b:%d (%d bits)\n", a, nba, b, nbb);
        if (nba + nbb > FXP_whole_bits) {
                // The product will for sure overflow just by
                // multiplying the whole parts.
                // Return appropriately signed infinity
                //printf("01. Overflowing!!!!!\n");
                return ((fxp1 >= 0 && fxp2 >= 0) \
                            || (fxp1 < 0 && fxp2 < 0))?
                            FXP_POS_INF: FXP_NEG_INF;
        }

        ab = a * b;
        // Frac parts, and nbits in them
        x = fxp_get_bin_frac(v1);
        y = fxp_get_bin_frac(v2);
        nbx = fxp_nbits(x);
        nby = fxp_nbits(y);
        // Compute pf1, pf2, pf3, and pfsum
        ay = a * y;
        bx = b * x;
        //printf("ab: %d\nay: %d\nbx: %d\n", ab, ay, bx);
        int pf1 = fxp_get_bin_frac(ay);
        int pf2 = fxp_get_bin_frac(bx);
        int pf3;

        // Calling mul_distrib to calculate pf3
        unsigned int ux = ((unsigned int) x) << FXP_whole_bits;
        unsigned int uy = ((unsigned int) y) << FXP_whole_bits;
        //printf("x:%X, ux:%X\n", x, ux);
        //printf("y:%X, uy:%X\n", y, uy);
        unsigned int pf3u = mul_distrib(ux, uy);
        int rbit = (pf3u >> FXP_whole_bits_m1) & 1;
        pf3 = (int) ((pf3u >> FXP_whole_bits) + rbit);

        //printf("pf3u: %u, pf3:%X\n", pf3u, pf3);;
        // We must sum safely because depending on the fxp config
        // there might not be enough whole bits to hold the
        // carry on from just the sum of these three frac pieces
        int pfsum = fxp_add(pf1, fxp_add(pf2, pf3));
        //printf("pf1:%X\npf2:%X\npf3:%x\n", pf1, pf2, pf3);
        //printf("pfsum:%X\n", pfsum);
        if (pfsum == FXP_POS_INF) {
                // Overflow. Return appropriately signed infinity
                return ((fxp1 >= 0 && fxp2 >= 0) || (fxp1 < 0 && fxp2 < 0))?
                            FXP_POS_INF: FXP_NEG_INF;
        }
        int pffrac = fxp_get_bin_frac(pfsum);

        // Compute remaining whole parts (pw2, pw3, pw4)
        int pw1 = fxp_get_whole_part(ay);
        int pw2 = fxp_get_whole_part(bx);
        int pw3 = fxp_get_whole_part(pfsum);
        // Sum all whole parts safely
        int pwsum = ab + pw1 + pw2 + pw3;
        if (pwsum > FXP_whole_max) {
                // Overflow.
                // Return appropriately signed infinity
                return ((fxp1 >= 0 && fxp2 >= 0) \
                            || (fxp1 < 0 && fxp2 < 0))?
                            FXP_POS_INF: FXP_NEG_INF;
        }
        //printf("pwsum:%d, pfsum_frac:%d\n", pwsum, fxp_get_bin_frac(pfsum));
        // No overflow, return the appropriately signed product
        int pproduct = fxp_bin(pwsum, pffrac);
        return ((fxp1 >= 0 && fxp2 >= 0) || (fxp1 < 0 && fxp2 < 0))?
                    pproduct: -pproduct;
}

/*
 * Safe implementation of fxp division using only integers.
 * This is software-based binary fxp division.
 * After 1st minuend it continues processing the dividend
 * bit by bit.
 * It works for systems in which sizeof(long) is not larger
 * than sizeof(int).
 */
int fxp_div(int fxp1, int fxp2)
{
        if (fxp1 == FXP_UNDEF || fxp2 == FXP_UNDEF)
                return FXP_UNDEF;
        if (fxp2 == 0)
                return (fxp1 == 0)? FXP_UNDEF:
                            (fxp1 > 0)? FXP_POS_INF: FXP_NEG_INF;

        // Positive values of the arguments
        int x = (fxp1 >= 0)? fxp1: -fxp1;
        int y = (fxp2 >  0)? fxp2: -fxp2;
        if (y == FXP_POS_INF)
                return (x == FXP_POS_INF)? FXP_UNDEF: 0;
        if (x == FXP_POS_INF)
                return ((fxp1 > 0 && fxp2 > 0) || \
                            (fxp1 < 0 && fxp2 < 0))?
                                FXP_POS_INF: FXP_NEG_INF;

        //printf("x:%d, mul upper lim:%d\n", x, fxp_mul(FXP_MAX, y));
        if (y <= FXP_frac_mask) {
                int mulcheck = fxp_mul(FXP_MAX, y);
                //printf("FXP_MAX : %10d, (x%X)\n", FXP_MAX, FXP_MAX);
                //printf("y       : %10d, (x%X)\n", y, y);
                //printf("mulcheck: %10d, (x%X)\n", mulcheck, mulcheck);
                // Due to rounding and border cases, using > is not
                // enough here, we need >=, otherwise another if
                // would be needed at end of division
                if (x >= mulcheck)
                    return ((fxp1 > 0 && fxp2 > 0) || \
                                (fxp1 < 0 && fxp2 < 0))?
                                    FXP_POS_INF: FXP_NEG_INF;
        }

        int bx = fxp_nbits(x);
        int by = fxp_nbits(y);

        int m, next_bit_mask;
        int bmax = bx + FXP_frac_bits;
        if (by == FXP_INT_BITS_M1) {
                // Border case, we need to shrink the divisor so that
                // the minuend can still have 1 more bit for the cases
                // when m < y
                y = (y >> 1);
                by--;
                //if (VERBOSE) printf("\n\t\tNew shrinked divisor %d\n\n", y);
                //qbits++;
                bmax--;
        }
        int bidx;
        int size_diff = bx - by;
        int shift = size_diff - 1;
        // Initialize starting minuend
        if (size_diff > 0) {
                int mmask = FXP_POS_INF >> (FXP_INT_BITS_M1 - by);
                m = (x & (mmask << size_diff)) >> size_diff;
                next_bit_mask = 1 << shift;
                bidx = by;
        } else {
                // Minor speedup, making the first minuend have as
                // many bits as the divisor
                m = x << -size_diff;
                next_bit_mask = 0;
                bidx = bx - size_diff;
        }
        //if (VERBOSE) printf("\td: starting m: %d (%d bits)\n", m, bidx);
        int q = 0;
        // bidx is (from left to right) the highest bit # we are
        // currently processing, so we loop till bidx exceeds the
        // right-most bit in the full (left-shifted) dividend
        while (bidx <= bmax) {
                //if (VERBOSE) printf("\tm: %d (hex:%x, %d bits) bidx:%d\n", m, m, bm, bidx);
                if (m >= y) {
                        // Append a 1 to the right of the quotient
                        q = (q << 1) | 1;
                        m = m - y;
                        //ba = 1;
                } else {
                        // Append a 0 to the right of the quotient
                        q = (q << 1);
                        //ba = 0;
                }
                //if (VERBOSE) {
                //    trace_fxp_div("div:", loops, fxp_frac_bits,
                //        bidx, x, y, q, ba, m, (m - ba * y));
                //}
                //m = difference;
                bidx++;
                if (next_bit_mask > 0) {
                        // Pull down next bit from the dividend,
                        // and append it to the right of the minuend
                        m = (m << 1) | \
                                    ((x & next_bit_mask) \
                                        >> shift);
                        next_bit_mask = (next_bit_mask >> 1);
                        shift--;
                } else {
                        // Pull down a 0 bit and append it to minuend
                        m = (m << 1);
                }
                //loops++;
        }
        // Return properly signed quotient
        return ((fxp1 >= 0 && fxp2 > 0) || (fxp1 < 0 && fxp2 < 0))?
                        q: -q;
}

inline unsigned int fxp_get_e()
{
        return FXP_shifted_e;
}

inline unsigned int fxp_get_pi()
{
        return FXP_shifted_pi;
}

inline unsigned int fxp_get_ln_2()
{
        return FXP_shifted_ln_2;
}

inline unsigned int fxp_get_lg10_2()
{
        return FXP_shifted_lg10_2;
}

/*
 * Calculates the characteristic c and full precision (non-shifted)
 * rounded mantissa m of lg2, leaving those two values separately
 * in arguments passed by reference.
 * Uses the BKM L-Mode algorithm (requires table of pre-calculated
 * values) and only ints.
 *
 * For more details on BKM:
 * https://en.wikipedia.org/wiki/BKM_algorithm
 */
static inline void fxp_lg2_sepcm(int fxp1, int * c, int * m)
{
        if (fxp1 <= 0) {
                *c = (fxp1 == 0)? FXP_NEG_INF: FXP_UNDEF;
                return;
        }
        if (fxp1 == FXP_POS_INF) {
                *c = FXP_POS_INF;
                return;
        }
        // Here fxp1 for sure > 0
        int clz = __builtin_clz((unsigned int) fxp1);
        // Notice clz > 0 since at least sign bit in fxp1 is 0
        int nbx = FXP_INT_BITS - clz;
        // Assign the characteristic to *c
        *c = ((fxp1 < FXP_one) || (fxp1 >= FXP_two))?
                (nbx - FXP_frac_bits - 1): 0;
        // Here we replicate what the lg2_l implementation does,
        // but simulating longs with two ints a and b, so
        // instead of a k of type long, two ints (ka, kb),
        // ka the most significant, kb the least one.
        // We do this for x, z, xs, y, and aux
        // (argumentb would always be 0, so we just need
        // the a part for argument)
        unsigned int argument = ((unsigned int) fxp1) << (clz - 1);
        unsigned int xa = FXP_BKM_ONE, xb = 0;      // x
        unsigned int za, zb;                        // z
        unsigned int xsa = 0, xsb = 0;              // xs
        unsigned int ya = 0, yb = 0;                // y
        unsigned int auxa, auxb;                    // aux
        unsigned int a_bits;
        unsigned int ab_mask = 1;
        //printf("c:%d argument:x%X\n", c, argument);
        for (int shift = 1; shift < FXP_INT_BITS; shift++) {
                //unsigned long xs = (x >> shift);
                // Here we must r-shift the combo (xa, xb) by shift units
                // leaving the results in combo (xsa, xsb)
                a_bits = (xa & ab_mask);
                xsa = (xa >> shift);
                xsb = (a_bits << (FXP_INT_BITS - shift)) \
                        | (xb >> shift);
                ab_mask = (ab_mask << 1) | 1;
                // z = x + xs;
                za = xa + xsa;
                if (xsb > FXP_BKM_MAXU - xb)
                        // Carry over from zb into za
                        za++;
                zb = xb + xsb;
                // if (z <= argument) {
                if ((za < argument) || \
                    ((za == argument) && (zb == 0))) {
                        // x = z;
                        xa = za;
                        xb = zb;
                        // Here we use the precalculated values
                        // y += FXP_BKM_LOGS_L[shift];
                        auxa = FXP_BKM_LOGS[shift];
                        auxb = FXP_BKM_LOGS_XTRA[shift];
                        ya += auxa;
                        if (yb > FXP_BKM_MAXU - auxb)
                                // Carry from yb into ya
                                ya++;
                        yb += auxb;
                        //printf("\tx:x%X-%X  Updating y:x%X-%X\n", \
                        //        xa, xb, ya, yb);
                }
        }
        // One final "iteration", as if shift == FXP_INT_BITS
        xsa = 0;
        xsb = xa;
        za = xa;
        if (xsb > FXP_BKM_MAXU - xb) za++;
        zb = xb + xsb;
        if ((za < argument) || \
            ((za == argument) && (zb == 0))) {
                xa = za;
                xb = zb;
                auxa = FXP_BKM_LOGS[FXP_INT_BITS];
                auxb = FXP_BKM_LOGS_XTRA[FXP_INT_BITS];
                ya += auxa;
                if (yb > FXP_BKM_MAXU - auxb) ya++;
                yb += auxb;
        }
        // Now (ya, yb) has the "long" calculated mantissa.
        // Assign the full precision (non-shifted) int-sized
        // rounded mantissa to *m
        int rbit = (yb >> FXP_INT_BITS_M1);
        *m = ya + rbit;
        return;
}

/*
 * Default log2 using BKM and only ints.
 */
int fxp_lg2(int fxp1)
{
        int c, y;
        // Get the separate characteristic and full mantissa
        fxp_lg2_sepcm(fxp1, &c, &y);
        if ((c == FXP_POS_INF) || (c == FXP_NEG_INF) \
                || (c == FXP_UNDEF)) {
                return c;
        }
        // Shift the mantissa (rounding it) for current fxp config
        int rbit = (y >> FXP_whole_bits_m2) & 1u;
        int m = (y >> FXP_whole_bits_m1) + rbit;
        //printf("Final mantissa:%d  (rounding bit was %d)\n", m, rbit);
        // Return the complete logarithm (c + m) as fxp
        if (c < FXP_whole_min) {
                if ((c + 1 == FXP_whole_min) && (m > 0)) {
                        // Special case: in spite of not enough whole
                        // bits available for c, there are for c + 1,
                        // which is what this logarithm will return as
                        // whole part
                        return INT_MIN + ((int) m);
                }
                // Overflow: the characteristic of the logarithm
                // will not fit in the whole bits part of our
                // current fxp settings
                return FXP_NEG_INF;
        }
        if (c >= 0)
                return (c << FXP_frac_bits) + ((int) m);
        else
                return (-(-c << FXP_frac_bits)) + ((int) m);
}

/*
 * Default implementation of ln()
 * Uses the default lg2() and only ints:
 * Calculates ln(x) as lg2(x) * ln(2)
 */
int fxp_ln(int fxp1)
{
        return 0;
        int c2, y2;
        fxp_lg2_sepcm(fxp1, &c2, &y2);
        if ((c2 == FXP_POS_INF) || (c2 == FXP_NEG_INF) \
                || (c2 == FXP_UNDEF)) {
                return c2;
        }

        // Calculate the characteristic for ln

        // Calculate mantissa for ln

/*
        // Return the final ln(x)
        // Shift the mantissa (rounding it) for current fxp config
        int rbit = (y >> FXP_whole_bits_m2) & 1u;
        int m = (y >> FXP_whole_bits_m1) + rbit;
        //printf("Final mantissa:%d  (rounding bit was %d)\n", m, rbit);
        // Return the complete logarithm (c + m) as fxp
        if (c < FXP_whole_min) {
                if ((c + 1 == FXP_whole_min) && (m > 0)) {
                        // Special case: in spite of not enough whole
                        // bits available for c, there are for c + 1,
                        // which is what this logarithm will return as
                        // whole part
                        return INT_MIN + ((int) m);
                }
                // Overflow: the characteristic of the logarithm
                // will not fit in the whole bits part of our
                // current fxp settings
                return FXP_NEG_INF;
        }
        if (c >= 0)
                return (c << FXP_frac_bits) + ((int) m);
        else
                return (-(-c << FXP_frac_bits)) + ((int) m);


        unsigned int f = 0; 
        f >>= FXP_whole_bits;
        int rw = 0, rf = 0;
        // TODO:
        //int w = mul_whole_frac(fxp_get_whole_part(fxp1), FXP_LN_2_I32, \
        //                        &rw, &rf);
        //return fxp_add(rw, rf);

        //return fxp_mul(fxp_lg2(fxp1), FXP_shifted_ln_2);
*/
}

/*
 * Default implementation of lg10()
 * Uses the default lg2() and only ints
 */
int fxp_lg10(int fxp1)
{
        return 0;
        return fxp_mul(fxp_lg2(fxp1), FXP_shifted_lg10_2);
}


/*
 * Analogous to fxp_pow2_l but using only ints
 */
int fxp_pow2(int fxp1)
{
        if (fxp1 == FXP_UNDEF) return FXP_UNDEF;
        if (fxp1 == FXP_NEG_INF) return 0;
        if (fxp1 == FXP_POS_INF) return FXP_POS_INF;
        int w = fxp_get_whole_part(fxp1);
        int f = fxp_get_bin_frac(fxp1);
        unsigned int pow2w, argument;
        if (fxp1 >= 0) {
                if (w >= FXP_whole_bits_m1) {
                        return FXP_POS_INF;
                }
                pow2w = FXP_one << (w + 1);
                // Argument will be in [0, 1)
                argument = ((unsigned int) f) << FXP_whole_bits_m1;
        } else {
                if (w <= FXP_INT_BITS_M1_NEG) {
                        return 0;
                }
                //if (w < 0)
                //        pow2w = FXP_one >> (-w - 1);
                //else
                //        pow2w = FXP_one << 1;
                if (w < 0)
                        pow2w = FXP_one >> (-w);
                else
                        pow2w = FXP_one;
                argument = (FXP_one + f) \
                                << FXP_whole_bits_m1;
        }
        //if (VERBOSE) printf("pow2w:%X  argument:%X\n", pow2w, argument);
        //unsigned int x = FXP_BKM_ONE, y = 0;
        unsigned int xa = FXP_BKM_ONE, xb = 0;
        unsigned int ya = 0, yb = 0;
        unsigned int auxa, auxb, za, zb, xsa, xsb, a_bits, ab_mask;
        for (int k = 0; k < FXP_INT_BITS; k++) {
                //unsigned int const  z = y + FXP_BKM_LOGS[k];
                auxa = FXP_BKM_LOGS[k];
                auxb = FXP_BKM_LOGS_XTRA[k];
                za = ya + auxa;
                if (yb > FXP_BKM_MAXU - auxb)
                        za++;
                zb = yb + auxb;

                //if (VERBOSE) printf("k:%d,  z:%X-%X\n", k, za, zb);
                //if (z <= argument) {
                if ((za < argument) \
                    || ((za == argument) && (zb == 0))) {
                        //y = z;
                        ya = za;
                        yb = zb;

                        //x = x + (x >> k);
                        // Equivalent to:
                        // xs = (x >> k);
                        // x += xs;
                        // Here we r-shift the combo (xa, xb) by k units
                        // leaving the results in combo (xsa, xsb)
                        //if (k < FXP_INT_BITS) {
                        ab_mask = (1u << k) - 1;
                        a_bits = (xa & ab_mask);
                        xsa = (xa >> k);
                        xsb = (a_bits << (FXP_INT_BITS - k)) \
                                | (xb >> k);
                        //} else {
                        //        xsa = 0;
                        //        xsb = (xa >> (k - FXP_INT_BITS));
                        //}
                        xa += xsa;
                        if (xsb > FXP_BKM_MAXU - xb)
                                xa++;
                        xb += xsb;
                        //if (VERBOSE) printf("\tUpdating y (%X-%X) und x (%X-%X)\n", \
                        //                ya, yb, xa, xb);
                }
        }
        unsigned int x;
        if ((xa >> FXP_INT_BITS_M1) == 0) {
                x = (xa << 1) | (xb >> FXP_INT_BITS_M1);
        } else {
                x = xa;
                pow2w = pow2w << 1;
        }
        //if (VERBOSE) printf("xa:%X, final x:%X\n", xa, x);
        unsigned int md = mul_distrib(pow2w, x);
        return (int) md;
}

int fxp_pow(int x, int y)
{
        // TODO
        return 0;
}

int fxp_exp(int fxp1)
{
        // TODO
        return 0;
}

/*
 * Square root implementation
 */
int fxp_sqrt(int fxp1)
{
        // TODO
        return 0;
}
