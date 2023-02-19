/*
 * fxp.c
 *
 * An implementation of Binary Fixed Point numbers
 * encoding them into integers, together with their
 * arithmetic operations +, -, *, and /, and more.
 *
 * Safe arithmetic operations for the fxp's compliant
 * with INT32-C, as in:
 * https://wiki.sei.cmu.edu/confluence/display/c/INT32-C.+Ensure+that+operations+on+signed+integers+do+not+result+in+overflow
 *
 * By Raul Saavedra, Bonn, Germany
 *
 */

#include "fxp.h"
#include "fxp_tconst.h"
//#include "fxp_aux.h"
//#include <stdio.h>
//#include <stdlib.h>
//#include <assert.h>

//Used while testing and debugging trying to optimize division
//#include "fxp_aux.h"
//#define VERBOSE 1

#define FXP_FRAC_BITS_MIN 4
#define FXP_FRAC_BITS_DEF 16

// Allowing for no whole part, so other than the sign bit,
// all bits used for fraction part. This use case represents
// the range [-0.999..., 0.999...], or equivalently: (-1, 1)
// (e.g. actual whole values -1 and 1 outside the valid range)
#define FXP_FRAC_BITS_MAX FXP_INT_BITS_M1

#define FXP_FRAC_MAX_DEC 9999999

// Default number of bits to use for the frac part.
int FXP_frac_bits = FXP_FRAC_BITS_DEF;
static int FXP_frac_bits_mod2 = FXP_FRAC_BITS_DEF % 2;
// Can be changed dynamically calling fxp_set_frac_bits()
// but it remains unchanged until calling that function
// again. All such static variables which will only change
// if modifying the number of frac bits are named here
// starting with "FXP_" (uppercase,) but the rest of the
// variable name in lowercase

// For improved-precision version of fxp_mul
static int FXP_frac_mshift = FXP_FRAC_BITS_DEF / 2;
static int FXP_frac_maskl = 4032;
static int FXP_frac_maskr = 63;

// Default number of bits for the whole (including sign) part
int FXP_whole_bits = FXP_INT_BITS - FXP_FRAC_BITS_DEF;
static int FXP_whole_bits_m1 = FXP_INT_BITS - FXP_FRAC_BITS_DEF - 1;

// FXP_FRAC_MASK should correspond to 2^FXP_FRAC_BITS - 1
int FXP_frac_mask = ((1 << FXP_FRAC_BITS_DEF) - 1);
int FXP_frac_max = ((1 << FXP_FRAC_BITS_DEF) - 1);

// Default max and min valid values for the whole part of the fxp's
int FXP_whole_max = FXP_MAX >> FXP_FRAC_BITS_DEF;
int FXP_whole_min = -(FXP_MAX >> FXP_FRAC_BITS_DEF);

// Default max and min valid values for the conversion types
float FXP_max_f = ((float) (FXP_MAX >> FXP_FRAC_BITS_DEF)) \
                        + (((float) ((1 << FXP_FRAC_BITS_DEF) - 1)) \
                            / (((unsigned int) 1) << FXP_FRAC_BITS_DEF));
float FXP_min_f = -((float) (FXP_MAX >> FXP_FRAC_BITS_DEF)) \
                        - (((float) ((1 << FXP_FRAC_BITS_DEF) - 1)) \
                            / (((unsigned int) 1) << FXP_FRAC_BITS_DEF));

double FXP_max_d = ((double) (FXP_MAX >> FXP_FRAC_BITS_DEF)) \
                        + (((double) ((1 << FXP_FRAC_BITS_DEF) - 1)) \
                            / (((unsigned long) 1) << FXP_FRAC_BITS_DEF));
double FXP_min_d = -((double) (FXP_MAX >> FXP_FRAC_BITS_DEF)) \
                        - (((double) ((1 << FXP_FRAC_BITS_DEF) - 1)) \
                            / (((unsigned long) 1) << FXP_FRAC_BITS_DEF));

long double FXP_max_ld = ((long double) (FXP_MAX >> FXP_FRAC_BITS_DEF)) \
                        + (((long double) ((1 << FXP_FRAC_BITS_DEF) - 1)) \
                            / (((unsigned long long) 1) << FXP_FRAC_BITS_DEF));

long double FXP_min_ld = -((long double) (FXP_MAX >> FXP_FRAC_BITS_DEF)) \
                        - (((long double) ((1 << FXP_FRAC_BITS_DEF) - 1)) \
                            / (((unsigned long long) 1) << FXP_FRAC_BITS_DEF));

// For the BKM log calculation when using only ints.
// Values in this array are: a[k] = log2(1 + 1/2^k) represented as
// unsigned fxp's (4 bytes) with 30 frac bits
static const unsigned int FXP_BKM_LOGS[] = {
        0x40000000, 0x2570068E, 0x149A784B, 0xAE00D1C,
        0x598FDBE,  0x2D75A6E,  0x16E7968,  0xB7F285,
        0x5C2711,   0x2E1F07,   0x171265,   0xB89EB,
        0x5C523,    0x2E29D,    0x17151,    0xB8A9,
        0x5C54,     0x2E2A,     0x1715,     0xB8A,
        0x5C5,      0x2E2,      0x171,      0xB8,
        0x5C,       0x2E,       0x17,       0xB,
        0x5,        0x2,        0x1,        0x0,
        0x0,        0x0,        0x0,        0x0
};

// 32 extra bits following each of the values in the previous array
static const unsigned int FXP_BKM_LOGS_XTRA[] = {
        0x00000000, 0x7EF5A27D, 0xCD1B8B50, 0xFDEB43FB,
        0xB244C5B5, 0xB1DFB0F1, 0x5C2D229E, 0xB778428E,
        0xB5EAB1DE, 0xFE14EACA, 0x3743F454, 0x17BCABE1,
        0xB0A86FF2, 0x623F4A6C, 0x93B17D35, 0x82801725,
        0xEF6A3E08, 0x833FB72C, 0x44828311, 0xA2F9EB95,
        0x51AB2053, 0xA8E11ACC, 0x5473700F, 0xAA3A70B1,
        0x551D6683, 0x2A8EBECC, 0x15476248, 0x8AA3B1DD,
        0xC551D91C, 0xE2A8EC99, 0x7154764F, 0xB8AA3B28,
        0x5C551D94, 0x2E2A8ECA, 0x17154765, 0xB8AA3B2
};

// fxp rounded representations of transcendental constants using
// the default number of frac bits
#define FXP_DEF_SHIFT (FXP_INT_BITS_M1 - FXP_FRAC_BITS_DEF)
int FXP_shifted_e = FXP_E_I32 >> (FXP_DEF_SHIFT - 2);
int FXP_shifted_pi = FXP_PI_I32 >> (FXP_DEF_SHIFT - 2);
int FXP_shifted_ln_2 = FXP_LN_2_I32 >> FXP_DEF_SHIFT;
int FXP_shifted_lg10_2 = FXP_LG10_2_I32 >> FXP_DEF_SHIFT;

#define FXP_BKM_PREC (FXP_INT_BITS - 2)
static const int FXP_BKM_ONE = (1 << FXP_BKM_PREC);
static const unsigned int FXP_BKM_HALF = (FXP_BKM_ONE >> 1);
static const unsigned int FXP_BKM_MSB1 = (1 << (FXP_INT_BITS - 1));
static const unsigned int FXP_BKM_IMSB1 = ~FXP_BKM_MSB1;
static const unsigned int FXP_BKM_NMSB1 = (FXP_BKM_MSB1 >> 1);
static const unsigned int FXP_BKM_MAXU = ~((unsigned int) 0);

// Auxiliary variables used in the lg2 implementations
int FXP_half = 1 << (FXP_FRAC_BITS_DEF - 1);
int FXP_one = 1 << FXP_FRAC_BITS_DEF;
int FXP_two = 1 << (FXP_FRAC_BITS_DEF + 1);
static int FXP_lg2_maxloops = FXP_FRAC_BITS_DEF;
static int FXP_lg2_yashift = FXP_BKM_PREC - FXP_FRAC_BITS_DEF;
static int FXP_lg2_ybmask = 0;
static int FXP_lg2_ybshift = 0;

// Default desired max frac decimal value
// (can be changed dynamically calling set_[auto_]frac_max_dec
int fxp_frac_max_dec = 9999;

// Auxiliary variables for the implementations that use longs (fxp_l.c)
const int FXP_BKM_PREC_L = ((sizeof(long) * 8) - 2);
const int FXP_BKM_L_CLZSHIFT = (FXP_BKM_PREC_L - FXP_BKM_PREC - 1);
const unsigned long FXP_BKM_ONE_L = \
                            ((unsigned long) 1) << FXP_BKM_PREC_L;
const unsigned long FXP_BKM_HALF_L = (FXP_BKM_ONE_L >> 1);
unsigned long FXP_max_lshifted = (FXP_MAX_L) << FXP_FRAC_BITS_DEF;
int FXP_lg2_l_maxloops = FXP_FRAC_BITS_DEF + 1;
int FXP_lg2_l_shift = FXP_BKM_PREC_L - FXP_FRAC_BITS_DEF;

/*
 * Given an fxp with x number of frac bits, returns
 * the rounded representation using y frac bits
 */
int fxp_change_nfracbits(int fxp, int x, int y)
{
        if ((fxp <= FXP_NEG_INF) || (fxp == FXP_POS_INF))
                return fxp;
        if (y < FXP_FRAC_BITS_MIN)
                y = FXP_FRAC_BITS_MIN;
        else if (y > FXP_FRAC_BITS_MAX)
                y = FXP_FRAC_BITS_MAX;
        int pfxp = (fxp > 0? fxp: -fxp);
        int shift = x - y;
        if (shift <= 0) {
                if (fxp_nbits(pfxp) - shift > FXP_INT_BITS_M1) {
                        return (fxp > 0)? FXP_POS_INF: FXP_NEG_INF;
                }
                return (fxp << (-shift));
        }
        int mask = (1 << shift) - 1;
        int frac_lost = pfxp & mask;
        int shifted = pfxp >> shift;
        // Round number by adding one to truncated value if the
        // magnitude in lost bits is >= half the largest number
        // representable in those lost bits
        int rounded = shifted + ((frac_lost >= ((mask + 1) / 2))? 1: 0);
        //printf("pfxp:%x, shift:%d, mask:%x, frac_lost:%x, shifted:%x, rounded:%x\n", \
        //        pfxp, shift, mask, frac_lost, shifted, rounded);
        return (fxp >= 0)? rounded: -rounded;
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
        FXP_frac_bits = (nfracbits < FXP_FRAC_BITS_MIN? FXP_FRAC_BITS_MIN:
                            (nfracbits > FXP_FRAC_BITS_MAX?
                                FXP_FRAC_BITS_MAX: nfracbits));
        FXP_frac_bits_mod2 = FXP_frac_bits % 2;

        FXP_whole_bits = FXP_INT_BITS - FXP_frac_bits;
        FXP_whole_bits_m1 = FXP_whole_bits - 1;

        // fxp_frac_mask should correspond to 2^FXP_FRAC_BITS - 1
        FXP_frac_mask = (1 << FXP_frac_bits) - 1;

        // Variables used in fxp_mul to process the
        // full multiplication of the frac parts avoiding
        // precision loss
        if (FXP_frac_bits == FXP_INT_BITS_M1) {
                FXP_frac_mshift = (FXP_INT_BITS_M1 - 1) / 2;
        } else {
                FXP_frac_mshift = (FXP_frac_bits / 2) + (FXP_frac_bits % 2);
        }
        FXP_frac_maskr = (1 << FXP_frac_mshift) - 1;
        FXP_frac_maskl = FXP_frac_maskr << FXP_frac_mshift;

        // When using all bits for frac (except for the sign bit),
        // then our max valid frac cannot be equal to frac_mask
        // (because in that case that value is already the largest
        // positive integer == POS_INF), so we must substract one
        // from the frac mask to get the largest valid frac value
        FXP_frac_max = FXP_frac_mask - (FXP_whole_bits == 1? 1: 0);

        // Auxiliary variables for fxp_l.c implementations
        FXP_max_lshifted = ((FXP_MAX_L) << FXP_frac_bits);

        // Max and min valid values for the whole part of the fxp's
        FXP_whole_max = FXP_MAX >> FXP_frac_bits;
        FXP_whole_min = (-FXP_whole_max);

        // Max and min floating-point conversion values
        FXP_max_f = ((float) FXP_whole_max) + \
                        ((float) FXP_frac_max) / \
                            ((float) (((unsigned int) 1) << FXP_frac_bits));
        FXP_max_d = ((double) FXP_whole_max) + \
                        ((double) FXP_frac_max) / \
                            ((double) (((unsigned int) 1) << FXP_frac_bits));
        FXP_max_ld = (long double) FXP_max_d;
        FXP_min_f = -FXP_max_f;
        FXP_min_d = -FXP_max_d;
        FXP_min_ld = -FXP_max_ld;

        // Adjust precision of e, pi, etc. to the frac bits in use
        FXP_shifted_e = fxp_change_nfracbits(FXP_E_I32, \
                        FXP_INT_BITS_M1 - 2, FXP_frac_bits);
        FXP_shifted_pi = fxp_change_nfracbits(FXP_PI_I32, \
                        FXP_INT_BITS_M1 - 2, FXP_frac_bits);
        FXP_shifted_ln_2 = fxp_change_nfracbits(FXP_LN_2_I32, \
                        FXP_INT_BITS_M1, FXP_frac_bits);
        FXP_shifted_lg10_2 = fxp_change_nfracbits(FXP_LG10_2_I32, \
                        FXP_INT_BITS_M1, FXP_frac_bits);

        // Auxiliary variables used for lg2 calculations
        FXP_one = fxp(1);
        FXP_half = FXP_one >> 1;
        FXP_two = FXP_one << 1;
        FXP_lg2_maxloops = FXP_frac_bits + 1;
        //FXP_lg2_l_maxloops = MIN(FXP_INT_BITS, FXP_frac_bits + 1);
        FXP_lg2_l_maxloops = FXP_lg2_maxloops;
        FXP_lg2_l_shift = FXP_BKM_PREC_L - FXP_frac_bits;
        FXP_lg2_yashift = FXP_BKM_PREC - FXP_frac_bits;
        if (FXP_lg2_yashift >= 0) {
            FXP_lg2_ybshift = 0;
            FXP_lg2_ybmask = 0;
        } else {
            // When lg2_shift is negative, lg2_ybmask is equal to
            // a number of 1's == -lg2_shift, shifted to the
            // left-most bits of an int
            FXP_lg2_ybshift = FXP_INT_BITS + FXP_lg2_yashift;
            FXP_lg2_ybmask = ((1 << (-FXP_lg2_yashift)) - 1) \
                    << FXP_lg2_ybshift;
        }

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
    return fxp_frac_max_dec;
}

/*
 * Return the max decimal value 99...9 to be used as
 * max posible decimal fraction value
 */
int fxp_get_frac_max_dec()
{
    return fxp_frac_max_dec;
}

/*
 * Create an fxp number given its whole part only
 */
int fxp(int whole)
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
 * Note last example frac is not 5 but 500. (5 would correspond to 24.005)
 * This decimal value gets scaled into the binary range available for frac.
 * For negative numbers with whole part=0, the frac value must be negative.
 * If frac is too large it will get truncated to its most
 * significant decimal digits until under the value of FXP_FRAC_MAX_DEC,
 * e.g. for a max of 999, a frac=987654 will be truncated to 987
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

int fxp_get_whole_part(int fxp)
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
int fxp_get_bin_frac(int fxp)
{
        if (fxp < 0)
                if (fxp <= FXP_NEG_INF)
                        // In the whole-part function above we are
                        // returning the same whole parts for both -INF
                        // and UNDEF. Here for -INF we return the value
                        // used, but for UNDEF something different
                        // and odd in order to enable differenciating
                        // -INF from UNDEF when checking their whole and
                        // frac constituents. It must be something that
                        // no other fxp could ever have, not even the
                        // +/-INF values. That can be some (alleged) frac
                        // bits that are actually invalid, e.g. outside
                        // the range of valid frac bits, like -frac_max -
                        // 1. Notice that this value logically
                        // "overflows" the frac bits, but in a way makes
                        // sense to mark and differenciate precisely only
                        // UNDEF this way. Conveniently, when the
                        // # of frac bits is 31, this scheme returns
                        // the full FXP_UNDEF as the frac part.
                        return -FXP_frac_mask - (fxp == FXP_UNDEF? 1: 0);
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
        long num, denom, ldivision;
        denom = (long) FXP_frac_max;
        int positive_frac;
        if (fxp < 0) {
                positive_frac = (-fxp) & FXP_frac_mask;
                if (positive_frac == 0) return 0;
                num = -(((long) positive_frac + 1) * \
                                ((long) fxp_frac_max_dec));
        }
        else {
                positive_frac = fxp & FXP_frac_mask;
                if (positive_frac == 0) return 0;
                num = ((long) positive_frac + 1)
                                * ((long) fxp_frac_max_dec);
        }
        ldivision = num / denom;
        if (ldivision > fxp_frac_max_dec) ldivision = fxp_frac_max_dec;
        return (int) ldivision;
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
int fxp_sub(int fxp1, int fxp2)
{
        // Conservative check to never attempt to change sign
        // to the true most negative int (aka. FXP_UNDEF)
        if (fxp2 == FXP_UNDEF) return FXP_UNDEF;
        return fxp_add(fxp1, -fxp2);
}

/*
Return number of bits used by x, that is,
most significant bit in x that is a 1
(returned value is between 0 and FXP_INT_BITS)
*/
int fxp_nbits(unsigned int x)
{
    if (x == 0) return 0;
    // Replacing with gcc builtin function that counts leading zeros
    return FXP_INT_BITS - __builtin_clz(x);
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
 * that already means the multiplication overflows.
 * Also of course, if the sum pw1 + pw2 + pw3 + pw4 > whole_max,
 * then overflow as well.
 *
 * Works for systems in which sizeof(long) is not larger
 * than sizeof(int).
 * Also it does not use divisions to check for overflows.
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
                return ((fxp1 >= 0 && fxp2 >= 0) || (fxp1 < 0 && fxp2 < 0))?
                            FXP_POS_INF: FXP_NEG_INF;
        }

        int pw1 = a * b;
        // Frac parts, and nbits in them
        x = fxp_get_bin_frac(v1);
        y = fxp_get_bin_frac(v2);
        nbx = fxp_nbits(x);
        nby = fxp_nbits(y);
        // Compute pf1, pf2, pf3, and pfsum
        ay = a * y;
        bx = b * x;
        int pf1 = fxp_get_bin_frac(ay);
        int pf2 = fxp_get_bin_frac(bx);
        int pf3;
        if (nbx + nby <= FXP_INT_BITS_M1) {
            // We have room to simply multiply x * y, then shift
            pf3 = (x * y) >> FXP_frac_bits;
        } else {
            /*
             * Below an improved method to compute pf3 when we do not
             * have enough room in an int to multiply x * y.
             * This method uses again a distributive scheme, not with
             * whole vs. frac parts (we are here all within fraction
             * parts after all,) but with left and right chunks of the
             * frac parts. For subshift = frac_bits / 2:
             * x = (xl << subshift) | xr
             * y = (yl << subshift) | yr
             *
             * rchunk = R part of qrsum
             *      where:
             *          qrsum = qr1 + qr2 + qr3
             *          qr1 = R part of (xl * yr)
             *          qr2 = R part of (yl * xr)
             *          qr3 = (xr * yr) >> subshift
             * lchunk: qlsum = ql1 + ql2 + ql3 + ql4
             *      where:
             *          ql1: (xl * yl)
             *          ql2: L part of (xl * yr)
             *          ql3: L part of (yl * xr)
             *          ql4: L part of qrsum
             *
             * Final pf3 will then  lchunk == qlsum, the most
             * significant part of the product of the frac parts
             */
            int submaskl, submaskr, subshift;
            int oddnfb;
            if (FXP_frac_bits == FXP_INT_BITS_M1) {
                // Eliminate the very lowest significant bit of each
                // operand to make sure none of the subchunk products
                // ever exceed FXP_INT_BITS_M1 bits.
                // (We must shift back this lost bit in the final pf3)
                oddnfb = -1;
                x = x >> 1;
                y = y >> 1;
            } else {
                oddnfb = FXP_frac_bits_mod2;
                x = x << oddnfb;
                y = y << oddnfb;
            }
            // Variables fxp_frac_maskl, _maskr, and _mshift are
            // all preinitialized for the default # of frac bits to use,
            // and readjusted when fxp_set_frac_bits() gets called
            int xl = (x & FXP_frac_maskl) >> FXP_frac_mshift;
            int xr = (x & FXP_frac_maskr);
            int yl = (y & FXP_frac_maskl) >> FXP_frac_mshift;
            int yr = (y & FXP_frac_maskr);
            int xlyr = xl * yr;
            int ylxr = yl * xr;
            int qr1 = xlyr & FXP_frac_maskr;
            int qr2 = ylxr & FXP_frac_maskr;
            int qr3 = (xr * yr) >> FXP_frac_mshift;
            int qrsum = qr1 + qr2 + qr3;
            int ql1 = xl * yl;
            int ql2 = (xlyr & FXP_frac_maskl) >> FXP_frac_mshift;
            int ql3 = (ylxr & FXP_frac_maskl) >> FXP_frac_mshift;
            int ql4 = (qrsum & FXP_frac_maskl) >> FXP_frac_mshift;
            pf3 = ql1 + ql2 + ql3 + ql4;
            // Now final adjustment shift depending on value of oddnfb
            if (oddnfb > 0) {
                pf3 = pf3 >> 1;
            } else if (oddnfb < 0) {
                pf3 = pf3 << 1;
            }
        }
        // We must sum safely since there might not be enough whole
        // bits to hold the carry on from just the sum of these three
        // frac pieces
        int pfsum = fxp_add(pf1, fxp_add(pf2, pf3));
        if (pfsum == FXP_POS_INF) {
                // Overflow. Return appropriately signed infinity
                return ((fxp1 >= 0 && fxp2 >= 0) || (fxp1 < 0 && fxp2 < 0))?
                            FXP_POS_INF: FXP_NEG_INF;
        }

        // Compute remaining whole parts (pw2, pw3, pw4)
        int pw2 = fxp_get_whole_part(ay);
        int pw3 = fxp_get_whole_part(bx);
        int pw4 = fxp_get_whole_part(pfsum);
        // Sum all whole parts safely
        int pwsum = pw1 + pw2 + pw3 + pw4;
        if (pwsum > FXP_whole_max) {
                // Overflow.
                // Return appropriately signed infinity
                return ((fxp1 >= 0 && fxp2 >= 0) || (fxp1 < 0 && fxp2 < 0))?
                            FXP_POS_INF: FXP_NEG_INF;
        }
        // No overflow, return the appropriately signed product
        int pproduct = fxp_bin(pwsum, fxp_get_bin_frac(pfsum));
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
                return ((fxp1 > 0 && fxp2 > 0) || (fxp1 < 0 && fxp2 < 0))?
                                FXP_POS_INF: FXP_NEG_INF;

        if ((y <= FXP_frac_mask) &&
                (x > fxp_mul(FXP_MAX, y))) {
                return ((fxp1 > 0 && fxp2 > 0) || (fxp1 < 0 && fxp2 < 0))?
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

int fxp_get_e()
{
        return FXP_shifted_e;
}

int fxp_get_pi()
{
        return FXP_shifted_pi;
}

int fxp_get_ln_2()
{
        return FXP_shifted_ln_2;
}

int fxp_get_lg10_2()
{
        return FXP_shifted_lg10_2;
}

/*
 * Default implementation of ln()
 * Uses the default lg2() and only ints
 */
int fxp_ln(int fxp1)
{
        return fxp_mul(fxp_lg2(fxp1), FXP_shifted_ln_2);
}

/*
 * Default implementation of lg10()
 * Uses the default lg2() and only ints
 */
int fxp_lg10(int fxp1)
{
        return fxp_mul(fxp_lg2(fxp1), FXP_shifted_lg10_2);
}

/*
 * Default log2 calculation using only ints and the BKM algorithm.
 * Requires table of pre-calculated values
 *
 * For more details on BKM:
 * https://en.wikipedia.org/wiki/BKM_algorithm
 */
int fxp_lg2(int fxp1)
{
        if (fxp1 <= 0) return ((fxp1 == 0)? FXP_NEG_INF: FXP_UNDEF);
        if (fxp1 == FXP_POS_INF) return FXP_POS_INF;
        // Here fxp1 for sure > 0
        int clz = __builtin_clz((unsigned int) fxp1);
        // clz > 0 since at least sign bit == 0
        int nbx = FXP_INT_BITS - clz;
        // c is characteristic
        int c = ((fxp1 < FXP_one) || (fxp1 >= FXP_two))?
                    (nbx - FXP_frac_bits - 1): 0;

        // Here we replicate what the lg2_l implementation does,
        // but simulating longs with two ints a and b, so
        // instead of a K of type long, two ints (ka, kb),
        // ka the most significant, kb the least one
        unsigned int za = ((unsigned int) fxp1) << (clz - 1);
        //unsigned int zb = 0; // <- zb not really needed, will remain 0

        unsigned int xa = FXP_BKM_ONE, xb = 0;      // x
        unsigned int zza, zzb;                      // zz
        unsigned int xsa = 0, xsb = 0;              // xs
        unsigned int ya = 0, yb = 0;                // y
        unsigned int auxa, auxb;                    // aux
        unsigned int a_bits;
        unsigned int ab_mask = 1;
        for (int shift=1; shift <= FXP_lg2_maxloops; shift++) {

                //unsigned long xs = (x >> shift);
                // Here we must r-shift the combo (xa, xb) by shift units
                // leaving the results in combo (xsa, xsb)
                if (shift < FXP_INT_BITS) {
                        a_bits = (xa & ab_mask);
                        xsa = (xa >> shift);
                        xsb = (a_bits << (FXP_INT_BITS - shift)) | (xb >> shift);
                        ab_mask = (ab_mask << 1) | 1;
                } else {
                        xsa = 0;
                        xsb = (xa >> (shift - FXP_INT_BITS));
                }

                // zz = x + xs;
                zza = xa + xsa;
                if (xsb > FXP_BKM_MAXU - xb)
                    // Carry over from zzb into zza
                    zza++;
                zzb = xb + xsb;

                // if (zz <= z) {
                if ((zza < za) || ((zza == za) && (zzb == 0))) {

                        // x = zz;
                        xa = zza;
                        xb = zzb;

                        // Here we use the precalculated values
                        // y += FXP_BKM_LOGS_L[shift];
                        auxa = FXP_BKM_LOGS[shift];
                        auxb = FXP_BKM_LOGS_XTRA[shift];
                        ya += auxa;
                        if (yb > FXP_BKM_MAXU - auxb)
                            // Carry from yb into ya
                            ya++;
                        yb += auxb;
                }
        }
        // Include any carry overs from summing the xtra bits together
        // Now y has the mantissa, shift it adjusting for final fxp
        int m;
        if (FXP_lg2_yashift >= 0)
                m = ya >> FXP_lg2_yashift;
        else {
                // Shifting the mantissa to the left, so we should pull
                // the left-most bits from yb
                m = (ya << (-FXP_lg2_yashift)) |
                        ((yb & FXP_lg2_ybmask) >> FXP_lg2_ybshift);
        }
        // Return the complete logarithm (c + m) as fxp
        if (c < FXP_whole_min) {
                if ((c + 1 == FXP_whole_min) && (m > 0)) {
                        // Special case: in spite of not enough whole
                        // bits available for c, there are for c + 1,
                        // which is what this logarithm will return as
                        // whole part
                        return INT_MIN + m;
                }
                // Overflow: the characteristic of the logarithm
                // will not fit in the whole bits part of our
                // current fxp settings
                return FXP_NEG_INF;
        }
        if (c >= 0)
            return (c << FXP_frac_bits) + m;
        else
            return -(-c << FXP_frac_bits) + m;
}

int fxp_pow2(int fxp1)
{
        // TODO:
        return 0;
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

