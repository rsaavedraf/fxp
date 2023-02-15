/*
 * fxp.c
 * An implementation of Binary Fixed Point numbers
 * encoding them into integers, together with their
 * arithmetic operations +, -, *, and /,
 * also log2, ln, and others.
 *
 * Safe arithmetic operations for the fxp's compliant
 * with INT32-C, as in:
 * https://wiki.sei.cmu.edu/confluence/display/c/INT32-C.+Ensure+that+operations+on+signed+integers+do+not+result+in+overflow
 *
 * By Raul Saavedra, Bonn, Germany
 *
 * v0.1: 2022-11-13
 * v0.2: 2023-01-08: runtime-modifiable number of frac bits to use.
 * v0.3: 2023-01-30: fxp_mul avoiding precision loss!
 * v0.4: 2023-02-09: first version of fxp_log2_l implemented
 */

#include "fxp.h"
#include "fxp_constants.h"
#include <stdio.h>
#include <stdlib.h>
//#include <assert.h>

//Used when testing and debugging when trying to optimize division
#include "fxp_aux.h"
//#define VERBOSE 1

#define FXP_FRAC_BITS_MIN 4
#define FXP_FRAC_BITS_DEF 16
// Allowing for no whole part, so other than the sign bit,
// all bits used for fraction part. This use case represents
// the range [-0.999..., 0.999...], or equivalently: (-1, 1)
// (so actual bits -1 and 1 not included)
#define FXP_FRAC_BITS_MAX FXP_INT_BITS_M1
#define FXP_FRAC_MAX_DEC 9999999

// Default number of bits to use for the frac part.
// Can be changed dynamically calling fxp_set_frac_bits()
static int fxp_frac_bits = FXP_FRAC_BITS_DEF;

// For improved-precision version of fxp_mul
static int fxp_frac_mshift = FXP_FRAC_BITS_DEF / 2;
static int fxp_frac_maskl = 4032;
static int fxp_frac_maskr = 63;

// Default number of bits for the whole (including sign) part
static int fxp_whole_bits = FXP_INT_BITS - FXP_FRAC_BITS_DEF;
static int fxp_whole_bits_m1 = FXP_INT_BITS - FXP_FRAC_BITS_DEF - 1;

// FXP_FRAC_MAX should correspond to 2^FXP_FRAC_BITS - 1
// Also used as mask for binary frac part of the int
static int fxp_frac_mask = ((1 << FXP_FRAC_BITS_DEF) - 1);
static int fxp_frac_max = ((1 << FXP_FRAC_BITS_DEF) - 1);

// Default desired max frac decimal value
// (can be changed dynamically calling set_[auto_]frac_max_dec
static int fxp_frac_max_dec = 9999;

static long int fxp_max_lshifted = (FXP_MAX_L) << FXP_FRAC_BITS_DEF;

// Max and min valid values for the whole part of the fxp's
static int fxp_whole_max = FXP_MAX >> FXP_FRAC_BITS_DEF;
static int fxp_whole_min = -(FXP_MAX >> FXP_FRAC_BITS_DEF);

/*
 * For the BKM log calculation when using longs.
 * Values in this array are: a[k] = log2(1 + 1/2^k) represented as
 * unsigned long fxp's (8 bytes,) with 62 frac bits
 */
static const unsigned long FXP_BKM_LOGS_L[] = {
        0x4000000000000000, 0x2570068E7EF5A27D, 0x149A784BCD1B8B50, 0xAE00D1CFDEB43FB,
        0x598FDBEB244C5B5, 0x2D75A6EB1DFB0F1, 0x16E79685C2D229E, 0xB7F285B778428E,
        0x5C2711B5EAB1DE, 0x2E1F07FE14EACA, 0x1712653743F454, 0xB89EB17BCABE1,
        0x5C523B0A86FF2, 0x2E29D623F4A6C, 0x1715193B17D35, 0xB8A982801725,
        0x5C54EF6A3E08, 0x2E2A833FB72C, 0x171544828311, 0xB8AA2F9EB95,
        0x5C551AB2053, 0x2E2A8E11ACC, 0x1715473700F, 0xB8AA3A70B1,
        0x5C551D6683, 0x2E2A8EBECC, 0x1715476248, 0xB8AA3B1DD,
        0x5C551D91C, 0x2E2A8EC99, 0x17154764F, 0xB8AA3B28,
        0x5C551D94, 0x2E2A8ECA, 0x17154765, 0xB8AA3B2,
//        0x5C551D9, 0x2E2A8EC, 0x1715476, 0xB8AA3B,      // <---- *
//        0x5C551D, 0x2E2A8E, 0x171547, 0xB8AA3,
//        0x5C551, 0x2E2A8, 0x17154, 0xB8AA,
//        0x5C55, 0x2E2A, 0x1715, 0xB8A,
//        0x5C5, 0x2E2, 0x171, 0xB8,
//        0x5C, 0x2E, 0x17, 0xB,
//        0x5, 0x2, 0x1
        // Starting with the row marked with the *, each entry is exactly
        // a 4 -bit right-shift of the value 4 positions earlier
};

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
#define FXP_DEF_SHIFT (FXP_INT_BITS_M1 - FXP_FRAC_BITS_DEF - 2)
static int fxp_e = FXP_E_I32 >> FXP_DEF_SHIFT;
static int fxp_pi = FXP_PI_I32 >> FXP_DEF_SHIFT;
static int fxp_lg2e = FXP_LOG2E_I32 >> (FXP_DEF_SHIFT + 1);
static int fxp_ln2 = FXP_LN2_I32 >> (FXP_DEF_SHIFT + 2);

#define FXP_BKM_PREC (FXP_INT_BITS - 2)
//#define FXP_BKM_MSHIFT (FXP_INT_BITS - 8)
//#define FXP_BKM_MMASK (0xF << FXP_BKM_MSHIFT)
#define FXP_BKM_PREC_L ((sizeof(long) * 8) - 2)
#define FXP_BKM_L_CLZSHIFT (FXP_BKM_PREC_L - FXP_BKM_PREC - 1)
static const int FXP_BKM_ONE = (1 << FXP_BKM_PREC);
static const unsigned long FXP_BKM_ONE_L = \
                            ((unsigned long) 1) << FXP_BKM_PREC_L;
static const unsigned int FXP_BKM_HALF = (FXP_BKM_ONE >> 1);
static const unsigned int FXP_BKM_MSB1 = (1 << (FXP_INT_BITS - 1));
static const unsigned int FXP_BKM_IMSB1 = ~FXP_BKM_MSB1;
static const unsigned int FXP_BKM_NMSB1 = (FXP_BKM_MSB1 >> 1);
static const unsigned int FXP_BKM_MAXU = ~((unsigned int) 0);
static const unsigned long FXP_BKM_HALF_L = (FXP_BKM_ONE_L >> 1);

// Auxiliary variables used in the lg2 implementations
static int fxp_half = 1 << (FXP_FRAC_BITS_DEF - 1);
static int fxp_one = 1 << FXP_FRAC_BITS_DEF;
static int fxp_two = 1 << (FXP_FRAC_BITS_DEF + 1);
static int fxp_lg2_l_maxloops = FXP_FRAC_BITS_DEF + 1;
static int fxp_lg2_l_shift = FXP_BKM_PREC_L - FXP_FRAC_BITS_DEF;
static int fxp_lg2_maxloops = FXP_FRAC_BITS_DEF;
static int fxp_lg2_yashift = FXP_BKM_PREC - FXP_FRAC_BITS_DEF;
static int fxp_lg2_ybmask = 0;
static int fxp_lg2_ybshift = 0;

int fxp_get_frac_bits() {
        return fxp_frac_bits;
}

int fxp_get_whole_bits()
{
        return fxp_whole_bits;
}

int fxp_get_frac_mask()
{
        return fxp_frac_mask;
}

int fxp_get_frac_max()
{
        return fxp_frac_max;
}

int fxp_get_whole_max()
{
        return fxp_whole_max;
}

int fxp_get_whole_min()
{
        return fxp_whole_min;
}

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
        fxp_frac_bits = (nfracbits < FXP_FRAC_BITS_MIN? FXP_FRAC_BITS_MIN:
                            (nfracbits > FXP_FRAC_BITS_MAX?
                                FXP_FRAC_BITS_MAX: nfracbits));

        fxp_whole_bits = FXP_INT_BITS - fxp_frac_bits;
        fxp_whole_bits_m1 = fxp_whole_bits - 1;

        // fxp_frac_mask should correspond to 2^FXP_FRAC_BITS - 1
        fxp_frac_mask = (1 << fxp_frac_bits) - 1;

        // Variables used in fxp_mul to process the
        // full multiplication of the frac parts avoiding
        // precision loss
        if (fxp_frac_bits == FXP_INT_BITS_M1) {
            fxp_frac_mshift = (FXP_INT_BITS_M1 - 1) / 2;
        } else {
            fxp_frac_mshift = (fxp_frac_bits / 2) + (fxp_frac_bits % 2);
        }
        fxp_frac_maskr = (1 << fxp_frac_mshift) - 1;
        fxp_frac_maskl = fxp_frac_maskr << fxp_frac_mshift;

        // When using all bits for frac (except for the sign bit),
        // then our max valid frac cannot be equal to frac_mask
        // (because in that case that value is already the largest
        // positive integer == POS_INF), so we must substract one
        // from the frac mask to get the largest valid frac value
        fxp_frac_max = fxp_frac_mask - (fxp_whole_bits == 1? 1: 0);

        fxp_max_lshifted = ((FXP_MAX_L) << fxp_frac_bits);

        // Max and min valid values for the whole part of the fxp's
        fxp_whole_max = FXP_MAX >> fxp_frac_bits;
        fxp_whole_min = (-fxp_whole_max);

        // Adjust precision of e, pi, etc. to the frac bits in use
        fxp_e = fxp_change_nfracbits(FXP_E_I32, \
                        FXP_INT_BITS - 3, fxp_frac_bits);
        fxp_pi = fxp_change_nfracbits(FXP_PI_I32, \
                        FXP_INT_BITS - 3, fxp_frac_bits);
        fxp_lg2e = fxp_change_nfracbits(FXP_LOG2E_I32, \
                        FXP_INT_BITS - 2, fxp_frac_bits);
        fxp_ln2 = fxp_change_nfracbits(FXP_LN2_I32, \
                        FXP_INT_BITS - 1, fxp_frac_bits);

        // Auxiliary variables used for lg2 calculations
        fxp_one = fxp(1);
        fxp_half = fxp_one >> 1;
        fxp_two = fxp_one << 1;
        fxp_lg2_l_maxloops = MIN(FXP_INT_BITS, fxp_frac_bits + 1);
        fxp_lg2_l_shift = FXP_BKM_PREC_L - fxp_frac_bits;
        fxp_lg2_maxloops = fxp_frac_bits + 1;
        fxp_lg2_yashift = FXP_BKM_PREC - fxp_frac_bits;
        if (fxp_lg2_yashift >= 0) {
            fxp_lg2_ybshift = 0;
            fxp_lg2_ybmask = 0;
        } else {
            // When lg2_shift is negative, lg2_ybmask is equal to
            // a number of 1's == -lg2_shift, shifted to the
            // left-most bits of an int
            fxp_lg2_ybshift = FXP_INT_BITS + fxp_lg2_yashift;
            fxp_lg2_ybmask = ((1 << (-fxp_lg2_yashift)) - 1) \
                    << fxp_lg2_ybshift;
        }

        return fxp_frac_bits;
}

/*
 * Automatically set the fractional max decimal to use.
 * The number of nines in fxp_frac_max_dec will be
 * floor(fxp_frac_bits / 4), e.g. one nine for every 4 bits
 */
int fxp_set_auto_frac_max_dec()
{
        int nnines = fxp_frac_bits / 4;
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
        if (whole > fxp_whole_max)
                return FXP_POS_INF;
        if (whole < fxp_whole_min)
                return FXP_NEG_INF;
        if (bin_frac > fxp_frac_max)
                bin_frac = fxp_frac_max;
        else if (bin_frac < -fxp_frac_max)
                bin_frac = -fxp_frac_max;
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
        int positive_fxp = (whole << fxp_frac_bits) | bin_frac;
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
        if (whole > fxp_whole_max)
                return FXP_POS_INF;
        if (whole < fxp_whole_min)
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
                                (long) fxp_frac_max) / \
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
                        // returning just what our whole-part bits of
                        // -INF (and also UNDEF) actually correspond to.
                        // However, see also the comment in
                        // fxp_get_bin_frac()
                        return -fxp_whole_max;
                else
                        return -((-fxp) >> fxp_frac_bits);
        else
                return (fxp >> fxp_frac_bits);
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
                        // and UNDEF. Here for -INF we return the actual
                        // bits used, but for UNDEF something different
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
                        return -fxp_frac_mask - (fxp == FXP_UNDEF? 1: 0);
                else
                        return -((-fxp) & fxp_frac_mask);
        else
                return (fxp & fxp_frac_mask);
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
        denom = (long) fxp_frac_max;
        int positive_frac;
        if (fxp < 0) {
                positive_frac = (-fxp) & fxp_frac_mask;
                if (positive_frac == 0) return 0;
                num = -(((long) positive_frac + 1) * \
                                ((long) fxp_frac_max_dec));
        }
        else {
                positive_frac = fxp & fxp_frac_mask;
                if (positive_frac == 0) return 0;
                num = ((long) positive_frac + 1)
                                * ((long) fxp_frac_max_dec);
        }
        ldivision = num / denom;
        if (ldivision > fxp_frac_max_dec) ldivision = fxp_frac_max_dec;
        int idivision = (int) ldivision;
        //printf("\n+frac is %d\n", positive_frac);
        //printf("Num   is %ld\n", num);
        //printf("Denom is %ld\n", denom);
        //printf("LDiv  is %ld\n", ldivision);
        //printf("iDiv  is %d\n", idivision);
        return idivision;
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
        return (int) (p >> fxp_frac_bits);
}

int fxp_unsafe_div(int fxp1, int fxp2)
{
        long n1 = ((long) fxp1) << fxp_frac_bits;
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
 * Safe implementation of fxp1 + fxp2 using longs.
 * Always ends up very slightly slower than fxp_add
 * at least on an intel i7, so there is no benefit
 * in using this at all over fxp_add. Commenting it out
 *
int fxp_add_l(int fxp1, int fxp2)
{
        // Check for undef or infinity arguments
        if (fxp1 == FXP_UNDEF || fxp2 == FXP_UNDEF)
                return FXP_UNDEF;
        if (fxp1 == FXP_POS_INF || fxp2 == FXP_POS_INF)
                return (fxp1 == FXP_NEG_INF || fxp2 == FXP_NEG_INF)?
                        FXP_UNDEF: FXP_POS_INF;
        if (fxp1 == FXP_NEG_INF || fxp2 == FXP_NEG_INF)
                return FXP_NEG_INF;

        long sum = ((long) fxp1) + fxp2;
        // Check for overflows
        if (sum > FXP_MAX) return FXP_POS_INF;
        if (sum < FXP_MIN) return FXP_NEG_INF;
        // No overflow, return the sum
        return ((int) sum);
}
*/

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
int fxp_sub_l(int fxp1, int fxp2)
{
        // Conservative check to never attempt to change sign
        // to the true most negative int (aka. FXP_UNDEF)
        if (fxp2 == FXP_UNDEF) return FXP_UNDEF;
        return fxp_add_l(fxp1, -fxp2);
}
*/

/*
 * Safe implementation of fxp multiplication using longs,
 * and no divisions.
 * Only applicable for systems in which sizeof(long) >= 2 * sizeof(int)
 */
int fxp_mul_l(int fxp1, int fxp2)
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
        int v1, v2;
        v1 = (fxp1 >= 0)? fxp1: -fxp1;
        v2 = (fxp2 >= 0)? fxp2: -fxp2;

        unsigned long product = ((unsigned long) v1) * v2;
        if (product > fxp_max_lshifted) {
                // Overflow, return infinity with the appropriate sign
                return ((fxp1 >= 0 && fxp2 >= 0) || \
                            (fxp1 < 0 && fxp2 < 0))?
                                FXP_POS_INF: FXP_NEG_INF;
        }
        // No overflow, return result as int with appropriate sign
        product = product >> fxp_frac_bits;
        return ((fxp1 >= 0 && fxp2 >= 0) || (fxp1 < 0 && fxp2 < 0))?
                        (int) product: -((int) product);
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
 * Notice that a*y and b*x can each never overflow the number of
 * bits in an int. Yet when frac_bits >= half the size of an int,
 * pf3 alone could overflow. But we process the calculation of pf3
 * itself with a distributive approach as well, not using whole vs.
 * frac parts, but left vs. right chunks.
 *
 * Notice that the final fraction part will be just the frac part
 * of the sum of the pfi's, because any carry over into the whole
 * part becomes pw4 (one of the components to build up the whole
 * part, see below).
 *
 * Notice also, when all bits (except sign one) are used for the
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
        if (nba + nbb > fxp_whole_bits) {
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
            pf3 = (x * y) >> fxp_frac_bits;
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
            int oddnfb = fxp_frac_bits % 2;
            if (fxp_frac_bits == FXP_INT_BITS_M1) {
                // Eliminate the very lowest significant bit of each
                // operand to make sure none of the subchunk products
                // ever exceed FXP_INT_BITS_M1 bits.
                // (We must shift back this lost bit in the final pf3)
                x = x >> 1;
                y = y >> 1;
                oddnfb = -1;
            } else {
                x = x << oddnfb;
                y = y << oddnfb;
            }
            // Variables fxp_frac_maskl, _maskr, and _mshift are
            // all preinitialized for the default # of frac bits to use,
            // and readjusted when fxp_set_frac_bits() gets called
            int xl = (x & fxp_frac_maskl) >> fxp_frac_mshift;
            int xr = (x & fxp_frac_maskr);
            int yl = (y & fxp_frac_maskl) >> fxp_frac_mshift;
            int yr = (y & fxp_frac_maskr);
            int xlyr = xl * yr;
            int ylxr = yl * xr;
            int qr1 = xlyr & fxp_frac_maskr;
            int qr2 = ylxr & fxp_frac_maskr;
            int qr3 = (xr * yr) >> fxp_frac_mshift;
            int qrsum = qr1 + qr2 + qr3;
            int ql1 = xl * yl;
            int ql2 = (xlyr & fxp_frac_maskl) >> fxp_frac_mshift;
            int ql3 = (ylxr & fxp_frac_maskl) >> fxp_frac_mshift;
            int ql4 = (qrsum & fxp_frac_maskl) >> fxp_frac_mshift;
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
        if (pwsum > fxp_whole_max) {
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

        if ((y <= fxp_frac_mask) &&
                (x > fxp_mul(FXP_MAX, y))) {
                return ((fxp1 > 0 && fxp2 > 0) || (fxp1 < 0 && fxp2 < 0))?
                                FXP_POS_INF: FXP_NEG_INF;
        }

        int bx = fxp_nbits(x);
        int by = fxp_nbits(y);

        int m, next_bit_mask;
        int bmax = bx + fxp_frac_bits;
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
        //printf("bidx processed while <= bmax: %d\n", bmax);
        //int loops = 0, difference = 0, ba;
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

/*
 * Safe implementation of fxp division using longs.
 * Only applicable for systems in which sizeof(long) >= 2 * sizeof(int)
 */
int fxp_div_l(int fxp1, int fxp2)
{
        if (fxp1 == FXP_UNDEF || fxp2 == FXP_UNDEF)
                return FXP_UNDEF;
        if (fxp2 == 0)
                return (fxp1 > 0)? FXP_POS_INF:
                                (fxp1 < 0)? FXP_NEG_INF: FXP_UNDEF;
        // Positive values of the arguments
        int x = (fxp1 >= 0)? fxp1: -fxp1;
        int y = (fxp2 >  0)? fxp2: -fxp2;
        if (y == FXP_POS_INF)
                return (x == FXP_POS_INF)? FXP_UNDEF: 0;
        if (x == FXP_POS_INF)
                return ((fxp1 > 0 && fxp2 > 0) || (fxp1 < 0 && fxp2 < 0))?
                                FXP_POS_INF: FXP_NEG_INF;
        // Compute positive division
        long numerator = ((long) x) << fxp_frac_bits;
        long division = numerator / y;

        if (division > FXP_MAX_L) {
                // Overflow -> Return properly signed infinity
                return ((fxp1 >= 0 && fxp2 > 0) || (fxp1 < 0 && fxp2 < 0))?
                                FXP_POS_INF: FXP_NEG_INF;
        }
        // No overflow -> return properly signed int
        return ((fxp1 >= 0 && fxp2 > 0) || (fxp1 < 0 && fxp2 < 0))?
                        (int) division: -((int) division);
}

int fxp_get_e()
{
        return fxp_e;
}

int fxp_get_pi()
{
        return fxp_pi;
}


int fxp_get_ln2()
{
        return fxp_ln2;
}

int fxp_get_lg2e()
{
        return fxp_lg2e;
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

/*
 * Default implementation of ln()
 * Uses the default log2() and only ints
 */
int fxp_ln(int fxp1)
{
        return fxp_mul(fxp_lg2(fxp1), fxp_ln2);
}

/*
 * Implementation of ln() using longs
 */
int fxp_ln_l(int fxp1)
{
        return fxp_mul_l(fxp_lg2_l(fxp1), fxp_ln2);
}

/*
 * Alternative lg2 of an fxp using multiplication.
 * Requires the fxp configuration to have 3 or more whole bits.
 * Needs no pre-calculated table of log values in any range, but
 * requires one multiplication per mantissa bit.
 *
 * Adapted to fxps from the general algorithm to calculate
 * binary logarithms explained by Clay Turner in IEEE Signal
 * Processing Magazine, Sep/2010.
 * D. E. Knuth's "The Art of Computer Programming Vol 2:
 * Seminumerical Algorithms", 2nd ed., 1981 (pages 441 - 446) is
 * the one and only reference in that short article, so that's the
 * effective reference.
 *
 */
int fxp_lg2_mul(int fxp1)
{
        if (fxp1 <= 0) return ((fxp1 == 0)? FXP_NEG_INF: FXP_UNDEF);
        if (fxp1 == FXP_POS_INF) return FXP_POS_INF;
        // A logarithm value consists of a characteristic (whole part)
        // + a mantissa (fraction part).
        // The mantissa calculation done by this algorithm assumes
        // 1 <= argument < 2, so we need to get there first. We do that
        // by calculating a normalized version of the original argument
        // with a succession of divides/multiplies by 2 (that is,
        // r- or l-shifts).
        // The count of those operations will actually correspond to
        // the log characteristic, with divides counted positively,
        // and multiplies negatively.
        int c = 0; // characteristic
        int nx;
        int nbx = fxp_nbits(fxp1);
        if (fxp1 < fxp_one) {
                c = nbx - fxp_frac_bits - 1;
                nx = fxp1 << (-c);
        } else if (fxp1 >= fxp_two) {
                c = nbx - fxp_frac_bits - 1;
                nx = fxp1 >> c; // notice: inevitably losing frac bits here
        } else {
                nx = fxp1;
        }
        // Mantissa calculation:
        // Here we have already calculated the log characteristic c, and
        // we have nx satisfying: 1 <= nx < 2
        int m = 0; // starting mantissa
        // b starts as 0.5 (corresponding to first frac bit)
        int b = fxp_half;
        while (b > 0) {
                // Here comes a squaring of nx, so one inevitable
                // multiplication required per bit of the mantissa. This
                // is by far the most expensive part of this entire log2
                // algorithm
                nx = fxp_mul(nx, nx);
                if (nx >= fxp_two) {
                        nx = nx >> 1;
                        // Notice the following sum is safe, cannot
                        // possibly ever overflow because m always
                        // remains < fxp_one
                        // Sets to 1 mantissa bit corresponding to b
                        m |= b;
                }
                b = b >> 1;
        }
        //int fxp_log = (c << fxp_frac_bits) | m;
        //printf("\nlog(%d) = c + m = %d\n", fxp1, fxp_log);
        // Return the calculated logarithm as fxp
        if (c < fxp_whole_min) {
                if ((c + 1 == fxp_whole_min) && (m > 0)) {
                    return INT_MIN + m;
                }
                return FXP_NEG_INF;
        }
        if (c >= 0)
            return (c << fxp_frac_bits) + m;
        else
            return -(-c << fxp_frac_bits) + m;
}

/*
 * Same lg2 implementation as above using multiplication, but
 * here a long multiplication. Only applicable
 * for systems in which sizeof(long) >= 2 * sizeof(int)
 * Requires the fxp configuration to have 3 or more whole bits.
 * Needs no pre-calculated table of log values in any range.
 */
int fxp_lg2_mul_l(int fxp1)
{
        if (fxp1 <= 0) return ((fxp1 == 0)? FXP_NEG_INF: FXP_UNDEF);
        if (fxp1 == FXP_POS_INF) return FXP_POS_INF;
        int c = 0; // characteristic
        int nx;
        int nbx = fxp_nbits(fxp1);
        if (fxp1 < fxp_one) {
                c = nbx - fxp_frac_bits - 1;
                nx = fxp1 << (-c);
        } else if (fxp1 >= fxp_two) {
                c = nbx - fxp_frac_bits - 1;
                nx = fxp1 >> c;
        } else {
                nx = fxp1;
        }
        // Mantissa calculation:
        int m = 0;
        int b = fxp_half;
        while (b > 0) {
                nx = fxp_mul_l(nx, nx);
                if (nx >= fxp_two) {
                        nx = nx >> 1;
                        m |= b;
                }
                b = b >> 1;
        }
        // Return the calculated logarithm as fxp
        if (c < fxp_whole_min) {
                if ((c + 1 == fxp_whole_min) && (m > 0)) {
                    return INT_MIN + m;
                }
                return FXP_NEG_INF;
        }
        if (c >= 0)
            return (c << fxp_frac_bits) + m;
        else
            return -(-c << fxp_frac_bits) + m;
}

/*
 * Default log2 calculation using only ints and the BKM algorithm.
 * Requires tables of pre-calculated values (FXP_BKM_LOGS[], and
 * FXP_BKM_LOGS_XTRA[]
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
        int c = ((fxp1 < fxp_one) || (fxp1 >= fxp_two))?
                    (nbx - fxp_frac_bits - 1): 0;
        //printf("fxp:%x, clz:%d, nbx:%d, c:%d\n", fxp1, clz, nbx, c);

        // Here simulating longs with two ints a and b, i.e.
        // a long K as two ints (ka, kb), ka being the most
        // significant int, and kb the least significant int
        unsigned int za = ((unsigned int) fxp1) << (clz - 1);
        //unsigned int zb = 0; // <- zb not really needed, will remain 0

        //printf("lg2   start:  c:%d,  fxp:%x,  z:%8X%8X\n", \
        //        c, fxp1, za, 0);
                                                    // Simulated long
        unsigned int xa = FXP_BKM_ONE, xb = 0;      // x
        unsigned int zza, zzb;                      // zz
        unsigned int xsa = 0, xsb = 0;              // xs
        unsigned int ya = 0, yb = 0;                // y
        unsigned int auxa, auxb;                    // aux
        unsigned int a_bits;
        unsigned int ab_mask = 1;
        for (int shift=1; shift <= fxp_lg2_maxloops; shift++) {

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

                //printf("lg2 %2d: x :%8X-%8X\n", shift, xa, xb);
                //printf("        xs:%8X-%8X\n", xsa, xsb);
                //printf("        zz:%8X-%8X\n", zza, zzb);

                // if (zz <= z) {
                if ((zza < za) || ((zza == za) && (zzb == 0))) {
                        //printf("\tUPDATING X AND M.\n");

                        // x = zz;
                        xa = zza;
                        xb = zzb;

                        // y += FXP_BKM_LOGS_L[shift];
                        auxa = FXP_BKM_LOGS[shift];
                        auxb = FXP_BKM_LOGS_XTRA[shift];
                        ya += auxa;
                        if (yb > FXP_BKM_MAXU - auxb)
                            // Carry from yb into ya
                            ya++;
                        yb += auxb;
                        //printf("\told_m_b:%u,  aux:%u,  new m_b:%u\n", \

                        //    old_m_b, aux, m_b);
                }
                //printf("lg2   %2d:  zz:%8X%8X  x:%8X-%8X  xs:%8X-%8X  m:%8X%8X  %d\n",
                //        shift, zz_a, zz_b, \
                //        x_a, x_b, xs_a, xs_b, \
                //        m_a, m_b, shift);
                //if (shift == 12) exit(-3);
        }
        // Include any carry overs from summing the xtra bits together
        // Now y has the mantissa, shift it adjusting for final fxp
        int m;
        if (fxp_lg2_yashift >= 0)
                m = ya >> fxp_lg2_yashift;
        else {
                // If shifted the mantissa to the left, we can pull the
                // left-most bits from yb
                m = (ya << (-fxp_lg2_yashift)) |
                        ((yb & fxp_lg2_ybmask) >> fxp_lg2_ybshift);
                //printf("ma:%x,  mb:%x ---> m:%x\n", ya, yb, m);
        }
        //printf("Final mantissa after %d shift from lg2: %x\n", fxp_lg2_shift, m);
        // Return the complete logarithm (c + m) as fxp
        if (c < fxp_whole_min) {
                if ((c + 1 == fxp_whole_min) && (m > 0)) {
                        // Special case: in spite of not enough whole
                        // bits available for c, there are for c + 1,
                        // which is what this logarithm will return as
                        // whole part
                        return INT_MIN + m;
                }
                // Overflow: the characteristic of the logarithm
                // will not fit in the whole bits part of our
                // current fxp settings
                //printf("\n\nReturning -INF here: c:%d, whole_min:%d, m:%d!!!!!\n\n", c, fxp_whole_min, m);
                return FXP_NEG_INF;
        }
        if (c >= 0)
            return (c << fxp_frac_bits) + m;
        else
            return -(-c << fxp_frac_bits) + m;
}

/*
 * log2 calculation using longs and the BKM algorithm
 */
int fxp_lg2_l(int fxp1)
{
        if (fxp1 <= 0) return ((fxp1 == 0)? FXP_NEG_INF: FXP_UNDEF);
        if (fxp1 == FXP_POS_INF) return FXP_POS_INF;
        // Here fxp1 for sure > 0
        int clz = __builtin_clz((unsigned int) fxp1);
        // clz > 0 since at least sign bit == 0
        int nbx = FXP_INT_BITS - clz;
        // c is characteristic
        int c = ((fxp1 < fxp_one) || (fxp1 >= fxp_two))?
                    (nbx - fxp_frac_bits - 1): 0;
        //printf("fxp:%x, clz:%d, nbx:%d, c:%d\n", fxp1, clz, nbx, c);

        //printf("Proceeding with BKM loop for mantissa:\n");
        // Here we have already calculated the log characteristic c
        // Now lshifting the original number to have it as unsigned
        // long with 61 frac bits (same as the pre-calculated values
        // in our bkm table of logs)
        unsigned long z = ((unsigned long) fxp1) << \
                            (clz + FXP_BKM_L_CLZSHIFT);
        //printf("lg2_l start:  c:%d,  fpx:%x,  z:%lx\n", c, fxp1, z);
        unsigned long x = FXP_BKM_ONE_L;
        unsigned long zz;
        unsigned long y = 0; // starting mantissa
        for (int shift=1; shift <= fxp_lg2_l_maxloops; shift++) {
                unsigned long xs = (x >> shift);
                zz = x + xs;
                //printf("lg2_l %2d: x :%16lX\n", shift, x);
                //printf("          xs:%16lX\n", xs);
                //printf("          zz:%16lX\n", zz);
                if (zz <= z) {
                        unsigned long aux = FXP_BKM_LOGS_L[shift];
                        //printf("\tUPDATING X AND Y, array: %lX\n", aux);
                        x = zz;
                        y += aux;
                }
                //printf("lg2_l %2d:  zz:%16lX  x:%16lX  xs:%16lX  m:%16lX  %d\n",
                //        shift, zz, x, xs, y, shift);
        }
        // Now y has the mantissa, shift it adjusting for final fxp
        long m;
        if (fxp_lg2_l_shift >= 0) {
                m = y >> fxp_lg2_l_shift;
        } else {
                m = y << (-fxp_lg2_l_shift);
        }
        int im = (int) m;
        //printf("Final mantissa from lg2_l: %x\n", im);
        // Return the complete logarithm (c + m) as fxp
        if (c < fxp_whole_min) {
                if ((c + 1 == fxp_whole_min) && (im > 0)) {
                        // Special case: in spite of not enough whole
                        // bits available for c, there are for c + 1,
                        // which is what this logarithm will return as
                        // whole part
                        return INT_MIN + im;
                }
                // Overflow: the characteristic of the logarithm
                // will not fit in the whole bits part of our
                // current fxp settings
                //printf("\n\nReturning -INF here: c:%d, whole_min:%d, m:%d!!!!!\n\n", c, fxp_whole_min, m);
                return FXP_NEG_INF;
        }
        if (c >= 0)
            return (c << fxp_frac_bits) + im;
        else
            return -(-c << fxp_frac_bits) + im;
}
