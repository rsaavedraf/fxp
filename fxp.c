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
//#include "fxp_aux.h"
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

// fxp rounded representations of transcendental constants using
// the default number of frac bits
#define FXP_DEF_SHIFT (FXP_INT_BITS_M1 - FXP_FRAC_BITS_DEF - 2)
static int fxp_e = FXP_E_I32 >> FXP_DEF_SHIFT;
static int fxp_pi = FXP_PI_I32 >> FXP_DEF_SHIFT;
static int fxp_log2e = FXP_LOG2E_I32 >> (FXP_DEF_SHIFT + 1);
static int fxp_ln2 = FXP_LN2_I32 >> (FXP_DEF_SHIFT + 2);

// Auxiliary variables used for log2 calculations
static int fxp_half = 1 << FXP_FRAC_BITS_DEF - 1;
static int fxp_one = 1 << FXP_FRAC_BITS_DEF;
static int fxp_two = 1 << FXP_FRAC_BITS_DEF + 1;

// Array of 31 pre-computed log values to calculate
// logs base 2 using the BKM algorithm. The values here are
// simply: a[k] = log2(1 + 1/2^k) represented as fxp's using
// 30 frac bits (not rounded.)
static const unsigned int FXP_BKM_LOGS[] = {
        0x40000000, 0x2570068e, 0x149a784b, 0xae00d1c,  0x598fdbe,
        0x2d75a6e,  0x16e7968,  0xb7f285,   0x5c2711,   0x2e1f07,
        0x171265,   0xb89eb,    0x5c523,    0x2e29d,    0x17151,
        0xb8a9,     0x5c54,     0x2e2a,     0x1715,     0xb8a,
        0x5c5,      0x2e2,      0x171,      0xb8,       0x5c,
        0x2e,       0x17,       0xb,        0x5,        0x2,
        0x1,        0x0
};

// 24 extra bits following each of the frac values
// above, stored here in lowest 3 bytes of ints.
// Exploring higher precision calculation of log2's with BKM,
// while still only using ints
static const unsigned int FXP_BKM_LOGSXTRA[] = {
    0x000000,   0x7ef5a2,   0xcd1b8b,   0xfdeb43,   0xb244c5,
    0xb1dfb0,   0x5c2d22,   0xb77842,   0xb5eab1,   0xfe14ea,
    0x3743f4,   0x17bcab,   0xb0a86f,   0x623f4a,   0x93b17d,
    0x828017,   0xef6a3e,   0x833fb7,   0x448283,   0xa2f9eb,
    0x51ab20,   0xa8e11a,   0x547370,   0xaa3a70,   0x551d66,
    0x2a8ebe,   0x154762,   0x8aa3b1,   0xc551d9,   0xe2a8ec,
    0x715476,   0xb8aa3b,
};

static const int FXP_BKM_LOGSX_SHIFT = (FXP_INT_BITS - 8);
static const unsigned int FXP_BKM_LOGSX_MASK = (255 << FXP_BKM_LOGSX_SHIFT);
static const int FXP_BKM_PREC = FXP_INT_BITS - 2;
static const int FXP_BKM_LOOPS = FXP_BKM_PREC + 1;
static const unsigned int FXP_BKM_ONE = (1 << FXP_BKM_PREC);
static const unsigned int FXP_BKM_HALF = (FXP_BKM_ONE >> 1);

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

        // Adjust precision of e, pi, and log2e to the frac bits in use
        fxp_e = fxp_change_nfracbits(FXP_E_I32, FXP_INT_BITS - 3, fxp_frac_bits);
        fxp_pi = fxp_change_nfracbits(FXP_PI_I32, FXP_INT_BITS - 3, fxp_frac_bits);
        fxp_log2e = fxp_change_nfracbits(FXP_LOG2E_I32, FXP_INT_BITS - 2, fxp_frac_bits);
        fxp_ln2 = fxp_change_nfracbits(FXP_LN2_I32, FXP_INT_BITS - 1, fxp_frac_bits);

        // Auxiliary variables used for log calculation
        fxp_one = fxp(1);
        fxp_half = fxp_one >> 1;
        fxp_two = fxp_one << 1;

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
 * bin_frac should be a value between 0 and fxp_frac_max
 */
int fxp_bin(int whole, int bin_frac)
{
        if (whole > fxp_whole_max)
                return FXP_POS_INF;
        if (whole < fxp_whole_min)
                return FXP_NEG_INF;
        int sign = 1;
        if ((whole == 0) && (bin_frac < 0)) {
                // Special case for negative numbers when whole part is zero,
                // then the fxp gets its sign from the frac
                sign = -1;
                bin_frac = -bin_frac;
        } else {
                // All other cases, fxp gets its sign from the whole part
                if (whole < 0) {
                        sign = -1;
                        whole = -whole;
                }
                if (bin_frac < 0) {
                        bin_frac = -bin_frac;
                }
        }
        if (bin_frac > fxp_frac_max) bin_frac = fxp_frac_max;
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
        int bin_frac = (int) (((long) trunc_frac * (long) fxp_frac_max) / \
                            (long) (fxp_frac_max_dec));

        return fxp_bin(whole, ((frac_sign == 1)? bin_frac: -bin_frac));
}

int fxp_get_whole_part(int fxp)
{
        if (fxp < 0)
                if (fxp <= FXP_NEG_INF)
                        // Technically, the whole part of -INF would still be -INF,
                        // and both the whole and the frac parts of UNDEF would still
                        // be UNDEF. But here returning just what our whole-part bits
                        // of -INF (and also UNDEF) actually correspond to. However,
                        // see also the comment in fxp_get_bin_frac()
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
                        // In the whole-part function above we are returning the same
                        // whole parts for both -INF and UNDEF. Here for -INF we return
                        // the actual bits used, but for UNDEF something different
                        // and odd in order to enable differenciating -INF from UNDEF when
                        // checking their whole and frac constituents. It must be something
                        // that no other fxp could ever have, not even the +/-INF values.
                        // That can be some (alleged) frac bits that are actually
                        // invalid, e.g. outside the range of valid frac bits, like
                        // -frac_mask - 1. Notice that this value "overflows" the frac bits,
                        // but in a way makes sense to mark and differenciate precisely
                        // only UNDEF this way. Conveniently, when the # of frac bits is 31,
                        // this scheme returns the full FXP_UNDEF as the alleged frac part.
                        return -fxp_frac_mask - (fxp == FXP_UNDEF? 1: 0);
                else
                        return -((-fxp) & fxp_frac_mask);
        else
                return (fxp & fxp_frac_mask);
}

/*
 * Get the frac part as decimal between 0 and fxp_frac_max_dec,
e.g. 0 .. 999
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
                // Sum of operands would overflow
                // Return infinity of the appropriate sign
                return (fxp2 > 0? FXP_POS_INF: FXP_NEG_INF);
        // No overflow danger, do sum
        return fxp1 + fxp2;
}

/*
 * Safe implementation of fxp1 + fxp2 using longs.
 * Only applicable for systems in which sizeof(long) > sizeof(int)
 */
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

int fxp_sub_l(int fxp1, int fxp2)
{
        // Conservative check to never attempt to change sign
        // to the true most negative int (aka. FXP_UNDEF)
        if (fxp2 == FXP_UNDEF) return FXP_UNDEF;
        return fxp_add_l(fxp1, -fxp2);
}

/*
 * Safe implementation of fxp multiplication using longs,
 * and no divisions.
 * Only applicable for systems in which sizeof(long) > sizeof(int)
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
                return ((fxp1 >= 0 && fxp2 >= 0) || (fxp1 < 0 && fxp2 < 0))?
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
 * Only applicable for systems in which sizeof(long) > sizeof(int)
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

int fxp_get_log2e()
{
        return fxp_log2e;
}


int fxp_exp(int fxp1)
{
        return 0;
}

/*
 * Square root implementation
 */
int fxp_sqrt(int fxp1)
{
        return 0;
}

/*
 * Implementation of ln() using longs (for the multiplication)
 */
int fxp_ln_l(int fxp1)
{
        return fxp_mul_l(fxp_log2(fxp1), fxp_ln2);
}

/*
 * Default implementation of ln()
 * Uses the default log2() and int-only multiplication
 */
int fxp_ln(int fxp1)
{
        return fxp_mul(fxp_log2(fxp1), fxp_ln2);
}

/*
 * Alternative log2 of an fxp using multiplication.
 * Requires the fxp configuration to have 3 or more whole bits.
 * Needs no pre-calculated table of log values in any range, but
 * requires one multiplication per mantissa bit, so slower than
 * the default. Operates within the current fxp's settings,
 * so # effective mantissa bits calculated == current # frac bits,
 * which compensates a bit the execution time for fewer
 * frac bit configurations.
 *
 * Adapted to fxps from a general algorithm to calculate
 * binary logarithms explained by Clay Turner in IEEE Signal
 * Processing Magazine, Sep/2010.
 * D. E. Knuth's "The Art of Computer Programming Vol 2:
 * Seminumerical Algorithms", 2nd ed., 1981 (pages 441 - 446) is
 * the one and only reference in that short article, so that's the
 * effective reference.
 *
 */
int fxp_log2_mul(int fxp1)
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
        int b = fxp_half; // b starts as 0.5 (corresponding to first frac bit)
        while (b > 0) {
                // Here comes a squaring of nx, so one inevitable multiplication
                // required per bit of the mantissa. This is by far the most
                // expensive part of this entire log2 algorithm
                nx = fxp_mul(nx, nx);
                if (nx >= fxp_two) {
                        nx = nx >> 1;
                        // Notice the following sum is safe, cannot
                        // possibly ever overflow because m always
                        // remains < fxp_one
                        m += b; // Sets to 1 mantissa bit corresponding to b
                }
                b = b >> 1;
        }
        //int fxp_log = (c << fxp_frac_bits) | m;
        //printf("\nlog(%d) = c + m = %d\n", fxp1, fxp_log);
        // Return the calculated logarithm as fxp

        //return ((c << fxp_frac_bits) | m);

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
 * Same log2 implementation as above using multiplication, but
 * here a long multiplication for speed. Therefore, only applicable
 * for systems in which sizeof(long) > sizeof(int)
 * Requires the fxp configuration to have 3 or more whole bits.
 */
int fxp_log2_mul_l(int fxp1)
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
                        m += b;
                }
                b = b >> 1;
        }
        //return ((c << fxp_frac_bits) | m);

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
 * Default log2 calculation using the BKM algorithm, and
 * only ints (no longs)
 * For more details on BKM:
 * https://en.wikipedia.org/wiki/BKM_algorithm
 */
int fxp_log2(int fxp1)
{
        if (fxp1 <= 0) return ((fxp1 == 0)? FXP_NEG_INF: FXP_UNDEF);
        if (fxp1 == FXP_POS_INF) return FXP_POS_INF;
        // Here fxp1 for sure > 0
        int clz = __builtin_clz(fxp1);
        // clz > 0 since at least sign bit == 0
        int nbx = FXP_INT_BITS - clz;
        // c is characteristic
        int c = ((fxp1 < fxp_one) || (fxp1 >= fxp_two))?
                    (nbx - fxp_frac_bits - 1): 0;
        //printf("clz:%d, nbx:%d\n", clz, nbx);
        //printf("c:%d  Proceeding with BKM loop for mantissa:\n", c);
        // Here we have already calculated the log characteristic c
        // Now lshifting the original number to have it as unsigned
        // fxp with 30 frac bits (same as the pre-calculated values
        // in our bkm table of logs)
        unsigned int z = ((unsigned int) fxp1) << (clz - 1);
        unsigned int x = FXP_BKM_ONE;
        //unsigned int s = FXP_BKM_HALF; // b starts as 0.5
        unsigned int zz;
        unsigned int y = 0; // starting mantissa
        unsigned int yxtra = 0; // sum of xtra bits
        for (int shift=1; shift <= FXP_BKM_LOOPS; shift++) {
                unsigned int xs = (x >> shift);
                zz = x + xs;
                //printf("z-zz:%d,  z:%u,  zz:%u,  x:%u,  xs:%u\n",
                //        ((int) z) - ((int) zz), z, zz, x, xs);
                if (zz <= z) {
                        x = zz;
                        y += FXP_BKM_LOGS[shift];
                        yxtra += FXP_BKM_LOGSXTRA[shift];
                        //printf("\t\ty:%u  yxtra:%d\n", y, yxtra);
                }
        }
        unsigned int co = ((yxtra & FXP_BKM_LOGSX_MASK) >> \
                                FXP_BKM_LOGSX_SHIFT);
        //printf("\nyxtra:%u  x%X  co:%u  x%X\n", yxtra, yxtra, co, co);
        //printf("y before:%d\ty after:%d\n\n", y, y + yxtra);
        y += co;
        // Now y has the mantissa with 30 frac bits
        // Shift it adjusting it for final fxp
        int m;
        int shift = FXP_BKM_PREC - fxp_frac_bits;
        if (shift >= 0) {
                m = ((int) y) >> shift;
                if (shift > 0) {
                    // Round last bit in m
                    shift--;
                    int roundbit = (y & (1 << shift)) >> shift;
                    //printf("shift:%d  roundbit:%d  m before:%d  m after rounding:%d\n", \
                    //        shift, roundbit, m, m + roundbit);
                    //m = m + roundbit;
                }
        } else {
                m = ((int) y) << (-shift);
        }
        // Return the complete logarithm (c + m) as fxp
        if (c < fxp_whole_min) {
                if ((c + 1 == fxp_whole_min) && (m > 0)) {
                    // Special case: in spite of not enough whole bits
                    // available for c, there are for c + 1, which is what
                    // this logarithm will have as whole part
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
