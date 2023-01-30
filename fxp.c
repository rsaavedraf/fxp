/*
 * fxp.c
 * An implementation of Binary Fixed Point numbers
 * encoding them into integers, together with
 * arithmetic operations +, -, *, and / for them.
 *
 * Safe arithmetic operations for the fxp's compliant
 * with INT32-C, as in:
 * https://wiki.sei.cmu.edu/confluence/display/c/INT32-C.+Ensure+that+operations+on+signed+integers+do+not+result+in+overflow
 *
 * By Raul Saavedra, Bonn, Germany
 *
 * v1: 2022-11-13
 * v2: 2023-01-08: runtime-modifiable number of frac bits to use.
 */

#include "fxp.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

//Used when testing and debugging when trying to optimize division
//#include "fxp_aux.h"
//#define VERBOSE 0

#define FXP_FRAC_BITS_MIN 4
#define FXP_FRAC_BITS_DEF 12
// Allowing for no whole part, so other than the sign bit,
// all bits used for fraction part. This use case represents
// the range [-0.999..., 0.999...], or equivalently: (-1, 1)
// (so actual bits -1 and 1 not included)
#define FXP_FRAC_BITS_MAX FXP_INT_BITS_M1
#define FXP_FRAC_MAX_DEC 9999999

// Default number of bits to use for the frac part
// (can be changed dynamically calling set_frac_bits)
static int fxp_frac_bits = FXP_FRAC_BITS_DEF;

// Default number of bits for the whole (and sign) part
static int fxp_whole_bits = FXP_INT_BITS - FXP_FRAC_BITS_DEF;
static int fxp_whole_bits_m1 = FXP_INT_BITS - FXP_FRAC_BITS_DEF - 1;

// FXP_FRAC_MAX should correspond to 2^FXP_FRAC_BITS - 1
// Also used as mask for binary frac part of the int
static int fxp_frac_mask = ((1 << FXP_FRAC_BITS_DEF) - 1);
static int fxp_frac_max = ((1 << FXP_FRAC_BITS_DEF) - 1);

// Default desired max frac decimal value
// (can be changed dynamically calling set_[auto_]frac_max_dec
static int fxp_frac_max_dec = 999;

//#define FXP_MAX_LSHIFTED ((FXP_MAX_L) << FXP_FRAC_BITS)
static long int fxp_max_lshifted = (FXP_MAX_L) << FXP_FRAC_BITS_DEF;

// Max and min valid values for the whole part of the fxp's
static int fxp_whole_max = FXP_MAX >> FXP_FRAC_BITS_DEF;

//#define FXP_WHOLE_MIN (-FXP_WHOLE_MAX)
static int fxp_whole_min = -(FXP_MAX >> FXP_FRAC_BITS_DEF);


int fxp_get_frac_bits()
{
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

        return fxp_frac_bits;
}

/*
 * Automatically set the fractional max decimal to use.
 * It's value will be the integer containing only nines
 * that is > 0 and < fxp_frac_max.
 * I.e., for fxp_frac_max = 4095 (frac with 12 bits,)
 * the fxp_frac_max_dec value will set to 999.
 * The number of nines in fxp_frac_max_dec will be
 * floor(fxp_frac_bits / 4), e.g. a nine for every 4 bits
 */
int fxp_set_auto_frac_max_dec()
{
        /*
        int current = 9;
        int next = 99;
        int maxnext = (FXP_MAX - 9) / 10;
        while (next <= fxp_frac_mask) {
                current = next;
                next = (next * 10) + 9;
                if (next >= maxnext) break;
        }
        // Too large a value here might overflow in the bin-to-dec
        // conversion of fracs if using there only ints and not longs
        // Just in case limit here up to FXP_FRAC_MAX_DEC
        fxp_frac_max_dec = (current > FXP_FRAC_MAX_DEC?
                                FXP_FRAC_MAX_DEC: current);
        */
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
 * bin_frac should be a value between 0 and FXP_MAX_FRAC
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
 * FXP_FRAC_MAX_DEC, e.g. between 0 and 999
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
                if (fxp == FXP_UNDEF)
                        return fxp_whole_min;
                else
                        return -((-fxp) >> fxp_frac_bits);
        else
                return fxp >> fxp_frac_bits;
}

/*
 * Get the frac part directly (binary)
 */
int fxp_get_bin_frac(int fxp)
{
        if (fxp < 0)
                if (fxp == FXP_UNDEF)
                        return -fxp_frac_mask;
                else
                        return -((-fxp) & fxp_frac_mask);
        else
                return (fxp & fxp_frac_mask);
}

/*
 * Get the frac part as decimal between 0 and FXP_FRAC_MAX_DEC,
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
                //return -(((-fxp) & fxp_frac_mask) * fxp_frac_max_dec) \
                //            / fxp_frac_max;
                positive_frac = (-fxp) & fxp_frac_mask;
                if (positive_frac == 0) return 0;
                num = -(((long) positive_frac + 1) * \
                                ((long) fxp_frac_max_dec));
        }
        else {
                //return ((fxp & fxp_frac_mask) * fxp_frac_max_dec) \
                //            / fxp_frac_max;
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

        /*
        // Equivalent to the code above, but checking first
        // if arguments have the same signs or not
        if ((fxp1 >= 0 && fxp2 >= 0) || (fxp1 < 0) && (fxp2 < 0)) {
                // Same signs, check for possible overflow
                if (fxp1 > 0)
                        return (FXP_MAX - fxp1 < fxp2)?
                                    FXP_POS_INF: fxp1 + fxp2;
                else
                        return (FXP_MAX + fxp1 < -fxp2)?
                                    FXP_NEG_INF: fxp1 + fxp2;
        };
        // Arguments with different sign -> No overflow danger, do sum
        return fxp1 + fxp2;
        */
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
        // There was no overflow
        return ((int) sum);
}

/*
 * Safe fxp1 - fxp2 just using the safe add functions :)
 */
int fxp_sub(int fxp1, int fxp2)
{
        return fxp_add(fxp1, -fxp2);
}

int fxp_sub_l(int fxp1, int fxp2)
{
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


// Note: both fxp_mul_l and fxp_mul_d use longs,
// but fxp_mul_l is clearly more efficient.
// No reason to choose fxp_mul_d over it.

/*
 * Safe implementation of fxp multiplication checking for
 * overflows using divisions, and using longs.
 * Only applicable for systems in which
 * sizeof(long) > sizeof(int)
 *
int fxp_mul_d(int fxp1, int fxp2)
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
        if (fxp1 > 0) {
                if (fxp2 > 0) {
                        // both positive
                        if ((fxp2 > fxp_frac_mask) &&
                                (fxp1 > fxp_unsafe_div(FXP_MAX, fxp2))) {
                                return FXP_POS_INF;
                        }
                } else {
                        // fxp2 non-positive
                        if ((fxp1 > fxp_frac_mask) &&
                                (fxp2 < fxp_unsafe_div(FXP_MIN, fxp1))) {
                                return FXP_NEG_INF;
                        }
                }
        } else {
                // fxp1 non-positive
                if (fxp2 > 0) {
                        // fxp2 positive
                        if ((fxp2 > fxp_frac_mask) &&
                                (fxp1 < fxp_unsafe_div(FXP_MIN, fxp2))) {
                                return FXP_NEG_INF;
                        }
                } else {
                        // both non-positive
                        if ((-fxp1 > fxp_frac_mask) &&
                                (fxp2 < fxp_unsafe_div(FXP_MAX, fxp1))) {
                                return FXP_POS_INF;
                        }
                }
        }
        long prod = ((long) fxp1) * fxp2;
        return (int) (prod >> fxp_frac_bits);
}
*/

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
    /*
    int onebitnum, nb = 1;
    //printf("max_bits: %d\n", max_bits);
    int chunksize = FXP_INT_BITS_HALF;
    while (chunksize > 0) {
            //printf("\tchunksize: %d\n", chunksize);
            onebitnum = 1 << chunksize;
            if (x >= onebitnum) {
                    nb += chunksize;
                    x = (x >> chunksize);
                    //printf("\t\tnb:%d, x:%d\n", nb, x);
            }
            chunksize /= 2;
    }
    //printf("\tFinal nb:%d\n\n", nb);
    return nb;
    */
}

/*
//Simpler code to count the # of bits, but less efficient
int fxp_nbits_v0(unsigned int x, int max_bits)
{
    if (x == 0) return 0;
    int nb = FXP_INT_BITS - 1;
    while (!(x & (1 << nb))) nb--;
    return nb + 1;
}
*/

/*
 * Safe implementation of fxp multiplication using only ints.
 *
 * A distributive approach: computes n1 * n2 as
 * (w1 + f1)*(w2 + f2) == w1*w2 + w1*f2 + f1*w2 + f1*f2,
 * with n1 = w1 + f1, and n2 = w2 + f2
 * (their corresponding whole and frac parts added)
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
        // Compute product v1*v2 == (w1 + f1) * (w2 + f2) as
        // w1*w2 + w1*f2 + f1*w2 + f1*w2 (using only positive values)
        int v1, v2, w1, w2, bw1, bw2, product;
        v1 = (fxp1 >= 0)? fxp1: -fxp1;
        v2 = (fxp2 >= 0)? fxp2: -fxp2;
        // Whole parts for both operands
        w1 = fxp_get_whole_part(v1);
        w2 = fxp_get_whole_part(v2);
        // Bits used in whole parts
        bw1 = fxp_nbits(w1);
        bw2 = fxp_nbits(w2);
        //printf("\nfxp1=%d  \tfxp2=%d\n", fxp1, fxp2);
        //printf("fxp1=%x    \tfxp2=%x\n", fxp1, fxp2);
        //printf("w1=%d, w2=%d, m1=%d\n", w1, w2, m1);
        //printf("w1=%x, w2=%x, m1=%x\n", \
        //        w1, w2, m1);
        if (bw1 + bw2 > fxp_whole_bits) {
                // The product will for sure overflow just by multiplying
                // the whole parts.
                // Return appropriately signed infinity
                return ((fxp1 >= 0 && fxp2 >= 0) || (fxp1 < 0 && fxp2 < 0))?
                            FXP_POS_INF: FXP_NEG_INF;
        }
        unsigned int m1 = (unsigned int) (w1 * w2);
        if (m1 > fxp_whole_max) {
                // Overflow just by multiplying the whole parts
                // Return appropriately signed infinity
                return ((fxp1 >= 0 && fxp2 >= 0) || (fxp1 < 0 && fxp2 < 0))?
                            FXP_POS_INF: FXP_NEG_INF;
        }
        int m2, m3, m4, f1, f2, bf1, bf2, bf1v0, bf2v0;
        f1 = fxp_get_bin_frac(v1);
        f2 = fxp_get_bin_frac(v2);
        // Whole x frac parts cannot overflow, so simply multiply them
        m2 = w1 * f2;
        m3 = f1 * w2;
        // Bits used in frac parts
        bf1 = fxp_nbits(f1);
        bf2 = fxp_nbits(f2);
        //printf("f1=%d, f2=%d, bf1=%d, bf2=%d\n", f1, f2, bf1, bf2);
        // Multiplying the frac parts together could overflow only if
        // FXP_FRAC_BITS > FXP_WHOLE_BITS, but we can always drop
        // enough least-significant bits in the frac parts before
        // multiplying them (at the inevitable cost of some precision)
        // to avoid the overflow
        int lostbits = 0;
        while (bf1 + bf2 > FXP_INT_BITS_M1) {
                if (f1 > f2) {
                        f1 = (f1 >> 1);
                        bf1--;
                } else {
                        f2 = (f2 >> 1);
                        bf2--;
                }
                lostbits++;
        }
        m4 = (f1 * f2) >> (fxp_frac_bits - lostbits);
        //printf("f1=%d, f2=%d, bf1=%d, bf2=%d\n", f1, f2, bf1, bf2);
        //printf("f1=%x, f2=%x\n", f1, f2);
        //printf("m1=%d, m2=%d, m3=%d, m4=%d\n", m1, m2, m3, m4);

        // Even if none of the components overflowed, their sum
        // could overflow. So sum them all safely
        product = fxp_add(
                        fxp_add((m1 << fxp_frac_bits),
                                fxp_add(m2, m3)),
                        m4);
        //printf("product=%u\n", product);
        if (product == FXP_POS_INF) {
                // Overflow.
                // Return appropriately signed infinity
                return ((fxp1 >= 0 && fxp2 >= 0) || (fxp1 < 0 && fxp2 < 0))?
                            FXP_POS_INF: FXP_NEG_INF;
        }
        // No overflow, return the appropriately signed product
        return ((fxp1 >= 0 && fxp2 >= 0) || (fxp1 < 0 && fxp2 < 0))?
                    product: -product;
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
        int dividend = (fxp1 >= 0)? fxp1: -fxp1;
        int divisor = (fxp2 >  0)? fxp2: -fxp2;
        if (divisor == FXP_POS_INF)
                return (dividend == FXP_POS_INF)? FXP_UNDEF: 0;
        if (dividend == FXP_POS_INF)
                return ((fxp1 > 0 && fxp2 > 0) || (fxp1 < 0 && fxp2 < 0))?
                                FXP_POS_INF: FXP_NEG_INF;
        if ((divisor <= fxp_frac_mask) &&
                (dividend > fxp_mul(FXP_MAX, divisor))) {
                return ((fxp1 > 0 && fxp2 > 0) || (fxp1 < 0 && fxp2 < 0))?
                                FXP_POS_INF: FXP_NEG_INF;
        }
        // Calculation of division of fixed-point nums
        // using only ints

        int bits_dividend = fxp_nbits(dividend);
        int bits_divisor = fxp_nbits(divisor);

        //if (VERBOSE) printf("d: x = %d (hex: %x, %d bits)\n   y = %d (hex: %x, %d bits)\n", \
        //            dividend, dividend, bits_dividend, \
        //            divisor, divisor, bits_divisor);
        int size_diff = bits_dividend - bits_divisor;
        int shift = size_diff - 1;
        int minuend, next_bit_mask, bminuend, ba;
        int difference = 0;
        int bidx;
        //int qbits = 0;
        // Initialize starting minuend
        if (size_diff > 0) {
                int mmask = FXP_POS_INF >> (FXP_INT_BITS - bits_divisor - 1);
                minuend = (dividend & (mmask << size_diff)) >> size_diff;
                next_bit_mask = 1 << shift;
                bidx = bits_divisor;
                bminuend = bits_divisor;
        } else {
                // Minor speedup, making the first minuend have as
                // many bits as the divisor
                minuend = dividend << -size_diff;
                next_bit_mask = 0;
                bidx = bits_dividend - size_diff;
                bminuend = bits_dividend - size_diff;
                //qbits = -size_diff;
        }
        //if (VERBOSE) printf("\td: starting m: %d (%d bits)\n", minuend, bidx);
        int quotient = 0;
        int bmax = bits_dividend + fxp_frac_bits;
        //printf("bidx processed while <= bmax: %d\n", bmax);
        //int loops = 0;
        // bidx is (from left to right) the highest bit # we are currently
        // processing, so we loop till bidx exceeds the right-most bit
        // in the full (left-shifted) dividend
        while (bidx <= bmax) {
                //if (VERBOSE) printf("\tm: %d (hex:%x, %d bits) bidx:%d\n", minuend, minuend, bminuend, bidx);
                if (minuend >= divisor) {
                        // Append a 1 to the right of the quotient
                        quotient = (quotient << 1) | 1;
                        ba = 1;
                        //difference = minuend - divisor;
                        minuend = minuend - divisor;
                        bminuend = fxp_nbits(minuend);
                } else {
                        // Append a 0 to the right of the quotient
                        quotient = (quotient << 1);
                        ba = 0;
                        //difference = minuend;
                }
                //qbits = fxp_nbits(quotient);
                //if (VERBOSE) {
                //    trace_fxp_div("div:", loops, fxp_frac_bits, bidx, \
                //        dividend, divisor, \
                //        quotient, ba, qbits, \
                //        minuend, \
                //        ba, \
                //        difference );
                //}
                bidx++;
                //minuend = difference;
                //bminuend = fxp_nbits(minuend);
                if (bminuend == FXP_INT_BITS_M1) {
                        if (divisor > 1) {
                                divisor = (divisor >> 1);
                                //if (VERBOSE) printf("\t\tnew shrinked divisor %d\n\n", divisor);
                        } else {
                                // Overflow -> Return properly signed infinity
                                return ((fxp1 >= 0 && fxp2 > 0) || \
                                        (fxp1 < 0 && fxp2 < 0))? \
                                                FXP_POS_INF: FXP_NEG_INF;
                        }
                } else {
                        if (next_bit_mask > 0) {
                                // Pull down next bit from the dividend,
                                // and append it to the right of the minuend
                                minuend = (minuend << 1) | \
                                            ((dividend & next_bit_mask) \
                                                >> shift);
                                next_bit_mask = (next_bit_mask >> 1);
                                shift--;
                        } else {
                                // Pull down a 0 bit and append it to minuend
                                minuend = (minuend << 1);
                        }
                        bminuend = fxp_nbits(minuend);
                }
                //loops++;
        }
        // Return properly signed quotient
        return ((fxp1 >= 0 && fxp2 > 0) || (fxp1 < 0 && fxp2 < 0))?
                        quotient: -quotient;
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
        int v1 = (fxp1 >= 0)? fxp1: -fxp1;
        int v2 = (fxp2 >  0)? fxp2: -fxp2;
        if (v2 == FXP_POS_INF)
                return (v1 == FXP_POS_INF)? FXP_UNDEF: 0;
        if (v1 == FXP_POS_INF)
                return ((fxp1 > 0 && fxp2 > 0) || (fxp1 < 0 && fxp2 < 0))?
                                FXP_POS_INF: FXP_NEG_INF;
        // Compute positive division
        long numerator = ((long) v1) << fxp_frac_bits;
        long division = numerator / v2;

        if (division > FXP_MAX_L) {
                // Overflow -> Return properly signed infinity
                return ((fxp1 >= 0 && fxp2 > 0) || (fxp1 < 0 && fxp2 < 0))?
                                FXP_POS_INF: FXP_NEG_INF;
        }
        // No overflow -> return properly signed int
        return ((fxp1 >= 0 && fxp2 > 0) || (fxp1 < 0 && fxp2 < 0))?
                        (int) division: -((int) division);
}
