/*
 * fxp.c
 * An implementation of Binary Fixed Point numbers
 * encoding them into integers
 * By Raul Saavedra, 2022-11-13
 */

#include "fxp.h"
//#include <stdio.h>

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
        if (whole > FXP_WHOLE_MAX)
                return FXP_POS_INF;
        if (whole < FXP_WHOLE_MIN)
                return FXP_NEG_INF;
        int sign = 1;
        if ((whole==0) && (bin_frac < 0)) {
                // Special case for negative numbers when whole part is zero,
                // then the fxp gets its sign from the frac
                sign = -1;
                bin_frac = -bin_frac;
        }
        else {
                // All other cases, fxp gets its sign from the whole part
                if (whole < 0) {
                    sign = -1;
                    whole = -whole;
                }
                if (bin_frac < 0)
                    bin_frac = -bin_frac;
        }
        int positive_fxp = (whole << FXP_FRAC_BITS) | bin_frac;
        //printf("fxp_from_bin_frac: frac is: %d\n", frac);
        return (sign==1)? positive_fxp: -positive_fxp;
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
        //printf("   fxp_from_dec_frac: original frac: %d\n", dec_frac);
        while (trunc_frac > FXP_FRAC_MAX_DEC) {
                //trimmed_frac = (trimmed_frac + 5)/10;
                trunc_frac = trunc_frac/10;
                //printf("   fxp_from_dec_frac: frac trimmed to: %d\n", trunc_frac);
        }
        int bin_frac = (trunc_frac * FXP_FRAC_MAX) / FXP_FRAC_MAX_DEC ;
        //printf("   fxp_from_dec_frac: bin_frac is: %d\n", bin_frac);
        return fxp_bin(whole, (frac_sign==1)? bin_frac: -bin_frac);
}

int fxp_get_whole_part(int fxp)
{
        if (fxp < 0)
                return -((-fxp) >> FXP_FRAC_BITS);
        else
                return fxp >> FXP_FRAC_BITS;
}

/*
 * Get the frac part directly (binary)
 */
int fxp_get_bin_frac(int fxp)
{
        if (fxp < 0)
                return -((-fxp) & FXP_FRAC_MAX);
        else
                return (fxp & FXP_FRAC_MAX);
}

/*
 * Get the frac part as decimal between 0 and FXP_FRAC_MAX_DEC, e.g. 0 .. 999
 */
int fxp_get_dec_frac(int fxp)
{
        if (fxp < 0)
                return -(((-fxp) & FXP_FRAC_MAX) * FXP_FRAC_MAX_DEC)/FXP_FRAC_MAX;
        else
                return ((fxp & FXP_FRAC_MAX) * FXP_FRAC_MAX_DEC)/FXP_FRAC_MAX;
}

/*
 * The simpler unsafe implementations of sum, sub, mul, and div.
 * No overflow checks, and minimal or no special considerations
 * for infinities or FXP_UNDEF. Efficient but dangerous
 */

int fxp_unsafe_sum(int fxp1, int fxp2)
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
        return (int) (p >> FXP_FRAC_BITS);
}

int fxp_unsafe_div(int fxp1, int fxp2)
{
        if (fxp2 == 0)
                return (fxp1 > 0)? FXP_POS_INF:
                                   (fxp1 < 0)? FXP_NEG_INF: FXP_UNDEF;
        long n1 = ((long) fxp1) << FXP_FRAC_BITS;
        long div = n1 / fxp2;
        return ((int) div);
}

/*
 * Safe implementation of fxp1 + fxp2
 */
int fxp_sum(int fxp1, int fxp2)
{
        // Check for any undef or infinity arguments
        if (fxp1 == FXP_UNDEF || fxp2 == FXP_UNDEF)
                return FXP_UNDEF;
        if (fxp1 == FXP_POS_INF || fxp2 == FXP_POS_INF)
                return (fxp1 == FXP_NEG_INF || fxp2 == FXP_NEG_INF)?
                            FXP_UNDEF: FXP_POS_INF;
        if (fxp1 == FXP_NEG_INF || fxp2 == FXP_NEG_INF)
                return FXP_NEG_INF;
        // Check if arguments have the same sign
        if ((fxp1 >= 0 && fxp2 >= 0) || (fxp1 < 0) && (fxp2 < 0)) {
                // Same signs, check for possible overflow
                if (fxp1 > 0)
                        return (FXP_MAX - fxp1 <=  fxp2)?
                                    FXP_POS_INF: fxp1 + fxp2;
                else
                        return (FXP_MAX + fxp1 <= -fxp2)?
                                    FXP_NEG_INF: fxp1 + fxp2;
        };
        // Arguments with different signs, no overflow danger, simply sum
        return fxp1 + fxp2;
}

/*
 * Safe implementation of fxp1 - fxp2 (just using the safe sum :)
 */
int fxp_sub(int fxp1, int fxp2)
{
        return fxp_sum(fxp1, -fxp2);
}

/*
Return number of bits used by the value in the argument
(returned value is between 0 and FXP_INT_BITS)
*/
unsigned int nbits(unsigned int x, unsigned int max_bits)
{
    if (x == 0) return 0;
    unsigned int nb = 1;
    unsigned int chunksize = max_bits / 2;
    while (chunksize > 0) {
        unsigned int onebitnum = 1 << chunksize;
        if (x >= onebitnum) {
            nb += chunksize;
            x = (x >> chunksize);
        }
        chunksize /= 2;
    }
    return nb;
}


/*
Simpler code to count the # of bits, but slightly less efficient
*/
unsigned int nbits_v0(unsigned int x)
{
    if (x == 0) return 0;
    int nb = FXP_INT_BITS - 1;
    while (!(x & (1 << nb))) nb--;
    return nb + 1;
}

unsigned int fxp_sub_mul(unsigned int v1,
                        unsigned int bitsv1,
                        unsigned int v2,
                        unsigned int bitsv2)
{
    if (bitsv1 + bitsv2 < FXP_INT_BITS) {
        return (v1 * v2);
        }
    else {
        return FXP_POS_INF;
    }
}

/*
 * Safe implementation of fxp1 * fxp2
 */
int fxp_mul(int fxp1, int fxp2)
{
        if (fxp1 == FXP_UNDEF || fxp2 == FXP_UNDEF)
                return FXP_UNDEF;
        if (fxp1 == FXP_POS_INF || fxp2 == FXP_POS_INF) {
                if (fxp1 == 0 || fxp2 == 0)
                        return FXP_UNDEF;
                return (fxp1 == FXP_NEG_INF || fxp2 == FXP_NEG_INF)?
                            FXP_NEG_INF: FXP_POS_INF;
        }
        if (fxp1 == FXP_NEG_INF || fxp2 == FXP_NEG_INF) {
                if (fxp1 == 0 || fxp2 == 0)
                        return FXP_UNDEF;
                return FXP_NEG_INF;
        }

        /* Implementation of multiplication supporting systems
        for which sizeof(long) == sizeof(int) */
        // Compute product v1*v2 == (w1 + f1) * (w2 + f2) as
        // w1*w2 + w1*f2 + f1*w2 + f1*w2 (using only positive values)
        unsigned int v1 = (fxp1 >= 0)? fxp1: -fxp1;
        unsigned int v2 = (fxp2 >= 0)? fxp2: -fxp2;
        // Whole and fraction parts
        unsigned int w1 = fxp_get_whole_part(v1);
        unsigned int f1 = fxp_get_bin_frac(v1);
        unsigned int w2 = fxp_get_whole_part(v2);
        unsigned int f2 = fxp_get_bin_frac(v2);
        // Bits used in each of them
        //unsigned int bw1 = nbits_v0(w1);
        //unsigned int bw2 = nbits_v0(w2);
        //unsigned int bf1 = nbits_v0(f1);
        //unsigned int bf2 = nbits_v0(f2);
        unsigned int bw1 = nbits(w1, FXP_WHOLE_BITS);
        unsigned int bw2 = nbits(w2, FXP_WHOLE_BITS);
        unsigned int bf1 = nbits(f1, FXP_FRAC_BITS);
        unsigned int bf2 = nbits(f2, FXP_FRAC_BITS);
        unsigned int m1 = fxp_sub_mul(w1, bw1, w2, bw2);
        unsigned int m2 = fxp_sub_mul(w1, bw1, f2, bf2);
        unsigned int m3 = fxp_sub_mul(f1, bf1, w2, bw2);
        unsigned int m4 = fxp_sub_mul(f1, bf1, f2, bf2);
        int product = fxp_sum(fxp(m1), fxp_sum(m2, fxp_sum(m3, m4)));
        if (product == FXP_POS_INF) {
            // Overflow. Return infinity with the appropriate sign
            return ((fxp1 >= 0 && fxp2 >= 0) || (fxp1 < 0 && fxp2 < 0))?
                        FXP_POS_INF: FXP_NEG_INF;
        }
        else {
            return ((fxp1 >= 0 && fxp2 >= 0) || (fxp1 < 0 && fxp2 < 0))?
                        product: -product;
        }

        /*
        // This code only applicable for systems where sizeof(long) > sizeof(int)
        long product = ((long) v1) * v2;
        //if (product <= (((long) FXP_MAX) << FXP_FRAC_BITS)) {
        if (product <= FXP_MAX_LSHIFTED) {
                // No overflow, return result as int with appropriate sign
                product = product >> FXP_FRAC_BITS;
                return ((fxp1 >= 0 && fxp2 >= 0) || (fxp1 < 0 && fxp2 < 0))?
                            (int) product: -((int) product);
        } else {
                // We overflowed. Return infinity with the appropriate sign
                return ((fxp1 >= 0 && fxp2 >= 0) || (fxp1 < 0 && fxp2 < 0))?
                            FXP_POS_INF: FXP_NEG_INF;
        }
        */
}

/*
 * Safe implementation of fxp1 / fxp2
 */
int fxp_div(int fxp1, int fxp2)
{
        if (fxp1 == FXP_UNDEF || fxp2 == FXP_UNDEF)
                return FXP_UNDEF;
        if (fxp2 == 0)
                return (fxp1 > 0)? FXP_POS_INF:
                                   (fxp1 < 0)? FXP_NEG_INF: FXP_UNDEF;
        // Compute positive division
        int v1 = (fxp1 >= 0)? fxp1: -fxp1;
        int v2 = (fxp2 >= 0)? fxp2: -fxp2;
        if (v2 == FXP_POS_INF)
                return (v1 == FXP_POS_INF)? FXP_UNDEF: 0;
        if (v1 == FXP_POS_INF)
                return ((fxp1 >= 0 && fxp2 >=0) || (fxp1 < 0 && fxp2 < 0))?
                            FXP_POS_INF: FXP_NEG_INF;
        // This code only applicable for systems where sizeof(long) > sizeof(int)
        long numerator = ((long) v1) << FXP_FRAC_BITS;
        long div = numerator / v2;
        if (div <= FXP_MAX_L)
                // No overflow -> return result as int appropriately signed
                return ((fxp1 >= 0 && fxp2 >= 0) || (fxp1 < 0 && fxp2 < 0))?
                            (int) div: -((int) div);
        else
                // Overflow -> Return signed infinity
                return ((fxp1 >= 0 && fxp2 >= 0) || (fxp1 < 0 && fxp2 < 0))?
                            FXP_POS_INF: FXP_NEG_INF;
}
