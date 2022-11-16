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
 */
int fxp_bin(int whole, int frac)
{
        if (whole > FXP_WHOLE_MAX)
                return FXP_POS_INF;
        if (whole < FXP_WHOLE_MIN)
                return FXP_NEG_INF;
        int sign = 1;
        if ((whole==0) && (frac < 0)) {
                // Special case for negative numbers when whole part is zero,
                // then the fxp gets its sign from the frac
                sign = -1;
                frac = -frac;
        }
        else {
                // All other cases, fxp gets its sign from the whole part
                if (whole < 0) {
                    sign = -1;
                    whole = -whole;
                }
                if (frac < 0)
                    frac = -frac;
        }
        int positive_fxp = (whole << FXP_FRAC_BITS) | frac;
        //printf("fxp_from_bin_frac: frac is: %d\n", frac);
        return (sign==1)? positive_fxp: -positive_fxp;
}

/*
 * Create an fxp number given its whole and (decimal) frac parts.
 * dec_frac is expected to be a decimal number between 0 and
 * FXP_FRAC_DECS-1, e.g. between 0 and 999
 * Usage examples:
 *     For fxp=16.001, you would invoke: fxp_from_dec_frac(16, 1)
 *     For 20.09: fxp_from_dec_frac(24, 90)
 *     For 24.5:  fxp_from_dec_frac(24, 500)
 * Note last example frac is not 5 but 500. (5 would correspond to 24.005)
 * This decimal value gets scaled into the binary range available for frac.
 * If frac is too large it will get trimmed-rounded to its most
 * significant digits until under the value of FXP_FRAC_DECS, e.g. for a
 * scale of 1000, if frac=987654 it will be trimmed to ~988
 */
int fxp_dec(int whole, int dec_frac)
{
        int frac_sign = 1;
        if (dec_frac < 0) {
                frac_sign = -1;
                dec_frac = -dec_frac;
        }
        int trimed_frac = dec_frac;
        while (trimed_frac >= FXP_FRAC_DECS)
                trimed_frac = (trimed_frac + 5)/10;
        int bin_frac = (trimed_frac * FXP_FRAC_MAX) / FXP_FRAC_DECS;
        //printf("fxp_from_dec_frac: bin_frac is: %d\n", bin_frac);
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
 * Get the frac part as decimal between 0 and FXP_FRAC_DECS, e.g. 0 .. 999
 */
int fxp_get_dec_frac(int fxp)
{
        if (fxp < 0)
                return -(((-fxp) & FXP_FRAC_MAX) * FXP_FRAC_DECS)/FXP_FRAC_MAX;
        else
                return ((fxp & FXP_FRAC_MAX) * FXP_FRAC_DECS)/FXP_FRAC_MAX;
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
        if (fxp2 == FXP_NEG_INF || fxp2 == FXP_NEG_INF)
                return FXP_NEG_INF;
        // Check if arguments have the same sign
        if ((fxp1 >= 0 && fxp2 >= 0) || (fxp1 < 0) && (fxp2 < 0)) {
                // Same signs, check for possible overflow
                if (fxp1 > 0)
                        return (FXP_MAX-fxp1 <=  fxp2)?
                                    FXP_POS_INF: fxp1 + fxp2;
                else
                        return (FXP_MAX+fxp1 <= -fxp2)?
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
        // Compute positive product
        int v1 = (fxp1 >= 0? fxp1: -fxp1);
        int v2 = (fxp2 >= 0? fxp2: -fxp2);
        long product = ((long) v1) * v2;
        if (product <= (((long) FXP_MAX) << FXP_FRAC_BITS)) {
                // No overflow, return result as int with appropriate sign
                product = product >> FXP_FRAC_BITS;
                return ((fxp1 >= 0 && fxp2 >= 0) || (fxp1 < 0 && fxp2 < 0))?
                            (int) product: -((int) product);
        } else {
                // We overflowed. Return infinity with the appropriate sign
                return ((fxp1 >= 0 && fxp2 >= 0) || (fxp1 < 0 && fxp2 < 0))?
                            FXP_POS_INF: FXP_NEG_INF;
        }
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
        int v1 = (fxp1 >= 0? fxp1: -fxp1);
        int v2 = (fxp2 >= 0? fxp2: -fxp2);
        if (v2 == FXP_POS_INF)
                return (v1 == FXP_POS_INF)? FXP_UNDEF: 0;
        if (v1 == FXP_POS_INF)
                return ((fxp1 >= 0 && fxp2 >=0) || (fxp1 < 0 && fxp2 < 0))?
                            FXP_POS_INF: FXP_NEG_INF;
        long numerator = ((long) v1) << FXP_FRAC_BITS;
        long div = numerator / v2;
        if (div <= (long) FXP_MAX)
                // No overflow -> return result as int appropriately signed
                return ((fxp1 >= 0 && fxp2 >= 0) || (fxp1 < 0 && fxp2 < 0))?
                            (int) div: -((int) div);
        else
                // Overflow -> Return signed infinity
                return ((fxp1 >= 0 && fxp2 >= 0) || (fxp1 < 0 && fxp2 < 0))?
                            FXP_POS_INF: FXP_NEG_INF;
}
