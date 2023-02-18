/*
 * fxp_l.c
 * Some of the functions for Binary Fixed Point numbers
 * using longs for speed.
 *
 * These can be used only when the target system's
 * sizeof(long) >= 2 * syzeof(int)
 *
 * By Raul Saavedra, Bonn, Germany
 */

#include "fxp_l.h"
#include "fxp_tconst.h"
//#include <stdio.h>
//#include <stdlib.h>
//#include <assert.h>

//Used while testing and debugging trying to optimize division
//#include "fxp_aux.h"
//#define VERBOSE 1

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

/**
 * Safe implementation of fxp multiplication using longs,
 * and no divisions.
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
        if (product > FXP_max_lshifted) {
                // Overflow, return infinity with the appropriate sign
                return ((fxp1 >= 0 && fxp2 >= 0) || \
                            (fxp1 < 0 && fxp2 < 0))?
                                FXP_POS_INF: FXP_NEG_INF;
        }
        // No overflow, return result as int with appropriate sign
        product = product >> FXP_frac_bits;
        return ((fxp1 >= 0 && fxp2 >= 0) || (fxp1 < 0 && fxp2 < 0))?
                        (int) product: -((int) product);
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
        long numerator = ((long) x) << FXP_frac_bits;
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

/*
 * Implementation of ln and lg10 using longs
 */
int fxp_ln_l(int fxp1)
{
        return fxp_mul_l(fxp_lg2_l(fxp1), FXP_shifted_ln_2);
}

int fxp_lg10_l(int fxp1) {
        return fxp_mul_l(fxp_lg2(fxp1), FXP_shifted_lg10_2);
}

/*
 * Alternative lg2 of an fxp using multiplication. Only applicable
 * for systems in which sizeof(long) >= 2 * sizeof(int)
 *
 * Requires the fxp configuration to have 3 or more whole bits.
 * Needs no pre-calculated table of log values in any range, but
 * requires one multiplication per mantissa bit (expensive.)
 *
 * Adapted to fxps from the general algorithm to calculate
 * binary logarithms explained by Clay Turner in IEEE Signal
 * Processing Magazine, Sep/2010.
 * D. E. Knuth's "The Art of Computer Programming Vol 2:
 * Seminumerical Algorithms", 2nd ed., 1981 (pages 441 - 446) is
 * the one and only reference in that short article, so that's the
 * effective reference.
 */
int fxp_lg2_mul_l(int fxp1)
{
        if (fxp1 <= 0) return ((fxp1 == 0)? FXP_NEG_INF: FXP_UNDEF);
        if (fxp1 == FXP_POS_INF) return FXP_POS_INF;
        int c = 0; // characteristic
        int nx;
        int nbx = fxp_nbits(fxp1);
        if (fxp1 < FXP_one) {
                c = nbx - FXP_frac_bits - 1;
                nx = fxp1 << (-c);
        } else if (fxp1 >= FXP_two) {
                c = nbx - FXP_frac_bits - 1;
                nx = fxp1 >> c;
        } else {
                nx = fxp1;
        }
        // Mantissa calculation:
        int m = 0;
        int b = FXP_half;
        while (b > 0) {
                nx = fxp_mul_l(nx, nx);
                if (nx >= FXP_two) {
                        nx = nx >> 1;
                        m |= b;
                }
                b = b >> 1;
        }
        // Return the calculated logarithm as fxp
        if (c < FXP_whole_min) {
                if ((c + 1 == FXP_whole_min) && (m > 0)) {
                    return INT_MIN + m;
                }
                return FXP_NEG_INF;
        }
        if (c >= 0)
            return (c << FXP_frac_bits) + m;
        else
            return -(-c << FXP_frac_bits) + m;
}

/*
 * log2 calculation using the BKM algorithm, and longs.
 * Only applicable for systems in which sizeof(long) >= 2 * sizeof(int)
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
        int c = ((fxp1 < FXP_one) || (fxp1 >= FXP_two))?
                    (nbx - FXP_frac_bits - 1): 0;

        // Here we have already calculated the log characteristic c
        // Now lshifting the original number to have it as unsigned
        // long with 61 frac bits (same as the pre-calculated values
        // in our bkm table of logs)
        unsigned long z = ((unsigned long) fxp1) << \
                            (clz + FXP_BKM_L_CLZSHIFT);
        unsigned long x = FXP_BKM_ONE_L;
        unsigned long zz;
        unsigned long y = 0; // starting mantissa
        for (int shift = 1; shift <= FXP_lg2_l_maxloops; shift++) {
                unsigned long xs = (x >> shift);
                zz = x + xs;
                if (zz <= z) {
                        // Here we use the precalculated values
                        unsigned long aux = FXP_BKM_LOGS_L[shift];
                        x = zz;
                        y += aux;
                }
        }
        // Now y has the mantissa, shift it adjusting for final fxp
        long m;
        if (FXP_lg2_l_shift >= 0) {
                m = y >> FXP_lg2_l_shift;
        } else {
                m = y << (-FXP_lg2_l_shift);
        }
        int im = (int) m;
        // Return the complete logarithm (c + m) as fxp
        if (c < FXP_whole_min) {
                if ((c + 1 == FXP_whole_min) && (im > 0)) {
                        // Special case: in spite of not enough whole
                        // bits available for c, there are for c + 1,
                        // which is what this logarithm will return as
                        // whole part
                        return INT_MIN + im;
                }
                // Overflow: the characteristic of the logarithm
                // will not fit in the whole bits part of our
                // current fxp settings
                return FXP_NEG_INF;
        }
        if (c >= 0)
            return (c << FXP_frac_bits) + im;
        else
            return -(-c << FXP_frac_bits) + im;
}
