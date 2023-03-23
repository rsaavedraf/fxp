/* SPDX-License-Identifier: MIT */
/*
 * fxp_l.c
 * Some of the functions for Binary Fixed Point numbers
 * using longs for speed (and also simplicity)
 *
 * These functions can/should only be used when the
 * target system's sizeof(long) >= 2 * sizeof(int)
 *
 * By Raul Saavedra, Bonn, Germany
 */

#include "fxp_l.h"
#include "fxp_extern.h"
#include "fxp_conv.h"
#include "fxp_aux.h"
#include <stdio.h>
#include <stdlib.h>
//#include <assert.h>

//#include "fxp_aux.h"
//#define VERBOSE 0

struct lgcharm_l {
        long characteristic;
        unsigned long mantissa;
};

// For the BKM lg2 calculation when using longs.
// Values in this array are: a[k] = lg2(1 + 1/2^k) represented as
// unsigned long fxp's (8 bytes,) with 63 frac bits, and
// one magnitude whole bit (notice, not a sign but a magnitude single
// whole bit, needed to be able to represent the first value in
// the table == 1)
static const unsigned long FXP_BKM_LOGS_L[] = {
        0x8000000000000000, 0x4AE00D1CFDEB4400, 0x2934F0979A371600, 0x15C01A39FBD68800,
        0xB31FB7D64898B00, 0x5AEB4DD63BF61C0, 0x2DCF2D0B85A4540, 0x16FE50B6EF08510,
        0xB84E236BD563B8, 0x5C3E0FFC29D594, 0x2E24CA6E87E8A8, 0x1713D62F7957C3,
        0xB8A476150DFE4, 0x5C53AC47E94D8, 0x2E2A32762FA6B, 0x1715305002E4A,
        0xB8A9DED47C11, 0x5C55067F6E58, 0x2E2A89050622, 0x171545F3D72B,
        0xB8AA35640A7, 0x5C551C23599, 0x2E2A8E6E01E, 0x1715474E163,
        0xB8AA3ACD06, 0x5C551D7D98, 0x2E2A8EC491, 0x17154763BA,
        0xB8AA3B239, 0x5C551D933, 0x2E2A8EC9F, 0x171547651,
        0xB8AA3B28, 0x5C551D94, 0x2E2A8ECA, 0x17154765,
        //0xB8AA3B2, 0x5C551D9, 0x2E2A8EC, 0x1715476,       // <---- *
        //0xB8AA3B, 0x5C551D, 0x2E2A8E, 0x171547,
        //0xB8AA3, 0x5C551, 0x2E2A8, 0x17154,
        //0xB8AA, 0x5C55, 0x2E2A, 0x1715,
        //0xB8A, 0x5C5, 0x2E2, 0x171,
        //0xB8, 0x5C, 0x2E, 0x17,
        //0xB, 0x5, 0x2, 0x1,
        //0x0
// Starting with the row marked with the *, each entry is exactly
// a 4-bit right-shift of the value 4 positions earlier
};


static inline int fxp_nbits_l(unsigned long x)
{
        if (x == 0) return 0;
        return FXP_LONG_BITS - __builtin_clzl(x);
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
        // No overflow, return rounded result as int with appropriate sign
        int rbit = (product >> FXP_frac_bits_m1) & 1;
        product = (product >> FXP_frac_bits) + rbit;
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
        int nx, c;
        int nbx = fxp_nbits(fxp1);
        if (fxp1 < FXP_one) {
                c = nbx - FXP_frac_bits - 1;
                nx = fxp1 << (-c);
        } else if (fxp1 >= FXP_two) {
                c = nbx - FXP_frac_bits - 1;
                nx = fxp1 >> c;
        } else {
                c = 0;
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
 * Internal auxiliary function that calculates the characteristic
 * and rounded mantissa of lg2, returning them separately in a
 * struct (lgcharm_l), both at their maximum precision:
 * The characteristic as long with 0 frac bits
 * The mantissa as unsigned long with 63 frac bits (exact same
 * configuration of the BKM array values)
 * Uses the BKM L-Mode algorithm (requires table of pre-calculated
 * values) and longs.
 *
 * For more details on BKM:
 * https://en.wikipedia.org/wiki/BKM_algorithm
 */
static inline struct lgcharm_l fxp_lg2_l_tuple(int fxp1)
{
        struct lgcharm_l result;
        // Assumes fxp1 is for sure > 0
        int clz = __builtin_clz((unsigned int) fxp1);
        // clz > 0 since at least sign bit == 0
        int nbx = FXP_INT_BITS - clz;
        // Assign the characteristic to *c
        result.characteristic = \
                ((fxp1 < FXP_one) || (fxp1 >= FXP_two))? \
                        (nbx - FXP_frac_bits - 1): 0;
        // Here we have already calculated the log characteristic c
        // Now lshifting the original number to have it as unsigned
        // long fxp with 2 whole bits (2 needed because z can
        // intermitently get a value > 2)
        // Watch out we are using two different shifted
        // configurations simultaneously, but independently:
        // - x aligned with FXP_BKM_ONE_L: 2 whole and 62 frac bits
        //   (again at least 2 whole bits for the BKM log algorithm).
        // - However argument, the BKM_LOGS array values, and the
        //   mantissa all have 1 unsigned whole and 63 frac bits
        //   for more accuracy
        unsigned long argument = ((unsigned long) fxp1) << \
                            (clz + FXP_INT_BITS_M1);
        unsigned long x = FXP_BKM_ONE_L;
        unsigned long z, xs;
        // The mantissa value will remain in the range of lg2(x),
        // x in [1, 2), meaning mantissa always in [0, 1)
        result.mantissa = 0lu;
        for (int shift = 1; shift <= FXP_INT_BITS; shift++) {
                //printf("shift:%d looping\n", shift);
                xs = (x >> shift);
                z = x + xs;
                if (z <= argument) {
                        x = z;
                        result.mantissa += FXP_BKM_LOGS_L[shift];
                        //printf("\tx:x%lX  Updating y:x%lX\n", x, y);
                }
        }
        //printf("characteristic:%ld (x%lX)  mantissa:%lu (x%lX)\n", \
        //            result.characteristic, result.characteristic, \
        //            result.mantissa, result.mantissa);
        //printf("Lf mantissa:%LE\n", ((long double) result.mantissa) / (~0ul >> 1));
        return result;
}

/*
 * lg2 using BKM and longs
 */
int fxp_lg2_l(int fxp1)
{
        if (fxp1 < 0) return FXP_UNDEF;
        if (fxp1 == 0) return FXP_NEG_INF;
        if (fxp1 == FXP_POS_INF) return FXP_POS_INF;
        // Get the separate characteristic and full mantissa
        struct lgcharm_l charm = fxp_lg2_l_tuple(fxp1);
        // Round and shift the mantissa
        unsigned long rbit = (charm.mantissa >> \
                                FXP_lg2_l_mshift_m1) & 1ul;
        long longm = (charm.mantissa >> FXP_lg2_l_mshift) + rbit;
        int m = (int) longm;
        //printf("Final mantissa:%d  (rounding bit was %d)\n", m, (int) rbit);
        // Return the complete logarithm (char + mant) as fxp
        if (charm.characteristic < FXP_whole_min) {
                if ((charm.characteristic == FXP_whole_min_m1) \
                        && (m > 0)) {
                        // Special case: in spite of not enough whole
                        // bits available for the characteristic, there
                        // are for characteristic + 1, which is what will
                        // get returned as whole part of this logarithm
                        return INT_MIN + m;
                }
                // Overflow: the characteristic of the logarithm
                // will not fit in the whole bits part of our
                // current fxp settings
                return FXP_NEG_INF;
        }
        return (charm.characteristic >= 0)?
                (int) ((charm.characteristic << FXP_frac_bits) + m):
                (int) (-(-charm.characteristic << FXP_frac_bits) + m);
}

/*
 * Analogougs to mul_distrib in fxp.c, but here using longs
 */
unsigned long mul_distrib_l( unsigned long x,
                             unsigned long y)
{
        unsigned long xa, xb, ya, yb, xaya, xayb, yaxb, xbyb;
        unsigned long qr1, qr2, qr3, qrsum, ql1, ql2, ql3, product;
        xa = x >> FXP_INT_BITS;
        xb = (x & FXP_RINT_MASK);
        ya = y >> FXP_INT_BITS;
        yb = (y & FXP_RINT_MASK);
        //if (VERBOSE) printf("\txa:%lX, xb:%lX, ya:%lX, yb:%lX\n", xa, xb, ya, yb);
        xayb = xa * yb;
        yaxb = ya * xb;
        qr1 = xayb & FXP_RINT_MASK;
        qr2 = yaxb & FXP_RINT_MASK;
        qr3 = xbyb >> FXP_INT_BITS;
        int rbit = ((unsigned int) (xbyb >> FXP_INT_BITS_M1)) & 1;
        qrsum = qr1 + qr2 + qr3 + rbit;
        xaya = xa * ya;
        xbyb = xb * yb;
        ql1 = xayb >> FXP_INT_BITS;
        ql2 = yaxb >> FXP_INT_BITS;
        ql3 = (qrsum >> FXP_INT_BITS);
        product = xaya + ql1 + ql2 + ql3;
        return product;
}

/*
 * Calculates log2 and then multiplies by the given factor
 */
static inline int lg2_x_factor_l(int fxp1, const unsigned long FACTOR)
{
        struct lgcharm_l charm = fxp_lg2_l_tuple(fxp1);
        int shift_for_c;
        unsigned long shifted_c;
        long cxf;
        if (charm.characteristic < 0) {
                if ((charm.characteristic < FXP_whole_min_m1) \
                        || ((charm.characteristic == FXP_whole_min_m1) \
                            && (charm.mantissa == 0))) {
                        return FXP_NEG_INF;
                }
                unsigned long posc = (-charm.characteristic);
                shift_for_c = __builtin_clzl(posc) - 1;
                shifted_c = posc << shift_for_c;
                cxf = -mul_distrib_l(shifted_c, FACTOR);

        } else {
                shift_for_c = __builtin_clzl(charm.characteristic) - 1;
                shifted_c = charm.characteristic << shift_for_c;
                cxf = mul_distrib_l(shifted_c, FACTOR);

        }
        unsigned long mxf = mul_distrib_l(charm.mantissa, FACTOR);
        int shift_for_m = FXP_LONG_BITS_M1 - shift_for_c;
        int rbit = (shift_for_m == 0)? 0: (mxf >> (shift_for_m - 1)) & 1ul;
        long shifted_mxf = (long) ((mxf >> shift_for_m) + rbit);
        long sum = cxf + shifted_mxf;
        long final_lg;
        // Finally round and shift for current fxp configuration
        if (sum < 0) {
                long psum = -sum;
                rbit = (psum >> (shift_for_c - FXP_frac_bits_m1)) & 1ul;
                final_lg =  -((psum >> (shift_for_c - FXP_frac_bits)) + rbit);
        } else {
                rbit = (sum >> (shift_for_c - FXP_frac_bits_m1)) & 1ul;
                final_lg =  (sum >> (shift_for_c - FXP_frac_bits)) + rbit;
        }
        return (int) final_lg;
}

/*
 * Implementation of ln_l using lg2
 */
int fxp_ln_l(int fxp1)
{
        if (fxp1 < 0) return FXP_UNDEF;
        if (fxp1 == 0) return FXP_NEG_INF;
        if (fxp1 == FXP_POS_INF) return FXP_POS_INF;
        return lg2_x_factor_l(fxp1, FXP_LN_2_I64);
}

/*
 * Implementation of lg10 using lg2
 */
int fxp_lg10_l(int fxp1) {
        if (fxp1 < 0) return FXP_UNDEF;
        if (fxp1 == 0) return FXP_NEG_INF;
        if (fxp1 == FXP_POS_INF) return FXP_POS_INF;
        return lg2_x_factor_l(fxp1, FXP_LG10_2_I64);
}

/*
 * pow2 calculation (2^fxp1) using the BKM (E-mode) algorithm
 * Analogous to implementation in bkm.c, but here tailored
 * for fxp's, and using longs.
 */
int fxp_pow2_l(int fxp1)
{
        if (fxp1 == FXP_UNDEF) return FXP_UNDEF;
        if (fxp1 == FXP_NEG_INF) return 0;
        if (fxp1 == FXP_POS_INF) return FXP_POS_INF;
        int w = fxp_get_whole_part(fxp1);
        int f = fxp_get_bin_frac(fxp1);
        unsigned long pow2w, argument;
        if (fxp1 >= 0) {
                if (w >= FXP_whole_bits_m1) {
                        return FXP_POS_INF;
                }
                pow2w = FXP_one_l << (w + 2);
                // Argument will be in [0, 1)
                argument = ((unsigned long) f) << FXP_lg2_l_mshift;
        } else {
                if (w <= FXP_INT_BITS_M1_NEG) {
                        return 0;
                }
                if (w < 0)
                        pow2w = FXP_one_l >> (-w - 1);
                else
                        pow2w = FXP_one_l << 1;
                // Notice argument is > 0, in (0, 1]
                argument = (FXP_one_l + f) \
                                << FXP_lg2_l_mshift;
        }
        //if (VERBOSE) printf("\npow2_l: pow2w:%lX  argument:%lX\n", pow2w, argument);
        unsigned long x = FXP_BKM_ONE_L, y = 0;
        for (int k = 0; k < FXP_INT_BITS; k++) {
                unsigned long const  z = y + FXP_BKM_LOGS_L[k];
                //if (VERBOSE) printf("k:%d,  z:%lX\n", k, z);
                if (z <= argument) {
                        y = z;
                        x = x + (x >> k);
                        //if (VERBOSE) printf("\tUpdating y (%lX) und x (%lX)\n", y, x);
                }
        }
        //if (VERBOSE) printf("final x:%lX\n", x);
        unsigned long md = mul_distrib_l(pow2w, x);
        return (int) md;
}
