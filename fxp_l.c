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
//#include "fxp_tconst.h"
#include "fxp_conv.h"
#include "fxp_aux.h"
#include <stdio.h>
#include <stdlib.h>
//#include <assert.h>

//Used while testing and debugging trying to optimize division
//#include "fxp_aux.h"
//#define VERBOSE 1


// For the BKM log calculation when using longs.
// Values in this array are: a[k] = log2(1 + 1/2^k) represented as
// unsigned long fxp's (8 bytes,) with 63 frac bits
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
 * log2 calculation using the BKM (L-Mode) algoritm,
 * and using longs.
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
        // long fxp with 2 whole bits (2 needed because zz can
        // intermitently get a value > 2)
        unsigned long z = ((unsigned long) fxp1) << \
                            //(clz + FXP_BKM_L_CLZSHIFT);
                            (clz + FXP_INT_BITS_M1);
        unsigned long x = FXP_BKM_ONE_L;
        unsigned long zz, xs;
        // y is the starting mantissa. Its value will remain in the
        // range of lg(x), x in [1, 2), so y will always be in [0, 1)
        unsigned long y = 0;
        //printf("z: %lu\n", z);
        //for (int shift = 1; shift <= FXP_lg2_maxloops; shift++) {
        for (int shift = 1; shift <= FXP_INT_BITS; shift++) {
                //printf("shift:%d looping\n", shift);
                xs = (x >> shift);
                zz = x + xs;
                if (zz <= z) {
                        x = zz;
                        y += FXP_BKM_LOGS_L[shift];
                        //printf("\tx:%lX  Updating y:%lX\n", x, y);
                }
        }
        // Here y has the calculated mantissa

        // TODO:
        // For ln() or lg10() we could apply the multiplication
        // with the corresponding factor at this point,
        // before rshifting the mantissa to the final frac bits)
        // so as to preserve precision. This would be a separate
        // internal multiplication possibly using 4 ints: w for
        // the whole and two for the frac parts of each operand

        // Shift y adjusting it for final fxp,
        // and round last bit
        //unsigned long mask = 1u << (FXP_lg2_l_mshift - 1);
        //unsigned long rbit = y & mask;
        //long m = (rbit > 0u)?
        //            (y >> FXP_lg2_l_mshift) + 1:
        //            (y >> FXP_lg2_l_mshift);

        unsigned long rbit = (y >> FXP_lg2_l_mshiftm1) & 1ul;
        long m = (y >> FXP_lg2_l_mshift) + rbit;
        int im = (int) m;
        // Return the complete logarithm (c + m) as fxp
        if (c < FXP_whole_min) {
                if ((c == FXP_whole_min_m1) && (im > 0)) {
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
        //printf("\txa:%X, xb:%X, ya:%X, yb:%X\n", xa, xb, ya, yb);
        xayb = xa * yb;
        yaxb = ya * xb;
        qr1 = xayb & FXP_RINT_MASK;
        qr2 = yaxb & FXP_RINT_MASK;
        qr3 = xbyb >> FXP_INT_BITS;
        //int rbit = ((unsigned int) (xbyb >> FXP_INT_BITS_M1)) & 1;
        qrsum = qr1 + qr2 + qr3; // + rbit;
        xaya = xa * ya;
        xbyb = xb * yb;
        ql1 = xayb >> FXP_INT_BITS;
        ql2 = yaxb >> FXP_INT_BITS;
        //rbit = ((unsigned int) (qrsum >> FXP_INT_BITS_M1)) & 1;
        ql3 = (qrsum >> FXP_INT_BITS); // + rbit;
        product = xaya + ql1 + ql2 + ql3;
        return product;

}

/*
 * pow2 calculation (2^fxp1) using the BKM (E-mode) algorithm
 * Analogous to implementation in bkm.c,
 * but here tailored for fxp's, and using longs.
 */
int fxp_pow2_l(int fxp1)
{
        if (fxp1 == FXP_UNDEF) return FXP_UNDEF;
        if (fxp1 == FXP_NEG_INF) return 0;
        if (fxp1 == FXP_POS_INF) return FXP_POS_INF;
        int w = fxp_get_whole_part(fxp1);
        int f = fxp_get_bin_frac(fxp1);
        //printf("pow2_l(%X)  w:%X  f:%X\n", fxp1, w, f);
        unsigned long pow2w, argument;
        if (fxp1 >= 0) {
                if (w >= FXP_whole_bits_m1) {
                        return FXP_POS_INF;
                }
                pow2w = FXP_one_l << w + 2;
                // Argument will be in [0, 1)
                argument = ((unsigned long) f) << FXP_lg2_l_mshift;
        } else {
                if (w <= FXP_INT_BITS_M1_NEG) {
                        return 0;
                }
                //pow2w = (long) (FXP_one >> (-w + 1));
                if (w < 0)
                        pow2w = FXP_one_l >> (-w - 1);
                else
                        pow2w = FXP_one_l << 1;
                // Notice argument is >= 0, in (0, 1]
                argument = ((unsigned long) (FXP_one + f)) \
                                << FXP_lg2_l_mshift;
        }

        printf("\nfxp:%d. %d(x%X), pow2w:%lX, argument:%lX\n", \
                    w, f, ((unsigned int) fxp1), pow2w, argument);

        // Watch out we are using two different shifted
        // configurations simultaneously, but independently:
        // - x aligned with FXP_BKM_ONE_L: 2 whole and 62 frac bits
        //   (we needed at least 2 whole bits for the BKM log algorithm).
        // - However argument, y and the BKM_LOGS array values
        //   have 1 whole and 63 frac bits for more accuracy
        unsigned long x = FXP_BKM_ONE_L, y = 0;
        //for (int k = 0; k < FXP_lg2_maxloops; k++) {
        for (int k = 0; k < FXP_INT_BITS; k++) {
                unsigned long const  z = y + FXP_BKM_LOGS_L[k];
                //printf("\tpow2_l iter. %d,  z:%lX,  x:%lX\n", \
                //            k, z, x);
                if (z <= argument) {
                        y = z;
                        x = x + (x >> k);
                        //printf("\t\tUpdating y:%lX  x:%lX\n", y, x);
                }
        }
        //unsigned int up2w = (unsigned int) pow2w;
        // Rounding x
        //unsigned int rbit = (unsigned int) ((x >> FXP_INT_BITS_M1) & 1l);
        //unsigned int ux = (unsigned int) (x >> FXP_INT_BITS) + rbit;
        //printf("up2w:%X  ux:%X\n", up2w, ux);
        //unsigned int muldist = mul_distrib(up2w, ux) << 2;
        //pow2w <<= 2;
        unsigned int muldist = mul_distrib_l(pow2w, x);
        //rbit = (x >> FXP_pow2_l_xshiftm1) & 1l;
        //x = (x >> FXP_pow2_l_xshift) + rbit;
        //unsigned long product = fxp_mul_l(pow2w, x);
        //int finalpow2 = (int) product; //+ rbit;
        //printf("product:%X\n", finalpow2);
        //printf("muldist:%X\n", muldist);
        //int pw = fxp_get_whole_part(finalpow2);
        //int pf = fxp_get_bin_frac(finalpow2);
        //printf("pow2() = %d. %d(/max frac)\n", pw, pf);
        //return finalpow2;
        return (int) muldist;
}
