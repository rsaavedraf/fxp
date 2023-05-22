
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
#include "fxp_conv.h"
#include "fxp_aux.h"
#include <stdio.h>
#include <stdlib.h>
#include "print_as_bits.h"
#include <assert.h>

//#define VERBOSE 1
//#define VERBOSE_BKM 1
#ifdef VERBOSE
#include "print_as_bits.h"
#endif

const unsigned long FXP_E_I64 = 0xADF85458A2BB4A9Bul;
const unsigned long FXP_PI_I64 = 0xC90FDAA22168C235ul;
const unsigned long FXP_LN_2_I64 = 0xB17217F7D1CF7C72ul;
const unsigned long FXP_LG10_2_I64 = 0x4D104D427DE7FD01ul;
const unsigned long FXP_LG2_E_I64 = 0xB8AA3B295C17F19Eul;
const unsigned long FXP_LG2_10_I64 = 0xD49A784BCD1B8B51ul;

const unsigned long ULONG_ALL_ONES = ~0ul;
const unsigned long ULONG_ALL_ONES_RS1 = ~0ul >> 1;

// For the BKM lg2 calculation when using longs.
// Values in this array are: a[k] = lg2(1 + 1/2^k) represented as
// unsigned long fxp's (8 bytes,) with 63 frac bits, and
// one magnitude whole bit. Notice, not a sign but a magnitude single
// whole bit, needed to be able to represent the first value in
// the table == 1.0. That value might be removed and handled
// implicitely later, for 1 extra bit of precision.
// 63 frac bits corresponds to a precision of 18.96 decimal digits
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
        0xB8AA3B2, 0x5C551D9, 0x2E2A8EC, 0x1715476,       // <---- * 40 entries
        //0xB8AA3B, 0x5C551D, 0x2E2A8E, 0x171547,
        //0xB8AA3, 0x5C551, 0x2E2A8, 0x17154,
        //0xB8AA, 0x5C55, 0x2E2A, 0x1715,         // <--- 52 entries
        //0xB8A, 0x5C5, 0x2E2, 0x171,
        //0xB8, 0x5C, 0x2E, 0x17,
        //0xB, 0x5, 0x2, 0x1,
        //0x0
// Starting with the row marked with the *, each entry is exactly
// a 4-bit right-shift of the value 4 positions earlier
};

// Auxiliary struct used internally for logarithms
typedef struct tuple_l {
        int ping;
        unsigned long pong;
} tuple_l;

// A super_fxp is not just a "bigger FXP," it also has its own
// fxp configuration (it's own number of whole vs. frac bits,)
// independent of the global fxp setting, always carrying this
// independent configuration within.
typedef struct super_fxp_l {
        int nwbits;             // <- Number of whole bits
        unsigned long number;   // <- The super fxp number
} super_fxp_l;

const super_fxp_l SFXP_ZERO_L = {0, 0ul};
const super_fxp_l SFXP_LG2_E_FACTOR_L = {1, FXP_LG2_E_I64};
const super_fxp_l SFXP_LG2_10_FACTOR_L = {2, FXP_LG2_10_I64};
const super_fxp_l SFXP_LN_2_FACTOR_L = {0, FXP_LN_2_I64};
const super_fxp_l SFXP_LG10_2_FACTOR_L = {0, FXP_LG10_2_I64};

static inline int fxp_nbits_l(unsigned long x)
{
        if (x == 0) return 0;
        return FXP_LONG_BITS - __builtin_clzl(x);
};

static inline unsigned long fxp_get_lshifted_frac_l(unsigned int fxp1)
{
        return ((unsigned long) (fxp1 & FXP_frac_mask)) << FXP_int_plus_whole_bits_m1;
}

static inline tuple_l tuple_l_from_fxp(int fxp1)
{
        tuple_l tup = {   fxp_get_whole_part(fxp1),
                          fxp_get_lshifted_frac_l(fxp1) };
        return tup;
}


static inline int super_fxp_l_get_nwbits(super_fxp_l x)
{
        return (x.nwbits < 0)? -x.nwbits: x.nwbits;
}

static inline unsigned long super_fxp_l_get_unumber(super_fxp_l x)
{
        return x.number;
}

static inline tuple_l negate_tuple_l(tuple_l tup)
{
        if (tup.pong == 0ul) {
                tuple_l negtup = { -tup.ping, 0ul };
                return negtup;
        } else {
                tuple_l negtup = { -tup.ping - 1, \
                                 FXP_BKM_A_ONE_L - tup.pong};
                return negtup;
        }
}

#ifdef VERBOSE
static inline int super_fxp_l_2_fxp(super_fxp_l x)
{
        if (x.nwbits >= 0) {
                int sfb = FXP_LONG_BITS - x.nwbits;
                int fbdiff = sfb - FXP_frac_bits;
                if (fbdiff >= 0) {
                        int rbit = (x.number >> (fbdiff - 1)) & 1u;
                        unsigned long num = (x.number >> fbdiff) + rbit;
                        unsigned long whole = num >> FXP_frac_bits;
                        if (whole > FXP_whole_max) return FXP_POS_INF;
                        return (int) num;
                } else {
                        unsigned long whole = x.number >> sfb;
                        if (whole > FXP_whole_max) return FXP_POS_INF;
                        int frac = (int) (x.number & (ULONG_ALL_ONES >> x.nwbits));
                        frac <<= -fbdiff;
                        return fxp_bin((int) whole, frac);
                }
        } else {
                int sfb = FXP_LONG_BITS + x.nwbits;
                int fbdiff = sfb - FXP_frac_bits;
                if (fbdiff >= 0) {
                        int rbit = (x.number >> (fbdiff - 1)) & 1u;
                        unsigned long num = (x.number >> fbdiff) + rbit;
                        unsigned long whole = num >> FXP_frac_bits;
                        if (whole > FXP_whole_max) return FXP_NEG_INF;
                        return -((int) num);
                } else {
                        unsigned long whole = x.number >> sfb;
                        if (whole > FXP_whole_max) return FXP_NEG_INF;
                        int frac = (int) (x.number & (ULONG_ALL_ONES >> x.nwbits));
                        frac <<= -fbdiff;
                        return -fxp_bin((int) whole, frac);
                }
        }
}

void print_super_fxp_l(super_fxp_l x)
{
        printf("Super_fxp_l: {x%X, x%lX}\n", x.nwbits, x.number);
        printf("Fxp equiv. : "); print_fxp(super_fxp_l_2_fxp(x)); printf("\n");
}

#endif

void print_tuple_l(char * msg, tuple_l tup)
{
        long double num = (long double) tup.ping \
                + ((long double) tup.pong) / FXP_BKM_A_ONE_L;
        printf("%s: ", msg);
        printf("{x%X, x%lX} == %.10Lf\n", \
                tup.ping, tup.pong, num);
}

static inline super_fxp_l super_fxp_l_create(unsigned int nwbits, \
                                             unsigned long num)
{
        nwbits = (nwbits > SFXP_MAX_WBITS)? SFXP_MAX_WBITS: nwbits;
        super_fxp_l x = { nwbits, num };
        return x;
}

static inline super_fxp_l super_fxp_l_from_tuple_l(tuple_l tup)
{
        unsigned long nwb = fxp_nbits(tup.ping);
        unsigned long number;
        if (nwb == 0) {
                number = tup.pong << 1;
        } else if (nwb == 1) {
                number = FXP_BKM_A_ONE_L | tup.pong;
        } else {
                number = (((unsigned long) tup.ping) << (FXP_LONG_BITS - nwb)) \
                                    | (tup.pong >> (nwb - 1));
        }
        return super_fxp_l_create(nwb, number);
}

// Commenting out these functions which aren't being used right now,

/*
static inline long super_fxp_l_get_signed_num(super_fxp_l x)
{
        return (x.nwbits < 0)? -((long) x.number): (long) x.number;
}

static inline unsigned long sfxp_l_get_poswhole_ul(super_fxp_l x)
{
        return x.number >> (FXP_LONG_BITS - x.nwbits);
}

static inline unsigned long super_fxp_l_get_lshifted_frac(super_fxp_l x)
{
        return (x.nwbits >= 0)? x.number << x.nwbits:
                                x.number << -x.nwbits;
}

static inline super_fxp_l super_fxp_l_negate(super_fxp_l x)
{
        super_fxp_l y = { -x.nwbits, x.number };
        return y;
}

static inline super_fxp_l posfxp_x_possuper_fxp_l(int x, super_fxp_l c)
{
        int clz_m1 = __builtin_clz(x) - 1;
        unsigned int shifted_x = x << clz_m1;
        unsigned long product = dmul_ulong_x_uint(c.number, x);
        int nwb = c.nwbits + FXP_whole_bits;
        super_fxp_l super_product = super_fxp_l_create_unsigned(nwb, product);
        return super_product;
}

static inline super_fxp_l get_fxp_x_super_fxp_l(int x, super_fxp_l c)
{
        int sx = (x >= 0);
        int sc = (c.nwbits >= 0);
        if (sx) {
                if (sc) {
                        // Both NON-negative
                        return posfxp_x_possuper_fxp_l(x, c);
                } else {
                        // Only c negative
                        return super_fxp_l_negate( \
                                    posfxp_x_possuper_fxp_l(x, \
                                            super_fxp_l_negate(c)));
                }
        } else {
                if (sc) {
                        // Only x negative
                        return super_fxp_l_negate( \
                                    posfxp_x_possuper_fxp_l(-x, c));
                } else {
                        // Both negative
                        return super_fxp_l_negate( \
                                    posfxp_x_possuper_fxp_l(-x, \
                                            super_fxp_l_negate(c)));
                }
        }
}
*/

/*
 * r-shift an unsigned long by shift, rounding its last bit
 */
static inline unsigned long rshift_ulong_rounding( \
                                    unsigned long n,
                                    unsigned int shift)
{
        unsigned long rbit = (shift == 0)? 0: (n >> (shift - 1)) & 1ul;
        return (n >> shift) + rbit;
}

static inline unsigned int rshift_ulong_into_uint_rounding( \
                                    unsigned long n)
{
        unsigned long rbit = (n >> FXP_INT_BITS_M1) & 1ul;
        return (unsigned int) ((n >> FXP_INT_BITS) + rbit);
}

/*
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
        int rbit = (product >> FXP_frac_bits_m1) & 1u;
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
        if (fxp1 <= FXP_almost1) {
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
 * struct (tuple), both at their maximum int-sized precision:
 * The characteristic as int with 0 frac bits,
 * the mantissa as unsigned int with 31 frac bits (exact same
 * alignment as hi int in the BKM array values)
 * Uses the BKM L-Mode algorithm (requires table of pre-calculated
 * values) and longs.
 *
 * For more details on BKM:
 * https://en.wikipedia.org/wiki/BKM_algorithm
 */
static inline tuple_l fxp_get_lg2_as_tuple_l(int fxp1, \
                                const int MAX_LOOPS)
{
        tuple_l result;
        // Assumes fxp1 is for sure > 0
        int clz = __builtin_clz((unsigned int) fxp1);
        // clz > 0 since at least sign bit == 0
        int nbx = FXP_INT_BITS - clz;
        // Assign the characteristic
        result.ping = \
                ((fxp1 <= FXP_almost1) || (fxp1 >= FXP_two))? \
                        (nbx - FXP_frac_bits - 1): 0;
        #ifdef VERBOSE
        printf("\nlg2_as_tuple: fxp1: x%X  clz: %d,  nbx: %d\n", fxp1, clz, nbx);
        printf("lg2_as_tuple: characteristic is : %d\n", result.ping);
        #endif

        // Here we have already calculated the log characteristic c
        // Now lshifting the original number to have it as unsigned
        // long fxp with 2 whole bits (2 needed because z can
        // intermitently get a value > 2)
        // Watch out we are using two different shifted
        // configurations simultaneously, but independently:
        // - "A" alignment (A for Array and argument)
        //   BKM_LOGS Array values, argument, and the
        //   mantissa all have 1 unsigned whole and 63 frac bits
        //   for more accuracy
        // - "X" alignment:
        //   x aligned with FXP_BKM_X_ONE_L: 2 whole and 62 frac bits
        //   (again at least 2 whole bits for the BKM log algorithm,
        //   given that z can occasionally exceed 2)
        unsigned long argument = ((unsigned long) fxp1) << \
                            (clz + FXP_INT_BITS_M1);
        unsigned long x = FXP_BKM_X_ONE_L;
        unsigned long z, xs, y;
        // The mantissa value will remain in the range of lg2(x),
        // x in [1, 2), meaning mantissa always in [0, 1)
        y = 0ul;
        for (int shift = 1; shift <= MAX_LOOPS; shift++) {
                #ifdef VERBOSE_BKM
                printf("bkm lg2 shift:%2d - ", shift);
                #endif
                xs = (x >> shift);
                z = x + xs;
                if (z <= argument) {
                        x = z;
                        y += FXP_BKM_LOGS_L[shift];
                        #ifdef VERBOSE_BKM
                        printf("Udated x:x%lX  Updated y:x%lX", x, y);
                        #endif
                }
                #ifdef VERBOSE_BKM
                printf("\n");
                #endif
        }
        //result.pong = rshift_ulong_into_uint_rounding(y);
        result.pong = y;
        #ifdef VERBOSE
        printf("lg2_as_tuple_l: mantissa is: x%lX\n", result.pong);
        #endif
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
        tuple_l tup = fxp_get_lg2_as_tuple_l(fxp1, FXP_lg2_maxloops);
        // Round and shift the mantissa
        int rbit = (tup.pong >> FXP_int_plus_whole_bits_m2) & 1u;
        int m = (int) (tup.pong >> FXP_int_plus_whole_bits_m1) + rbit;
        #ifdef VERBOSE
        printf("\nlg2 as tuple result: {%d, %lu} == {x%X,x%lX}\n", \
                        tup.ping, tup.pong, tup.ping, tup.pong);
        printf("Final mantissa: x%X  (rounding bit was %d)\n", m, (int) rbit);
        #endif
        // Return the complete logarithm (charact + mantissa) as fxp
        if (tup.ping < FXP_whole_min) {
                if ((tup.ping == FXP_whole_min_m1) \
                        && (m > 0ul)) {
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
        return (tup.ping >= 0)?
                        ((tup.ping << FXP_frac_bits) + m):
                        (-(-tup.ping << FXP_frac_bits) + m);
}

/*
 * Analogougs to distributive multiplication done in
 * ulongy_from_dmul in ulongy.c, but here using longs and
 * returning only the rounded higher long from the result.
 * Identifying xa, xb, ya and yb as the inner "uints"
 * (e.g. half ulongs) of the unsigned long arguments
 *      x == (xa << INT_BITS) | xb
 *      y == (ya << INT_BITS) | yb
 * Calculates the product as:
 *      x * y = (xa * ya) + ql1 + ql2 + ql3, where
 *      ql1 = lint of xa * yb
 *      ql2 = lint of ya * xb
 *      ql3 = lint of qrsum
 * and qrsum = qr1 + qr2 + qr3, where
 *      qr1 = rint of xa * yb
 *      qr2 = rint of ya * xb
 *      qr3 = lint of xb * yb
 */
unsigned long dmul_ulongs(unsigned long x, unsigned long y)
{
        unsigned long xa, xb, ya, yb, xaya, xayb, yaxb, xbyb;
        unsigned long qr1, qr2, qr3, qrsum, ql1, ql2, ql3, product;
        xa = (x >> FXP_INT_BITS);
        xb = (x & FXP_RINT_MASK);
        ya = (y >> FXP_INT_BITS);
        yb = (y & FXP_RINT_MASK);
        xayb = xa * yb;
        yaxb = ya * xb;
        xaya = xa * ya;
        xbyb = xb * yb;
        qr1 = (xayb & FXP_RINT_MASK);
        qr2 = (yaxb & FXP_RINT_MASK);
        qr3 = (xbyb >> FXP_INT_BITS);
        unsigned int rbit1 = ((xbyb >> FXP_INT_BITS_M1) & 1u);
        qrsum = qr1 + qr2 + qr3 + rbit1;
        ql1 = (xayb >> FXP_INT_BITS);
        ql2 = (yaxb >> FXP_INT_BITS);
        unsigned int rbit2 = ((qrsum >> FXP_INT_BITS_M1) & 1u);
        ql3 = (qrsum >> FXP_INT_BITS);
        product = xaya + ql1 + ql2 + ql3 + rbit2;
        return product;
}

/*
 * Equivalent to dmul_ulong, but simplified because yb == zero
 */
unsigned long dmul_ulong_x_uint(unsigned long x, unsigned int ya)
{
        unsigned long xa, xb, xaya, yaxb;
        unsigned long ql2, qr2, ql3, product;
        xa = (x >> FXP_INT_BITS);
        xb = (x & FXP_RINT_MASK);
        yaxb = ya * xb;
        xaya = xa * ya;
        qr2 = (yaxb & FXP_RINT_MASK);
        ql2 = (yaxb >> FXP_INT_BITS);
        ql3 = (qr2 >> FXP_INT_BITS);
        unsigned int rbit2 = ((qr2 >> FXP_INT_BITS_M1) & 1u);
        product = xaya + ql2 + ql3 + rbit2;
        return product;
}

static inline int neg_lg2_x_factor_l(tuple_l tup, super_fxp_l factor)
{
        if ((tup.ping < FXP_whole_min_m1) \
                || ((tup.ping == FXP_whole_min_m1) \
                    && (tup.pong == 0ul))) {
                return FXP_NEG_INF;
        }
        unsigned int posc = (-tup.ping);
        int shift_for_c = __builtin_clz(posc) - 1;
        unsigned int shifted_posc = posc << shift_for_c;
        unsigned long fval = super_fxp_l_get_unumber(factor);
        long cxf = -dmul_ulong_x_uint(fval, shifted_posc);

        unsigned long mxf = dmul_ulongs(fval, tup.pong);
        int shift_for_m = FXP_INT_BITS_M1 - shift_for_c;
        long shifted_mxf = rshift_ulong_rounding(mxf, shift_for_m);
        long sum = cxf + shifted_mxf;

                #ifdef VERBOSE
                printf("\tShifted c (by %d): x%X,  b", shift_for_c, shifted_posc);
                print_uint_as_bin(shifted_posc); printf("\n");
                printf("\tcxf: x%16lX\n\tb", cxf);
                print_ulong_as_bin(cxf); printf("\n");
                printf("\tmxf: x%16lX\n\tb", mxf);
                print_ulong_as_bin(mxf); printf("\n");
                printf("\tR-shifted mxf (by %d): x%lX\n\tb", shift_for_m, shifted_mxf);
                print_ulong_as_bin((unsigned long) shifted_mxf); printf("\n");
                printf("\tSum: x%lX\n\tb", sum);
                print_ulong_as_bin((unsigned long) sum); printf("\n");
                #endif

        int final_rshift = FXP_whole_bits + shift_for_c;
        // Finally round and shift for current fxp configuration
        long psum = -sum;
        int final_lg = (int) -(rshift_ulong_rounding(psum, final_rshift));

                #ifdef VERBOSE
                printf("\tfinal_rshift: %d\n", final_rshift);
                printf("\tFinal lg: x%X,  b", final_lg);
                print_int_as_bin(final_lg, 0); printf("\n");
                #endif

        return final_lg;
}

static inline int pos_lg2_x_factor_l(tuple_l tup, super_fxp_l factor)
{
        int shift_for_c;
        unsigned int shifted_c;
        long cxf;
        unsigned long fval = super_fxp_l_get_unumber(factor);
        if (tup.ping > 0) {
                shift_for_c = __builtin_clz(tup.ping) - 1;
                shifted_c = tup.ping << shift_for_c;
                cxf = dmul_ulong_x_uint(fval, shifted_c);
        } else {
                shift_for_c = FXP_INT_BITS_M1;
                shifted_c = 0;
                cxf = 0l;
        }

        unsigned long mxf = dmul_ulongs(fval, tup.pong);
        int shift_for_m = FXP_INT_BITS_M1 - shift_for_c;
        long shifted_mxf = rshift_ulong_rounding(mxf, shift_for_m);
        long sum = cxf + shifted_mxf;

                #ifdef VERBOSE
                printf("\tShifted c (by %d): x%X,  b", shift_for_c, shifted_c);
                print_uint_as_bin(shifted_c); printf("\n");
                printf("\tcxf: x%16lX\n\tb", cxf);
                print_ulong_as_bin(cxf); printf("\n");
                printf("\tmxf: x%16lX\n\tb", mxf);
                print_ulong_as_bin(mxf); printf("\n");
                printf("\tR-shifted mxf (by %d): x%lX\n\tb", shift_for_m, shifted_mxf);
                print_ulong_as_bin((unsigned long) shifted_mxf); printf("\n");
                printf("\tSum: x%lX\n\tb", sum);
                print_ulong_as_bin((unsigned long) sum); printf("\n");
                #endif

        int final_rshift = FXP_whole_bits + shift_for_c;
        // Finally round and shift for current fxp configuration
        int final_lg = (int) rshift_ulong_rounding(sum, final_rshift);

                #ifdef VERBOSE
                printf("\tfinal_rshift: %d\n", final_rshift);
                printf("\tFinal lg: x%X,  b", final_lg);
                print_int_as_bin(final_lg, 0); printf("\n");
                #endif

        return final_lg;
}

/*
 * Calculates lg2 and then multiplies by the given factor
 * using longs.
 */
static inline int lg2_x_factor_l(int fxp1, super_fxp_l factor)
{
        tuple_l tup = fxp_get_lg2_as_tuple_l(fxp1, FXP_LOGX_LOOPS);

        // The lg2(fxp1) is equals to c + m:
        // c == characteristic (tup.ping),
        // m == mantissa (tup.pong)
        // We multiply each part by FACTOR separately keeping maximum
        // possible precision (l-shifting them first to the max,)
        // then we align the results properly to sum them together

                #ifdef VERBOSE
                printf("\nlg2_x_factor_l:\n");
                printf("\tlg2 c : %d (x%X,  b", tup.ping, tup.ping);
                print_int_as_bin(tup.ping, 0); printf(")\n");
                printf("\tlg2 m : x%16lX\n", tup.pong);
                printf("\tFACTOR: ");
                print_super_fxp_l(factor);
                #endif

        return (tup.ping < 0)?
                neg_lg2_x_factor_l(tup, factor):
                pos_lg2_x_factor_l(tup, factor);
}

/*
 * Implementation of ln_l using lg2_l
 */
int fxp_ln_l(int fxp1)
{
        if (fxp1 < 0) return FXP_UNDEF;
        if (fxp1 == 0) return FXP_NEG_INF;
        if (fxp1 == FXP_POS_INF) return FXP_POS_INF;
        return lg2_x_factor_l(fxp1, SFXP_LN_2_FACTOR_L);
}

/*
 * Implementation of lg10_l using lg2_l
 */
int fxp_lg10_l(int fxp1) {
        if (fxp1 < 0) return FXP_UNDEF;
        if (fxp1 == 0) return FXP_NEG_INF;
        if (fxp1 == FXP_POS_INF) return FXP_POS_INF;
        return lg2_x_factor_l(fxp1, SFXP_LG10_2_FACTOR_L);
}

#ifdef VERBOSE
static long double print_ulong_as_ld(unsigned long x, int nfbits)
{
        unsigned long long twopower = 1llu << nfbits;
        long double ldfrac = 0.0L;
        unsigned long frac = x;
        while (frac > 0) {
                if (frac & 1u) ldfrac += ((long double) 1.0L) / twopower;
                frac = frac >> 1;
                twopower = twopower >> 1;
        }
        return ldfrac;
}
#endif

/*
 * Implementation of the BKM E-Mode algorithm for 2^argument calculation.
 * Analogous to the inner BKM E-Mode loop inside the my_pow2_bkm function
 * in bkm.c, but here using longs. Just as in the BKM function for lg2,
 * here the argument must come BKM A-aligned with the values in the log
 * table, that is: 1 unsigned whole, and 63 frac bits
 */
static inline unsigned long fxp_bkm_emode_l(unsigned long argument, \
                                            const int MAX_LOOPS)
{
        #ifdef VERBOSE_BKM
        printf("\nbkm_emode_l running with argument: x%lX\n", argument);
        print_ulong_as_bin(argument); printf("\n");
        #endif
        //unsigned long argument = ((unsigned long) a) << FXP_INT_BITS;
        unsigned long x = FXP_BKM_X_ONE_L, y = 0;
        // TODO: check optimizing here possibly starting at k = 1
        for (int k = 0; k < MAX_LOOPS; k++) {
                unsigned long const z = y + FXP_BKM_LOGS_L[k];
                #ifdef VERBOSE_BKM
                printf("k:%d,  z:%lX,  z as ld:%.19Lf\n", \
                        k, z, print_ulong_as_ld(z, 63));
                #endif
                if (z <= argument) {
                        y = z;
                        x = x + (x >> k);
                        #ifdef VERBOSE_BKM
                        printf("\tUpdating y (%lX) und x (%lX)\n", y, x);
                        printf("\tx as long double: %.19Lf\n", \
                                        print_ulong_as_ld(x, 62));
                        #endif
                }
        }
        #ifdef VERBOSE_BKM
        printf("bkm_emode_l returning x = %lX  (%.19Lf)\n", \
                        x, print_ulong_as_ld(x, 62));
        #endif
        return x;
}

/*
 * pow2 calculation (2^fxp1) of a NON-negative argument,
 * using the BKM (E-mode) algorithm and using longs.
 * Argument comes as a {whole,frac} tuple, where the
 * frac component must already have the same alignment
 * as the BKM array values (1 whole bit).
 * Input value for 2^x here is x = whole + frac
 */
static inline int fxp_pow2_wpos_tuple_l(tuple_l tup, \
                                        const int MAX_LOOPS)
{
        #ifdef VERBOSE
        print_tuple_l("\npow2_wpos_tuple_l for ", tup);
        #endif
        if (tup.ping >= FXP_whole_bits_m1) return FXP_POS_INF;
        if ((tup.ping == 0) && (tup.pong == 0ul)) return FXP_one;
        // Argument tup.pong will be in [0, 1)
        unsigned long x = fxp_bkm_emode_l(tup.pong, MAX_LOOPS);
        // When one of the operands is a shifted 1, we don't really
        // need to call the (expensive) dmul operation, since the
        // result will be identical to the other operand anyway,
        // just shifted appropriately
        int shift = FXP_int_plus_whole_bits_m2 - tup.ping;
        int shifted = (int) rshift_ulong_rounding(x, shift);
        #ifdef VERBOSE
        printf("x from bkm_e: x%lX\n", x);
        print_ulong_as_bin(x); printf("\n");
        printf("shift is %d\n", shift);
        print_ulong_as_bin(x); printf("\n");
        print_uint_as_bin(shifted); printf("\n");
        #endif
        return (shifted != FXP_UNDEF)? shifted: FXP_POS_INF;
}

/*
 * pow2 calculation (2^fxp1) of a negative argument,
 * using the BKM (E-mode) algorithm and longs
 * Argument comes as a {whole,frac} tuple, where
 * both whole and frac are positive, both corresponding
 * to minus the corresponding parts of the original
 * negative fxp (see how the argument gets built
 * in fxp_pow2_l.) Frac component must already have the same
 * alignment as the BKM array values (1 whole bit).
 *
 * TODO:
 * These functions for fxp_pow2_l were written with performance
 * in mind but we have an abstraction issue here:
 * we can't pass in just any tuple as argument.
 * Whenever we need to call fxp_pow2_wpos_tuple_l or
 * fxp_pow2_wneg_tuple_l, we will need to know what the
 * exact sign and alignment requirements for the members in
 * those tuples are, and build them accordingly and precisely.
 * A tuple can hold values matching these requirements, but
 * these requirements really transcend the tuple abstraction.
 * In any case, consider some refactoring, maybe even
 * eliminating tuples altogether? Using only ulongys/super_fxp's?
 * See also long comment in sqrt_l
 */
static inline int fxp_pow2_wneg_tuple_l(tuple_l tup, \
                                        const int MAX_LOOPS)
{
        #ifdef VERBOSE
        print_tuple_l("\npow2_wneg_tuple_l for (-)", tup);
        #endif
        if (tup.ping > FXP_frac_bits) return 0;
        if ((tup.ping == 0) && (tup.pong == 0ul)) return FXP_one;
        // Notice argument a for bkm_emode will be in (0, 1]
        unsigned long a = FXP_BKM_A_ONE_L - tup.pong;
        unsigned long x = fxp_bkm_emode_l(a, MAX_LOOPS);
        int shift = FXP_int_plus_whole_bits_m1 + tup.ping;
        int shifted = (int) rshift_ulong_rounding(x, shift);
        #ifdef VERBOSE
        printf("x from bkm_e: x%lX\n", x);
        printf("shift is %d\n", shift);
        print_ulong_as_bin(x); printf("\n");
        print_uint_as_bin(shifted); printf("\n");
        #endif
        return (shifted != FXP_UNDEF)? shifted: FXP_POS_INF;
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
        #ifdef VERBOSE
        printf("\n\npow2_l running for x%X (int = %d)\n", fxp1, fxp1);
        #endif
        if (fxp1 >= 0) {
                tuple_l tup = tuple_l_from_fxp(fxp1);
                return fxp_pow2_wpos_tuple_l( tup, FXP_POWX_LOOPS );
        } else {
                tuple_l ntup = tuple_l_from_fxp(-fxp1);
                return fxp_pow2_wneg_tuple_l( ntup, FXP_POWX_LOOPS );
        }
}

// Internal function to calculate the product of
// x times C and return it as an {int, ulong} tuple
// with the ulong frac part having the same alignment
// as the BKM array values (1 whole bit)
static inline tuple_l get_xc_as_tuple_l(int x, super_fxp_l factorc)
{
        unsigned long fc, wfc, ffc, ffmask;
        unsigned int sf;

        // First calculating frac(x) * C
        // fraction part of x l-shifted so as
        // to leave just 1 whole bit
        unsigned long cval = super_fxp_l_get_unumber(factorc);
        unsigned int c_nwbits = super_fxp_l_get_nwbits(factorc);
        sf = fxp_get_bin_frac(x) << FXP_whole_bits_m1;
        // sf times the factor C
        fc = dmul_ulong_x_uint(cval, sf);
        // whole and frac parts of that product fxC
        ffmask = ULONG_ALL_ONES_RS1 >> c_nwbits;
        wfc = fc >> (FXP_LONG_BITS_M1 - c_nwbits);
        ffc = fc & ffmask;

                #ifdef VERBOSE
                printf("\n\tf (1 whole bit) and factor C (%d whole bits):\n\t", \
                                c_nwbits);
                print_uint_as_bin(sf);
                printf(" (%LE)\n\t", ((long double) sf / (1ul << FXP_INT_BITS_M1)));
                print_ulong_as_bin(cval);
                printf(" (%LE)\n\t", ((long double) cval) \
                                        / (1ul << FXP_LONG_BITS - c_nwbits));
                printf("Product fxC (%d whole bits) is:\n\t", c_nwbits + 1);
                print_ulong_as_bin(fc);
                printf(" (%LE)\n\t", ((long double) fc) \
                                        / (1ul << FXP_LONG_BITS_M1 - c_nwbits));
                printf("w_fxC: %lu\n\t", wfc);
                printf("f_fxC (%d whole bits):\n\t", c_nwbits + 1);
                print_ulong_as_bin(ffc);
                printf(" (%LE)\n\t", ((long double) ffc) \
                                        / (1ul << FXP_LONG_BITS_M1 - c_nwbits));
                #endif

        // Now calculating whole(x) * factorc
        unsigned long wc, wwc, wfmask, fwc;
        int wx_nbits, wx_clz, w_margin;
        unsigned int wx, swx;
        wx = fxp_get_whole_part(x);
        wx_clz = (wx == 0)? FXP_INT_BITS: __builtin_clz(wx);
        wx_nbits = FXP_INT_BITS - wx_clz;
        // swx is whole(x) l-shifted all the way
        swx = wx << wx_clz;
        // times factorc
        wc = dmul_ulong_x_uint(cval, swx);

        // Whole and frac parts of that product wxC
        w_margin = wx_nbits + c_nwbits;
        wfmask = ULONG_ALL_ONES >> w_margin;
        wwc = wc >> (FXP_LONG_BITS - w_margin);
        fwc = wc & wfmask;

                #ifdef VERBOSE
                printf("\n\tw (%d whole bits) and factor C (%d whole bits):\n\t", \
                            wx_nbits, c_nwbits);
                print_uint_as_bin(swx);
                printf(" (%LE)\n\t", ((long double) swx / (1ul << wx_clz)));
                print_ulong_as_bin(cval);
                printf(" (%LE)\n\t", ((long double) cval) \
                                        / (1ul << FXP_LONG_BITS - c_nwbits));
                printf("Product wxC (%d whole bits) is:\n\t", w_margin);
                print_ulong_as_bin(wc);
                printf(" (%LE)\n\t", ((long double) wc) \
                                        / (1ul << (FXP_LONG_BITS - w_margin)));

                printf("w_wxC: %lu\n\t", wwc);
                printf("f_wxC (%d whole bits):\n\t", w_margin);
                print_ulong_as_bin(fwc);
                printf(" (%LE)\n\t", ((long double) fwc) \
                                        / (1ul << FXP_LONG_BITS - w_margin));
                #endif

        // Left-align both f_fxC and f_wxC leaving 1 whole bit
        // in order to add them up
        ffc <<= c_nwbits;
        fwc <<= w_margin - 1;

        // Sum of frac parts
        unsigned long fplusf = ffc + fwc;

        // Carry over from frac sum into whole part
        unsigned long w_fplusf = fplusf >> FXP_LONG_BITS_M1;

        int final_whole = (int) wwc + wfc + w_fplusf;
        unsigned long final_frac = fplusf & ULONG_ALL_ONES_RS1;

        // Tuple to return with final whole and frac parts
        tuple_l xc = { final_whole, final_frac };

                #ifdef VERBOSE
                printf("\n\tfplusf:\n\t");
                print_ulong_as_bin(fplusf);
                printf("\n\tw_fplusf is: %lu", w_fplusf);
                printf(" (%u)\n\t", xc.ping);
                printf("final fsum:\n\t");
                print_ulong_as_bin(xc.pong);
                printf(" (%LE)\n", ((long double) xc.pong) \
                                        / (1ul << FXP_LONG_BITS_M1));
                #endif

        return xc;
}

// Calculate the pow2 of a NON-negative argument x
// multiplied by a factor c
static inline int fxp_pow2_pos_arg_xfactor_l(int x, \
                                        super_fxp_l factorc, \
                                        const int MAX_LOOPS)
{
        tuple_l xc = get_xc_as_tuple_l(x, factorc);
                #ifdef VERBOSE
                print_tuple_l("Pos tuple is: ", xc);
                #endif
        return fxp_pow2_wpos_tuple_l( xc, MAX_LOOPS );
}

// Calculate the pow2 of a negative argument x
// multiplied by a factor c
static inline int fxp_pow2_neg_arg_xfactor_l(int x, \
                                        super_fxp_l factorc, \
                                        const int MAX_LOOPS)
{
        tuple_l xc = get_xc_as_tuple_l(x, factorc);
                #ifdef VERBOSE
                print_tuple_l("Neg tuple is: ", xc);
                #endif
        return fxp_pow2_wneg_tuple_l( xc, MAX_LOOPS );
}

// Implementation of exp_l(x) using pow2_l():
// e^x  == 2^( lg2(e^x) )
//      == 2^( x * lg2(e) )
int fxp_exp_l(int fxp1)
{
        if (fxp1 == FXP_UNDEF) return FXP_UNDEF;
        if (fxp1 == FXP_NEG_INF) return 0;
        if (fxp1 == FXP_POS_INF) return FXP_POS_INF;
        return (fxp1 >= 0)?
                fxp_pow2_pos_arg_xfactor_l(fxp1, \
                                    SFXP_LG2_E_FACTOR_L, \
                                    FXP_POWX_LOOPS):
                fxp_pow2_neg_arg_xfactor_l( -fxp1, \
                                    SFXP_LG2_E_FACTOR_L, \
                                    FXP_POWX_LOOPS);
}

// Implementation of pow10_l(x) pow2_l():
// 10^x == 2^( lg2(10^x) )
//      == 2^( x * lg2(10) )
int fxp_pow10_l(int fxp1)
{
        if (fxp1 == FXP_UNDEF) return FXP_UNDEF;
        if (fxp1 == FXP_NEG_INF) return 0;
        if (fxp1 == FXP_POS_INF) return FXP_POS_INF;
        return (fxp1 >= 0)?
                fxp_pow2_pos_arg_xfactor_l(fxp1, \
                                    SFXP_LG2_10_FACTOR_L, \
                                    FXP_POWX_LOOPS):
                fxp_pow2_neg_arg_xfactor_l( -fxp1, \
                                    SFXP_LG2_10_FACTOR_L, \
                                    FXP_POWX_LOOPS);
}

// Square root implementation as:
// sqrt(x) = 2^( 0.5 * lg2(x) )
int fxp_sqrt_l(int fxp1)
{
        if (fxp1 < 0) return FXP_UNDEF;
        if (fxp1 == 0) return 0;
        if (fxp1 == FXP_POS_INF) return FXP_POS_INF;
        // First get the lg2 of the argument
        tuple_l tup = fxp_get_lg2_as_tuple_l(fxp1, FXP_SQRT_LOOPS);
        // TODO:
        // Div by 2 on a normal fxp or a super_fxp would be a simple
        // matter of r-shifting 1 bit + adding rounding bit.
        // But converting to long/super_fxp, r-shifting there,
        // then back, seems quite expensive.
        // Later will try not using tuples at all, only super_fxp's
        // comparing possible performance and simplicity gains/losses,
        // while also keeping an eye on potential precision loss
        unsigned long notrounded, rbit = tup.pong & 1ul; // rounding bit for frac
        if (tup.ping >= 0) {
                // Halve the non-negative tuple before passing it to pow2_wpos
                notrounded = (tup.ping % 2)? \
                                // odd characteristic, add 0.5 to the rshifted mantissa
                                FXP_BKM_A_POINT5_L + (tup.pong >> 1): \
                                // even characteristic, simply r-shift the mantissa
                                tup.pong >> 1;
                tup.pong = notrounded + rbit;
                tup.ping = (tup.ping >> 1) + (tup.pong < notrounded);
                #ifdef VERBOSE
                print_tuple_l("sqrt_l: halved pos tuple for pow2: ", tup);
                #endif
                return fxp_pow2_wpos_tuple_l(tup, FXP_SQRT_LOOPS);
        } else {
                // Negative tuple: negate and halve before passing to pow2_wneg
                tuple_l ntup = negate_tuple_l(tup);
                notrounded = (ntup.ping % 2)? \
                                FXP_BKM_A_POINT5_L + (ntup.pong >> 1): \
                                ntup.pong >> 1;
                ntup.pong = notrounded + rbit;
                ntup.ping = (ntup.ping >> 1) + (ntup.pong < notrounded);
                #ifdef VERBOSE
                print_tuple_l("\nsqrt_l: 'halved neg' tuple for pow2: ", ntup);
                #endif
                return fxp_pow2_wneg_tuple_l(ntup, FXP_SQRT_LOOPS);
        }
}

/*
Implementation of powxy_l(x) using pow2_l() and lg2_l():
     x^y == 2^( y * lg2(x) )
Expected output given the combination of inputs x and y:

  x     y:  Und. -INF [MIN,-1)  -1  (-1,0)   0   (0,1)   1  (1,MAX] +INF
< 0         Und.  Und.   Und.   Und.  Und.  Und.  Und.  Und.  Und.   Und.
  0         Und.  +INF  +INF   +INF  +INF   Und.   0     0     0      0
(0,1)       Und.   0      *      *     *     1     *     *     *      0
  1         Und.  Und.    1      1     1     1     1     1     1     Und.
(1,MAX]     Und.   0      *      *     *     1     *     *     *     +INF
 +INF       Und.   0      0      0     0    Und.  +INF  +INF  +INF   +INF

The * means doing the actual calculations for x and y
*/
int fxp_powxy_l(int x, int y)
{
        if ((x < 0) || (y == FXP_UNDEF)) {
                return FXP_UNDEF;
        }
        if (x == 0) {
                return (y < 0)? FXP_POS_INF : \
                            (y == 0)? FXP_UNDEF: 0;
        }
        if (x == FXP_POS_INF) {
                return (y < 0)? 0: \
                            (y == 0)? FXP_UNDEF: FXP_POS_INF;
        }
        if (y == FXP_NEG_INF) {
                return (x == FXP_one)? FXP_UNDEF: 0;
        }
        if (y == FXP_POS_INF) {
                return (x < FXP_one)? 0: \
                            (x == FXP_one)? FXP_UNDEF: FXP_POS_INF;
        }
        // First get the lg2 of x
        tuple_l tup = fxp_get_lg2_as_tuple_l(x, FXP_POWXY_LG_LOOPS);
        #ifdef VERBOSE
        print_tuple_l("powxy_l: lg2 tuple", tup);
        #endif
        // Create superfxp from the tuple and return
        // appropriately signed 2^( y*superfxp )
        if (tup.ping >= 0) {
                super_fxp_l superposlog = super_fxp_l_from_tuple_l(tup);
                // Now compute number of whole bits we will get in the product
                if (y >= 0) {
                        int pnwbits = superposlog.nwbits + fxp_nbits(y) \
                                        - FXP_frac_bits;
                        #ifdef VERBOSE
                        printf("\ncase 0: +ping +y, product would have %d whole bits\n", pnwbits);
                        printf("Whole bit counts:  factor: %d,  y: %d\n", \
                                        superposlog.nwbits, \
                                        fxp_nbits(y) - FXP_frac_bits);
                        #endif
                        return (pnwbits >= FXP_whole_bits)? \
                                FXP_POS_INF: \
                                fxp_pow2_pos_arg_xfactor_l(y, \
                                                    superposlog, \
                                                    FXP_POWXY_POW_LOOPS);
                } else {
                        int posy = -y;
                        int pnwbits = superposlog.nwbits + fxp_nbits(posy) \
                                        - FXP_frac_bits;
                        #ifdef VERBOSE
                        printf("\ncase 1: +ping -y, product would have %d whole bits\n", pnwbits);
                        printf("Whole bit counts:  factor: %d,  y: %d\n", \
                                        superposlog.nwbits, \
                                        fxp_nbits(posy) - FXP_frac_bits);
                        #endif
                        return (pnwbits >= FXP_whole_bits)? \
                                0: \
                                fxp_pow2_neg_arg_xfactor_l(posy, \
                                                    superposlog, \
                                                    FXP_POWXY_POW_LOOPS);
                }
        } else {
                super_fxp_l superposlog = super_fxp_l_from_tuple_l( \
                                                    negate_tuple_l(tup));
                if (y >= 0) {
                        #ifdef VERBOSE
                        int pnwbits = superposlog.nwbits + fxp_nbits(y) \
                                        - FXP_frac_bits;
                        printf("\ncase 2: -ping +y, product would have %d whole bits\n", pnwbits);
                        printf("Whole bit counts:  factor: %d,  y: %d\n", \
                                        superposlog.nwbits, fxp_nbits(y) - FXP_frac_bits);
                        #endif
                        return fxp_pow2_neg_arg_xfactor_l(y, \
                                                    superposlog, \
                                                    FXP_POWXY_POW_LOOPS);
                } else {
                        int posy = -y;
                        int pnwbits = superposlog.nwbits + fxp_nbits(posy) - FXP_frac_bits;
                        #ifdef VERBOSE
                        printf("\ncase 3: -ping -y, product would have %d whole bits\n", pnwbits);
                        printf("Whole bit counts:  factor: %d,  y: %d\n", \
                                        superposlog.nwbits, fxp_nbits(posy) - FXP_frac_bits);
                        #endif
                        return (pnwbits >= FXP_whole_bits)? \
                                FXP_POS_INF: \
                                fxp_pow2_pos_arg_xfactor_l(posy, \
                                                    superposlog, \
                                                    FXP_POWXY_POW_LOOPS);
                }
        }
}
