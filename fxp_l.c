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

const unsigned long FXP_PI_I64 = 0x6487ED5110B4611Aul;
const unsigned long FXP_LG2_10_I64 = 0x6A4D3C25E68DC57Ful;
const unsigned long FXP_LG2_E_I64 = 0x5C551D94AE0BF85Eul;
const unsigned long FXP_LN_2_I64 = 0x58B90BFBE8E7BCD6ul;
const unsigned long FXP_LG10_2_I64 = 0x268826A13EF3FDE6ul;

const unsigned long ULONG_ALL_ONES = ~0ul;
const unsigned long ULONG_ALL_ONES_RS1 = ~0ul >> 1;
const unsigned long ULONG_SIGN = ~ULONG_ALL_ONES_RS1;

// For the BKM lg2 calculation when using longs.
// Values in this array are: a[k] = lg2(1 + 1/2^k) represented as
// unsigned long fxp's (8 bytes,) with 63 frac bits, and
// one magnitude whole bit. Notice, not a sign but a magnitude single
// whole bit, needed to be able to represent the first value in
// the table == 1.0. That value might be removed and handled
// implicitely later, for 1 extra bit of precision.
// 63 frac bits corresponds to a precision of 18.96 decimal digits
static const unsigned long FXP_BKM_LOGS_L[] = {
        0x8000000000000000ul, 0x4AE00D1CFDEB43CFul, 0x2934F0979A3715FCul, 0x15C01A39FBD6879Ful,
        0xB31FB7D64898B3Eul, 0x5AEB4DD63BF61CCul, 0x2DCF2D0B85A4531ul, 0x16FE50B6EF08517ul,
        0xB84E236BD563BAul, 0x5C3E0FFC29D593ul, 0x2E24CA6E87E8A8ul, 0x1713D62F7957C3ul,
        0xB8A476150DFE4ul, 0x5C53AC47E94D8ul, 0x2E2A32762FA6Bul, 0x1715305002E4Aul,
        0xB8A9DED47C11ul, 0x5C55067F6E58ul, 0x2E2A89050622ul, 0x171545F3D72Bul,
        0xB8AA35640A7ul, 0x5C551C23599ul, 0x2E2A8E6E01Eul, 0x1715474E163ul,
        0xB8AA3ACD06ul, 0x5C551D7D98ul, 0x2E2A8EC491ul, 0x17154763BAul,
        0xB8AA3B239ul, 0x5C551D933ul, 0x2E2A8EC9Ful, 0x171547651ul,     //  32 entries
        0xB8AA3B28ul, 0x5C551D94ul, 0x2E2A8ECAul, 0x17154765ul,
        0xB8AA3B2ul, 0x5C551D9ul, 0x2E2A8ECul, 0x1715476ul,     // 40 entries *
        /*
        0xB8AA3B, 0x5C551D, 0x2E2A8E, 0x171547,
        0xB8AA3, 0x5C551, 0x2E2A8, 0x17154,
        0xB8AA, 0x5C55, 0x2E2A, 0x1715,
        0xB8A, 0x5C5, 0x2E2, 0x171,
        0xB8, 0x5C, 0x2E, 0x17,
        0xB, 0x5, 0x2, 0x1,
        0x0
        // Starting with the row marked with the *, each entry is exactly
        // a 4-bit right-shift of the value 4 positions earlier
        */
};

// CORDIC angles as fxp's with 62 frac bits
// Angles in this array have all a tangent that is exactly an inverse power of 2:
// tan(angle[i]) = 1/(2^-i).
// So all tangents of these angles have a single bit on, and the tangent value
// for angle[i] is the same as right-shifting by 1 the tangent of angle[i-1],
// with tan(angle[0]) == 1 (angle[0] is pi/4 == 45º in radians == 0.78539...)
// For more information on CORDIC:
// https://en.wikipedia.org/wiki/CORDIC
// https://en.wikibooks.org/wiki/Trigonometry/For_Enthusiasts/The_CORDIC_Algorithm
static const unsigned long FXP_CORDIC_ANGLES_L[] = {
0x3243F6A8885A308Dul, 0x1DAC670561BB4F68ul, 0x0FADBAFC96406EB1ul, 0x07F56EA6AB0BDB71ul, // 4
0x03FEAB76E59FBD39ul, 0x01FFD55BBA97624Aul, 0x00FFFAAADDDB94D5ul, 0x007FFF5556EEEA5Cul, // 8
0x003FFFEAAAB7776Eul, 0x001FFFFD5555BBBBul, 0x000FFFFFAAAAADDEul, 0x0007FFFFF555556Ful, // 12
0x0003FFFFFEAAAAABul, 0x0001FFFFFFD55555ul, 0x0000FFFFFFFAAAAAul, 0x00007FFFFFFF5555ul, // 16
0x00003FFFFFFFEAAAul, 0x00001FFFFFFFFD55ul, 0x00000FFFFFFFFFAAul, 0x000007FFFFFFFFF5ul, // 20
0x000003FFFFFFFFFEul, 0x0000020000000000ul, 0x0000010000000000ul, 0x0000008000000000ul, // 24
0x0000004000000000ul, 0x0000002000000000ul, 0x0000001000000000ul, 0x0000000800000000ul, // 28
0x0000000400000000ul, 0x0000000200000000ul, 0x0000000100000000ul, 0x0000000080000000ul  // 32
/*
0x0000000040000000ul, 0x0000000020000000ul, 0x0000000010000000ul, 0x0000000008000000ul,
0x0000000004000000ul, 0x0000000002000000ul, 0x0000000001000000ul, 0x0000000000800000ul,
0x0000000000400000ul, 0x0000000000200000ul, 0x0000000000100000ul, 0x0000000000080000ul,
0x0000000000040000ul, 0x0000000000020000ul, 0x0000000000010000ul, 0x0000000000008000ul,
0x0000000000004000ul, 0x0000000000002000ul, 0x0000000000001000ul, 0x0000000000000800ul,
0x0000000000000400ul, 0x0000000000000200ul, 0x0000000000000100ul, 0x0000000000000080ul,
0x0000000000000040ul, 0x0000000000000020ul, 0x0000000000000010ul, 0x0000000000000008ul,
0x0000000000000004ul, 0x0000000000000002ul, 0x0000000000000001ul, 0x0000000000000000ul,
*/
};

// Cordic scaling factor: 0.607252935008881256169446752504929
// Cordic scaling factor as fxp with 62 frac bits
static const unsigned long FXP_CORDIC_KFACTOR_L = 0x26DD3B6A10D7969A;

// Auxiliary struct used internally for lg2
typedef struct lg2tuple_l {
        int characteristic;
        unsigned long mantissa;
} lg2tuple_l;

// A super_fxp is not just a "bigger FXP," it also has its own
// fxp configuration (it's own number of whole vs. frac bits,)
// independent of the global fxp setting, always carrying this
// individual configuration within.
const super_fxp_l SFXP_ZERO_L = {0, 1, 0ul};
const super_fxp_l SFXP_LG2_10_FACTOR_L = {0, 3, FXP_LG2_10_I64};
const super_fxp_l SFXP_LG2_E_FACTOR_L = {0, 2, FXP_LG2_E_I64};
const super_fxp_l SFXP_LN_2_FACTOR_L = {0, 1, FXP_LN_2_I64};
const super_fxp_l SFXP_LG10_2_FACTOR_L = {0, 1, FXP_LG10_2_I64};

static inline int fxp_nbits_l(unsigned long x)
{
        if (x == 0) return 0;
        return FXP_LONG_BITS - __builtin_clzl(x);
};

// Functions for super_fxp_l structs

static inline super_fxp_l sfxp_l_create(unsigned int sign, \
                                        unsigned int nwbits, \
                                        unsigned long num)
{
        nwbits = (nwbits > SFXP_MAX_WBITS)? SFXP_MAX_WBITS: nwbits;
        super_fxp_l x = { (sign != 0), nwbits, num };
        return x;
}

inline super_fxp_l sfxp_l_from_fxp(int fxp1)
{
        super_fxp_l sup;
        if (fxp1 >= 0) {
                int w = fxp_get_whole_part(fxp1);
                int f = fxp_get_frac_part_bin(fxp1);
                sup.sign = 0;
                sup.nwbits = fxp_nbits(w) + 1;
                sup.number = (((unsigned long) w) << (FXP_LONG_BITS - sup.nwbits)) \
                                | (((unsigned long) f) << (FXP_long_mfrac_bits - sup.nwbits));
        } else {
                int posf = -fxp1;
                int w = fxp_get_whole_part(posf);
                int f = fxp_get_frac_part_bin(posf);
                int nb = fxp_nbits(w) + 1;
                long pnum = (((long) w) << (FXP_LONG_BITS - nb)) \
                                | (((long) f) << (FXP_long_mfrac_bits - nb));
                sup.sign = 1;
                sup.nwbits = nb;
                sup.number = ((unsigned long) pnum);
        }
        return sup;
}

inline int sfxp_l_2_fxp(super_fxp_l x)
{
        int sfb = FXP_LONG_BITS - x.nwbits;
        int fbdiff = sfb - FXP_frac_bits;
        if (x.sign) {
                if (fbdiff >= 0) {
                        // As many or more frac bits than in fxp's
                        int rbit = (x.number >> (fbdiff - 1)) & 1u;
                        unsigned long num = (x.number >> fbdiff) + rbit;
                        unsigned long whole = num >> FXP_frac_bits;
                        if (whole > FXP_whole_max) return FXP_NEG_INF;
                        return -((int) num);
                } else {
                        // Fewer frac bits than in fxp's
                        unsigned long whole = x.number >> sfb;
                        if (whole > FXP_whole_max) return FXP_NEG_INF;
                        int frac = (int) (x.number & (ULONG_ALL_ONES >> x.nwbits));
                        frac <<= -fbdiff;
                        return -fxp_bin((int) whole, frac);
                }
        } else {
                if (fbdiff >= 0) {
                        // As many or more frac bits than in fxp's
                        int rbit = (x.number >> (fbdiff - 1)) & 1u;
                        unsigned long num = (x.number >> fbdiff) + rbit;
                        unsigned long whole = num >> FXP_frac_bits;
                        if (whole > FXP_whole_max) return FXP_POS_INF;
                        return (int) num;
                } else {
                        // Fewer frac bits than in fxp's
                        unsigned long whole = x.number >> sfb;
                        if (whole > FXP_whole_max) return FXP_POS_INF;
                        int frac = (int) (x.number & (ULONG_ALL_ONES >> x.nwbits));
                        frac <<= -fbdiff;
                        return fxp_bin((int) whole, frac);
                }
        }
}

void print_sfxp_l(char * msg, super_fxp_l x)
{
        printf("\n%s", msg);
        printf("Sfxp_l : {%d, %d, x%lX}\n", x.sign, x.nwbits, x.number);
        printf(" Fxp eq: "); print_fxp(sfxp_l_2_fxp(x)); printf("\n");
}

static inline unsigned long sfxp_l_get_poswhole(super_fxp_l x)
{
        return (x.nwbits)? x.number >> (FXP_LONG_BITS - x.nwbits): 0ul;
}

static inline unsigned long get_sfxp_frac_for_bkme_l(super_fxp_l x)
{
        return (x.number << x.nwbits) >> 1;
}

static inline super_fxp_l posfxp_x_possfxp_l(int x, super_fxp_l c)
{
        int clz_m1 = __builtin_clz(x) - 1;
        unsigned int shifted_x = x << clz_m1;
        unsigned long product = dmul_ulong_x_uint(c.number, x);
        int nwb = c.nwbits + FXP_whole_bits;
        super_fxp_l super_product = sfxp_l_create(0, nwb, product);
        return super_product;
}

static inline super_fxp_l get_fxp_x_sfxp_l(int x, super_fxp_l c)
{
        super_fxp_l product;
        if (x >= 0) {
                product = posfxp_x_possfxp_l(x, c);
                product.sign = c.sign;
        } else {
                product = posfxp_x_possfxp_l(-x, c);
                product.sign = !c.sign;
        }
        return product;
}

/*
 * r-shift a long by shift, rounding its last bit
 */
static inline int rshift_long_rounding(long n, unsigned int shift)
{
        if (n >= 0) {
                int rbit = shift? (int) ((((unsigned long) n) >> (shift - 1)) & 1ul): 0;
                int res = (int) ((n >> shift) + rbit);
                //#ifdef VERBOSE
                //printf("n: %016lX (%lu),  shift: %d,  rbit: %d,  result: %08X (%d)\n", \
                //        n, n, shift, rbit, res, res);
                //#endif
                return res;
        } else {
                unsigned long posn = -n;
                unsigned long rbit = shift? (int) ((posn >> (shift - 1)) & 1ul): 0;
                return -((int) ((posn >> shift) + rbit));
        }
}

static inline unsigned int rshift_ulong_into_uint_rounding( \
                                    unsigned long n)
{
        unsigned long rbit = (n >> FXP_INT_BITS_M1) & 1ul;
        return (unsigned int) ((n >> FXP_INT_BITS) + rbit);
}

static inline int rshift_long_into_int_rounding(long n)
{
        if (n >= 0) {
                unsigned long posn = n;
                unsigned long rbit = (posn >> FXP_INT_BITS_M1) & 1ul;
                return (int) ((posn >> FXP_INT_BITS) + rbit);
        } else {
                unsigned long posn = -n;
                unsigned long rbit = (posn >> FXP_INT_BITS_M1) & 1ul;
                return -((int) ((posn >> FXP_INT_BITS) + rbit));
        }
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

static inline unsigned long fxp_bkm_lmode_l(unsigned long argument, \
                                            const int MAX_LOOPS)
{
        #ifdef VERBOSE_BKM
        printf("bkm_lmode_l argument: x%lX\n", argument);
        #endif
        unsigned long x = FXP_BKM_X_ONE_L;
        unsigned long z, xs, y = 0ul;
        // The mantissa value will remain in the range of lg2(x),
        // x in [1, 2), meaning mantissa always in [0, 1)
        // Watch out we are using two different shifted
        // configurations simultaneously, but independently:
        // - "A" alignment (A for Array and argument)
        //   BKM_LOGS Array values, argument, and the
        //   mantissa all have 1 unsigned whole and 63 frac bits
        //   for more accuracy
        // - "X" alignment:
        //   x aligned with FXP_BKM_X_ONE_L: 2 whole and 62 frac bits
        //   (at least 2 whole bits for the BKM log algorithm,
        //   given that z can occasionally exceed 2)
        for (int shift = 1; shift <= MAX_LOOPS; shift++) {
                #ifdef VERBOSE_BKM
                printf("bkm_lmode_l shift:%2d - ", shift);
                #endif
                xs = (x >> shift);
                z = x + xs;
                if (z <= argument) {
                        x = z;
                        y += FXP_BKM_LOGS_L[shift];
                        #ifdef VERBOSE_BKM
                        printf("Updated x:x%16lX  Updated y:x%16lX", x, y);
                        #endif
                }
                #ifdef VERBOSE_BKM
                printf("\n");
                #endif
        }
        #ifdef VERBOSE
        printf("bkm_lmode_l final mantissa is: x%lX\n", y);
        #endif
        return y;
}

/*
 * Internal auxiliary function that calculates the characteristic
 * and rounded mantissa of lg2, returning them separately in a
 * struct (tuple):
 * The characteristic as int with 0 frac bits,
 * the mantissa as unsigned long frac bits (exact same
 * alignment as the BKM array values: 63 frac bits)
 * Uses the BKM L-Mode algorithm (requires table of pre-calculated
 * values) and longs.
 *
 * For more details on BKM:
 * https://en.wikipedia.org/wiki/BKM_algorithm
 */
static inline lg2tuple_l fxp_get_lg2tuple_l(int fxp1, \
                                    const int MAX_LOOPS)
{
        lg2tuple_l result;
        // Assumes fxp1 is for sure > 0
        int clz = __builtin_clz((unsigned int) fxp1);
        // clz > 0 since at least sign bit == 0
        int nbx = FXP_INT_BITS - clz;
        // Assign the characteristic
        result.characteristic = \
                ((fxp1 <= FXP_almost1) || (fxp1 >= FXP_two))? \
                        (nbx - FXP_frac_bits - 1): 0;
        #ifdef VERBOSE
        printf("\nget_lg2tuple_l: fxp1: x%X  clz: %d,  nbx: %d\n", fxp1, clz, nbx);
        printf("get_lg2tuple_l: characteristic is : %d\n", result.characteristic);
        #endif

        // Prepare the argument for bkm_lmode
        // Here we have already calculated the log characteristic c
        // Now lshifting the original number to have it as unsigned
        // long fxp with 2 whole bits (2 needed because z in the
        // bkm_lmode algorithm can intermitently get a value > 2)
        unsigned long argument = ((unsigned long) fxp1) << \
                            (clz + FXP_INT_BITS_M1);
        // Get the mantissa with bkm_lmode
        result.mantissa = fxp_bkm_lmode_l(argument, MAX_LOOPS);
        // Return result as tuple {characteristic, mantissa}
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
        lg2tuple_l tup = fxp_get_lg2tuple_l(fxp1, FXP_lg2_maxloops);
        int rbit = (tup.mantissa >> FXP_int_plus_whole_bits_m2) & 1u;
        int m = (int) (tup.mantissa >> FXP_int_plus_whole_bits_m1) + rbit;
        #ifdef VERBOSE
        printf("\nlg2_l as tuple result: {%d, %lu} == {x%X,x%lX}\n", \
                        tup.characteristic, tup.mantissa, \
                        tup.characteristic, tup.mantissa);
        printf("Final int mantissa: x%X  (rounding bit was %d)\n", m, (int) rbit);
        #endif
        // Return the complete logarithm (charact + mantissa) as fxp
        if (tup.characteristic < FXP_whole_min) {
                if ((tup.characteristic == FXP_whole_min_m1) \
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
        return (tup.characteristic >= 0)?
                        ((tup.characteristic << FXP_frac_bits) + m):
                        (-(-tup.characteristic << FXP_frac_bits) + m);
}

static inline super_fxp_l sfxp_l_from_cm(int c, \
                                         unsigned long m)
{
        #ifdef VERBOSE
        printf("c: %d, m: %lX\n", c, m);
        printf("characteristic:  "); print_int_as_bin(c, 0); printf("\n");
        printf("mantissa      :  "); print_ulong_as_bin(m); printf("\n");
        #endif
        if (c >= 0) {
                int nb = fxp_nbits(c);
                unsigned long num = (((unsigned long) c) << (FXP_LONG_BITS_M1 - nb));
                int rbit = (m >> (nb - 1)) & 1;
                m = (m >> nb) + rbit;
                num = num + m;
                super_fxp_l x = { 0, nb + 1, num };
                #ifdef VERBOSE
                printf("sfxp num      :  "); print_ulong_as_bin(num); printf("\n");
                #endif
                return x;
        } else {
                int pc = -c;
                int nb = fxp_nbits(pc);
                unsigned long num = (((unsigned long) pc) << (FXP_LONG_BITS_M1 - nb));
                int rbit = (m >> (nb - 1)) & 1;
                m = (m >> nb) + rbit;
                num = num - m;
                super_fxp_l x = { 1, nb + 1, num };
                #ifdef VERBOSE
                printf("sfxp num      : -"); print_ulong_as_bin(num); printf("\n");
                #endif
                return x;
        }
}

/*
 * Internal auxiliary function that calculates the characteristic
 * and rounded mantissa of lg2, returning them already summed up
 * into a super_fxp_l
 * Uses the BKM L-Mode algorithm (requires table of pre-calculated
 * values) and longs.
 *
 * For more details on BKM:
 * https://en.wikipedia.org/wiki/BKM_algorithm
 */
static inline super_fxp_l fxp_get_lg2_as_sfxp_l(int fxp1, \
                                        const int MAX_LOOPS)
{
        int clz = __builtin_clz((unsigned int) fxp1);
        int nbx = FXP_INT_BITS - clz;
        // Assign the characteristic
        int c = ((fxp1 <= FXP_almost1) || (fxp1 >= FXP_two))? \
                        (nbx - FXP_frac_bits_p1): 0;
        #ifdef VERBOSE
        printf("\nlg2_as_sfxp_l: fxp1: x%X  clz: %d,  nbx: %d\n", fxp1, clz, nbx);
        printf("lg2_as_sfxp_l: characteristic is : %d\n", c);
        #endif
        // Prepare argument for bkm_lmode
        unsigned long argument = ((unsigned long) fxp1) << \
                            (clz + FXP_INT_BITS_M1);
        // Get the mantissa with bkm_lmode
        unsigned long m = fxp_bkm_lmode_l(argument, MAX_LOOPS);
        // Return result as a super_fxp_l
        return sfxp_l_from_cm(c, m);
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

static inline super_fxp_l get_sfxp_x_sfxp_l(super_fxp_l x, super_fxp_l y)
{
        // Sign of the product assuming arguments are never both negative
        super_fxp_l prod = { x.sign | y.sign, \
                             x.nwbits + y.nwbits,
                             dmul_ulongs(x.number, y.number) };
        return prod;
}

/*
 * Calculates lg2 and then multiplies by the given factor
 * using longs.
 */
static inline int lg2_x_factor_l(int fxp1, super_fxp_l factor)
{
        super_fxp_l lg2x = fxp_get_lg2_as_sfxp_l(fxp1, FXP_LOGX_LOOPS);
        super_fxp_l prod = get_sfxp_x_sfxp_l(factor, lg2x);
        return sfxp_l_2_fxp(prod);
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
                printf("k:%2d,  z:%16lX,  z as ld:%.19Lf\n", \
                        k, z, print_ulong_as_ld(z, 63));
                #endif
                if (z <= argument) {
                        y = z;
                        x = x + (x >> k);
                        #ifdef VERBOSE_BKM
                        printf("\tUpdating y (x%16lX) und x (x%16lX)\n", y, x);
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
 * The frac ulong argument must already have the same alignment
 * as the BKM array values (1 whole bit).
 * Input value for 2^x here is x = whole + frac
 */
static inline int fxp_pow2_wpos_l(int whole, \
                                  unsigned long frac, \
                                  const int MAX_LOOPS)
{
        #ifdef VERBOSE
        printf("pow2_wpos_l for x%X.%lX\n", whole, frac);
        #endif
        if (whole >= FXP_whole_bits_m1) return FXP_POS_INF;
        if ((whole == 0) && (frac == 0ul)) return FXP_one;
        // Argument frac will be in [0, 1)
        unsigned long x = fxp_bkm_emode_l(frac, MAX_LOOPS);
        // When one of the operands is a shifted 1, we don't really
        // need to call the (expensive) dmul operation, since the
        // result will be identical to the other operand anyway,
        // just shifted appropriately
        int shift = FXP_int_plus_whole_bits_m2 - whole;
        int shifted = (int) rshift_ulong_rounding(x, shift);
        #ifdef VERBOSE
        printf("x from bkm_e: x%lX\n", x);
        print_ulong_as_bin(x); printf("\n");
        printf("shift is %d\n", shift);
        printf("pow2_wpos_l result is ");
        print_uint_as_bin(shifted); printf("\n");
        #endif
        return (shifted != FXP_UNDEF)? shifted: FXP_POS_INF;
}

/*
 * pow2 calculation (2^fxp1) of a negative argument,
 * using the BKM (E-mode) algorithm and longs
 * Argument whole and frac are positive, both corresponding
 * to minus the corresponding parts of the original
 * negative fxp (see how the arguments get built
 * in fxp_pow2_l.) Frac component must already have the same
 * alignment as the BKM array values (1 whole bit).
 */
static inline int fxp_pow2_wneg_l(int whole, \
                                  unsigned long frac, \
                                  const int MAX_LOOPS)
{
        #ifdef VERBOSE
        printf("pow2_wneg_l for -x%X.%lX\n", whole, frac);
        #endif
        if (whole > FXP_frac_bits) return 0;
        if ((whole == 0) && (frac == 0ul)) return FXP_one;
        // Notice argument a for bkm_emode will be in (0, 1]
        unsigned long a = FXP_BKM_A_ONE_L - frac;
        unsigned long x = fxp_bkm_emode_l(a, MAX_LOOPS);
        int shift = FXP_int_plus_whole_bits_m1 + whole;
        int shifted = (int) rshift_ulong_rounding(x, shift);
        #ifdef VERBOSE
        printf("x from bkm_e: x%lX\n", x);
        printf("shift is %d\n", shift);
        print_ulong_as_bin(x); printf("\n");
        print_uint_as_bin(shifted); printf("\n");
        #endif
        return (shifted != FXP_UNDEF)? shifted: FXP_POS_INF;
}

static inline unsigned long get_fxp_frac_for_bkme_l(int frac)
{
        #ifdef VERBOSE
        printf("get_fxp_frac_for_bkme_l: Fracbits is %d\n", FXP_frac_bits);
        printf("Frac is x%X (d: %d)\n", frac, frac);
        printf("Shift for frac is: %d\n", FXP_int_plus_whole_bits_m1);
        #endif
        return ((unsigned long) frac) << FXP_int_plus_whole_bits_m1;
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
                int w = fxp_get_whole_part(fxp1);
                unsigned long bkmearg = get_fxp_frac_for_bkme_l( \
                                                fxp_get_frac_part_bin(fxp1));
                return fxp_pow2_wpos_l( w, bkmearg, FXP_POWX_LOOPS );
        } else {
                int pfxp1 = -fxp1;
                int w = fxp_get_whole_part(pfxp1);
                unsigned long bkmearg = get_fxp_frac_for_bkme_l( \
                                                fxp_get_frac_part_bin(pfxp1));
                return fxp_pow2_wneg_l( w, bkmearg, FXP_POWX_LOOPS );
        }
}

// Calculate the pow2 of a NON-negative argument x
// multiplied by a factor c
static inline int fxp_pow2_pos_arg_xfactor_l(int x, \
                                        super_fxp_l factorc, \
                                        const int MAX_LOOPS)
{
        super_fxp_l xc = get_fxp_x_sfxp_l(x, factorc);
        #ifdef VERBOSE
        print_sfxp_l("Pos xc (sfxp_l) is: ", xc);
        #endif
        int w = (int) sfxp_l_get_poswhole(xc);
        unsigned long bkmearg = get_sfxp_frac_for_bkme_l(xc);
        return fxp_pow2_wpos_l( w, bkmearg, MAX_LOOPS );
}

// Calculate the pow2 of a negative argument x
// multiplied by a factor C
static inline int fxp_pow2_neg_arg_xfactor_l(int x, \
                                        super_fxp_l factorc, \
                                        const int MAX_LOOPS)
{
        super_fxp_l xc = get_fxp_x_sfxp_l(x, factorc);
        #ifdef VERBOSE
        print_sfxp_l("Neg xc (sfxp_l) is: ", xc);
        #endif
        int w = (int) sfxp_l_get_poswhole(xc);
        unsigned long bkmearg = get_sfxp_frac_for_bkme_l(xc);
        return fxp_pow2_wneg_l( w, bkmearg, MAX_LOOPS );
}

// Implementation of exp_l(x) using pow2_l():
// e^x == 2^( x * lg2(e) )
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

// Implementation of pow10_l(x) using pow2_l() and lg2_l:
// 10^x == 2^( x * lg2(10) )
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
        super_fxp_l slg = fxp_get_lg2_as_sfxp_l(fxp1, FXP_SQRT_LG_LOOPS);
        // Halving that value
        slg.nwbits--;
        // parameters for pow2
        int w = (int) sfxp_l_get_poswhole(slg);
        unsigned long bkmearg = get_sfxp_frac_for_bkme_l(slg);
        return (slg.sign)?
                fxp_pow2_wneg_l(w, bkmearg, FXP_sqrt_pw_loops):
                fxp_pow2_wpos_l(w, bkmearg, FXP_sqrt_pw_loops);
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
        super_fxp_l slg2x = fxp_get_lg2_as_sfxp_l(x, FXP_POWXY_LG_LOOPS);
        #ifdef VERBOSE
        printf("powxy_l:\n");
        printf("      x: "); print_fxp(x); printf("\n");
        printf("      y: "); print_fxp(y);
        print_sfxp_l(" lg2(x): ", slg2x);
        #endif
        // Return appropriately signed 2^( y * lg2x )
        if (slg2x.sign) {
                slg2x.sign = 0; // <- negating in place
                if (y >= 0) {
                        #ifdef VERBOSE
                        int pnwbits = slg2x.nwbits + fxp_nbits(y) - FXP_frac_bits;
                        printf("\ncase 2: -lgx +y, product would have %d whole bits\n", pnwbits);
                        printf("Whole bit counts:  factor: %d,  y: %d\n", \
                                        slg2x.nwbits, fxp_nbits(y) - FXP_frac_bits);
                        #endif
                        return fxp_pow2_neg_arg_xfactor_l(y, slg2x, FXP_POWXY_PW_LOOPS);
                } else {
                        int posy = -y;
                        int pnwbits = slg2x.nwbits + fxp_nbits(posy) - FXP_frac_bits_p1;
                        #ifdef VERBOSE
                        printf("\ncase 3: -lgx -y, product would have %d whole bits\n", pnwbits);
                        printf("Whole bit counts:  factor: %d,  y: %d\n", \
                                        slg2x.nwbits, fxp_nbits(posy) - FXP_frac_bits);
                        #endif
                        return (pnwbits > FXP_whole_bits)? \
                                        FXP_POS_INF: \
                                        fxp_pow2_pos_arg_xfactor_l(posy, slg2x, \
                                                                FXP_POWXY_PW_LOOPS);
                }
        } else {
                if (y >= 0) {
                        int pnwbits = slg2x.nwbits + fxp_nbits(y) - FXP_frac_bits;
                        #ifdef VERBOSE
                        printf("\ncase 0: +lgx +y, product would have %d whole bits\n", pnwbits);
                        printf("Whole bit counts:  factor: %d,  y: %d\n", \
                                        slg2x.nwbits, \
                                        fxp_nbits(y) - FXP_frac_bits);
                        #endif
                        return (pnwbits > FXP_whole_bits)?
                                        FXP_POS_INF: \
                                        fxp_pow2_pos_arg_xfactor_l(y, slg2x, \
                                                                FXP_POWXY_PW_LOOPS);
                } else {
                        int posy = -y;
                        int pnwbits = slg2x.nwbits + fxp_nbits(posy) - FXP_frac_bits;
                        #ifdef VERBOSE
                        printf("\ncase 1: +lgx -y, product would have %d whole bits\n", pnwbits);
                        printf("Whole bit counts:  factor: %d,  y: %d\n", \
                                        slg2x.nwbits, \
                                        fxp_nbits(posy) - FXP_frac_bits);
                        #endif
                        return (pnwbits > FXP_whole_bits)?
                                        0: \
                                        fxp_pow2_neg_arg_xfactor_l(posy, slg2x, \
                                                                FXP_POWXY_PW_LOOPS);
                }
        }
}


/*
 *=======================================
 * Trigonometric functions
 *=======================================
 * Note: these initial implementations requires +/-PI to be
 * representable in the current fxp configuration, so at least
 * 3 whole bits are needed.
 * TODO: consider using sfxp's to work around this, enabling
 * trig functions for all frac bit configs up.
 */

/*
 * Simultaneous calculation of cos(x) and sin(x)
 * using the CORDIC algorithm.
 *
 * For more details on CORDIC:
 * https://en.wikipedia.org/wiki/CORDIC
 * https://en.wikipedia.org/wiki/BKM_algorithm
 */
static inline fxptuple_l fxp_cordic_cossin_l(long x)
{
        // argument assumed to be already in range [-pi/2, pi/2],
        // and l-shifted like the Cordic factor and angles,
        // so with 62 frac bits
        long lc = FXP_CORDIC_KFACTOR_L;
        long ls = 0l;
        long newc;
        long a = 0;
        #ifdef VERBOSE
        long double adenom = (1ul << (FXP_LONG_BITS - 2));
        long double cdenom = (1ul << (FXP_LONG_BITS - 2));
        printf("Running CORDIC for %.4Lf (%.4Lfº)", x/adenom, x*180.0L/(adenom*acosl(-1.0L)));
        #endif
        for (int i = 0; i < FXP_cordic_loops; i++) {
                if (a < x) {
                        a += FXP_CORDIC_ANGLES_L[i];
                        newc = lc - (ls >> i);
                        ls += (lc >> i);
                        lc = newc;
                        #ifdef VERBOSE
                        printf("\n%02d: + rot: ", i);
                        #endif
                } else {
                        a -= FXP_CORDIC_ANGLES_L[i];
                        newc = lc + (ls >> i);
                        ls -= (lc >> i);
                        lc = newc;
                        #ifdef VERBOSE
                        printf("\n%02d: - rot: ", i);
                        #endif
                }
                #ifdef VERBOSE
                printf("x%016lX (%7.4Lf),  c:x%016lX (%7.4Lf),  s:x%016lX (%7.4Lf)", \
                            a, a/adenom, lc, lc/cdenom, ls, ls/cdenom);
                #endif
        }
        fxptuple_l tup_l = { lc, ls };
        return tup_l;
}

static inline fxptuple_l fxp_cossin_arg_inrange_l(int fxp1)
{
        // argument must be already in range [-pi/2, pi/2]
        #ifdef VERBOSE
        printf("\nAngle 0x%08X (%d, %.4Lfº) in [-pi/2, pi/2]\n", \
                fxp1, fxp1, fxp2ld(fxp1)*180.0L/acosl(-1.0L));
        #endif
        if (fxp1 >= 0)
                return fxp_cordic_cossin_l( \
                            ((long) fxp1) << FXP_int_plus_whole_bits_m2);
        else {
                return fxp_cordic_cossin_l( \
                            -(((long) -fxp1) << FXP_int_plus_whole_bits_m2) );
        }
}

static inline fxptuple_l fxp_cossin_arg_above_phalfpi_l(int fxp1)
{
        #ifdef VERBOSE
        printf("\nAngle 0x%08X (%d, %.4Lfº) > pi/2\n", \
                fxp1, fxp1, fxp2ld(fxp1)*180.0L/acosl(-1.0L));
        #endif
        // Get an equivalent angle within [-pi/2, pi/2]
        unsigned int times = 0;
        const long angle0 = ((unsigned long) fxp1) << FXP_INT_BITS;
        long angle = angle0;
        do {
                times++;
                angle -= FXP_shifted_pi_l;
        } while (angle > FXP_shifted_phalfpi_l);
        #ifdef VERBOSE
        int iangle = rshift_long_rounding(angle, FXP_INT_BITS);
        printf("Original num: 0x%016lX\n", angle0);
        printf("shifted pi_l: 0x%016lX\n", FXP_shifted_pi_l);
        printf("Final angle : 0x%016lX  ", angle);
        printf("~(%d, %.4Lfº) in [-pi/2, pi/2]\n", \
                        iangle, fxp2ld(iangle)*180.0L/acosl(-1.0L));
        #endif
        // l-shift angle as required by the cordic function
        angle = (angle >= 0)? angle << FXP_whole_bits_m2
                            : -((-angle) << FXP_whole_bits_m2);
        fxptuple_l tup_l = fxp_cordic_cossin_l(angle);
        if (times & 1u) {
                // Adjust signs of results
                tup_l.a = -tup_l.a;
                tup_l.b = -tup_l.b;
        }
        return tup_l;
}

static inline fxptuple_l fxp_cossin_arg_below_nhalfpi_l(int fxp1)
{
        #ifdef VERBOSE
        printf("\nAngle 0x%08X (%d, %Lfº) < -pi/2\n", \
                fxp1, fxp1, fxp2ld(fxp1)*180.0L/acosl(-1.0L));
        #endif
        // Get an equivalent angle within [-pi/2, pi/2]
        unsigned int times = 0;
        const long angle0 = -((long) (((unsigned long) -fxp1) << FXP_INT_BITS));
        long angle = angle0;
        do {
                times++;
                angle += FXP_shifted_pi_l;
        } while (angle < FXP_shifted_nhalfpi_l);
        #ifdef VERBOSE
        printf("Original num:\t0x%016lX\n", angle0);
        printf("shifted pi_l:\t0x%016lX\n", FXP_shifted_pi_l);
        printf("Final angle :\t0x%016lX  ", angle);
        int iangle = rshift_long_rounding(angle, FXP_INT_BITS);
        printf("~(%d, %.4Lfº) in [-pi/2, pi/2]\n", \
                        iangle, fxp2ld(iangle)*180.0L/acosl(-1.0L));
        #endif
        // l-shift angle as required by the cordic function
        angle = (angle >= 0)? angle << FXP_whole_bits_m2
                            : -((-angle) << FXP_whole_bits_m2);
        fxptuple_l tup_l = fxp_cordic_cossin_l(angle);
        if (times & 1u) {
                // Adjust signs of results
                tup_l.a = -tup_l.a;
                tup_l.b = -tup_l.b;
        }
        return tup_l;
}

// Auxiliary function to process depending on where the
// argument lands with respect to range [-pi/2, pi/2]
static inline fxptuple_l fxp_cos_and_sin_l_aux(int fxp1)
{
        if (fxp1 > FXP_shifted_phalfpi)
                return fxp_cossin_arg_above_phalfpi_l(fxp1);
        if (fxp1 < FXP_shifted_nhalfpi)
                return fxp_cossin_arg_below_nhalfpi_l(fxp1);
        return fxp_cossin_arg_inrange_l(fxp1);
}

// If needing both cos(x) and sin(x), this is the function to call
// to get them simultaneously most efficiently (basically at half
// the cost of calling cos(x), and then calling sin(x))
fxptuple fxp_cos_and_sin_l(int fxp1)
{
        // Argument here can be any fxp value, even outside [-pi/2, pi/2].
        if ((fxp1 == FXP_UNDEF) || (fxp1 == FXP_NEG_INF) || (fxp1 == FXP_POS_INF))
                return FXP_TUPLE_UNDEFS;
        fxptuple_l tup_l = fxp_cos_and_sin_l_aux(fxp1);
        fxptuple tup = { rshift_long_rounding(tup_l.a, FXP_int_plus_whole_bits_m2), \
                         rshift_long_rounding(tup_l.b, FXP_int_plus_whole_bits_m2) };
        return tup;
}

// If needing only cos(x), but not interested in sin(x)
int fxp_cos_l(int fxp1)
{
        if ((fxp1 == FXP_UNDEF) || (fxp1 == FXP_NEG_INF) || (fxp1 == FXP_POS_INF))
                return FXP_UNDEF;
        fxptuple_l tup_l = fxp_cos_and_sin_l_aux(fxp1);
        return rshift_long_rounding(tup_l.a, FXP_int_plus_whole_bits_m2);
}

// If needing only sin(x), but not interested in cos(x)
int fxp_sin_l(int fxp1)
{
        if ((fxp1 == FXP_UNDEF) || (fxp1 == FXP_NEG_INF) || (fxp1 == FXP_POS_INF))
                return FXP_UNDEF;
        fxptuple_l tup_l = fxp_cos_and_sin_l_aux(fxp1);
        return rshift_long_rounding(tup_l.b, FXP_int_plus_whole_bits_m2);
}

int fxp_tan_l(int fxp1)
{
        return 0;
}

int fxp_asin_l(int fxp1)
{
        return 0;
}

int fxp_acos_l(int fxp1)
{
        return 0;
}

int fxp_atan_l(int fxp1){
        return 0;
}
