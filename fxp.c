/* SPDX-License-Identifier: MIT */
/*
 * fxp.c
 *
 * An implementation of Binary Fixed Point numbers
 * encoding them into integers, together with their
 * arithmetic operations +, -, *, and /, and more.
 *
 * By Raul Saavedra, Bonn, Germany
 */

#include "fxp.h"
#include <stdio.h>
#include <stdlib.h>
#include "ulongy.h"

//#define VERBOSE 1
#ifdef VERBOSE
#include "print_as_bits.h"
#endif

const int FXP_INT_BITS = ((int) sizeof(int)) * 8;
const int FXP_INT_BITS_M1 = FXP_INT_BITS - 1;
const int FXP_INT_BITS_M1_NEG = -FXP_INT_BITS_M1;
const int FXP_INT_BITS_M2 = FXP_INT_BITS - 2;
const int FXP_POS_INF = INT_MAX;
const int FXP_NEG_INF = -INT_MAX;
const int FXP_MAX = FXP_POS_INF - 1;
const long FXP_MAX_L = (long) FXP_MAX;
const int FXP_MIN = FXP_NEG_INF + 1;
// Saving the true most negative int for "undefined"
const int FXP_UNDEF = INT_MIN;

// Default # frac bits: half the size of ints
const int FXP_FRAC_BITS_DEF = FXP_INT_BITS / 2;
// Table of some frac bits and equivalent decimal digit precision:
// Fbits  Dec.Precision = FBits*ln(2)/ln(10)
//	60  	18.0618
//	50  	15.0515
//	30  	 9.0309
//	20  	 6.0206
//  16       4.8165
//	10  	 3.0103
//	7   	 2.1072
//	4   	 1.2041
// For 32-bit ints -> default is 16 frac bits,
// offering better than 4 decimal digits of precision.
// As a minimum at least 4 frac bits, offering
// at least 1 decimal digit of precision
const int FXP_FRAC_BITS_MIN = 4;

const int FXP_WORD_BITS = FXP_INT_BITS >> 1;
const int FXP_WORD_BITS_M1 = FXP_WORD_BITS - 1;
const unsigned int FXP_RWORD_MASK = ((1u << FXP_WORD_BITS) - 1);
const unsigned int FXP_LWORD_MASK = FXP_RWORD_MASK \
                                    << FXP_WORD_BITS;
const unsigned long FXP_RINT_MASK = ((1lu << FXP_INT_BITS) - 1);
const unsigned long FXP_LINT_MASK = FXP_RINT_MASK \
                                    << FXP_INT_BITS;
const int FXP_LONG_BITS = ((int) sizeof(long)) * 8;
const int FXP_LONG_BITS_M1 = FXP_LONG_BITS - 1;

// Allowing for no whole part, so other than the sign bit,
// all bits used for fraction part. This use case represents
// the range [-0.999..., 0.999...], or equivalently: (-1, 1)
// (e.g. actual whole values -1 and 1 outside the valid range)
const int FXP_FRAC_BITS_MAX = FXP_INT_BITS_M1;
const int FXP_FRAC_MAX_DEC = 9999999;

// Default number of bits to use for the frac part.
int FXP_frac_bits = FXP_FRAC_BITS_DEF;
// FXP_frac_bits can be changed dynamically calling
// fxp_set_frac_bits(), but it remains unchanged until calling
// that function again.
// All such static variables which will only change
// if modifying the number of frac bits are named here
// starting with "FXP_" (uppercase,) but the rest of the
// variable name in lowercase

int FXP_frac_bits_m1 = FXP_FRAC_BITS_DEF - 1;

// Default number of bits for the whole part (includes sign bit)
int FXP_whole_bits = FXP_INT_BITS - FXP_FRAC_BITS_DEF;
int FXP_whole_bits_m1 = FXP_INT_BITS_M1 - FXP_FRAC_BITS_DEF;
int FXP_whole_bits_m2 = FXP_INT_BITS_M1 - FXP_FRAC_BITS_DEF - 1;

// FXP_FRAC_MASK should correspond to 2^FXP_FRAC_BITS - 1
unsigned int FXP_frac_mask = ((1u << FXP_FRAC_BITS_DEF) - 1);
unsigned int FXP_frac_max = ((1u << FXP_FRAC_BITS_DEF) - 1);
int FXP_frac_max_int = ((1u << FXP_FRAC_BITS_DEF) - 1);
int FXP_frac_max_neg = -((1u << FXP_FRAC_BITS_DEF) - 1);
unsigned int FXP_frac_max_p1 = (1u << FXP_FRAC_BITS_DEF);

// Default max and min valid values for the whole part of the fxp's
int FXP_whole_max = FXP_MAX >> FXP_FRAC_BITS_DEF;
int FXP_whole_min = -(FXP_MAX >> FXP_FRAC_BITS_DEF);
int FXP_whole_min_m1 = -(FXP_MAX >> FXP_FRAC_BITS_DEF) - 1;

// Default max and min valid values for the conversion types
float FXP_min_fx = \
                -((float) (FXP_MAX >> FXP_FRAC_BITS_DEF) \
                    + ((float) (FXP_MAX & ((1uL << FXP_FRAC_BITS_DEF) - 1))) \
                        / (1uL << (FXP_FRAC_BITS_DEF))) \
                - (((float) 1uL) / (1uL << FXP_FRAC_BITS_DEF)) / 2.0;
float FXP_min_f = \
                -((float) (FXP_MAX >> FXP_FRAC_BITS_DEF) \
                    + ((float) (FXP_MAX & ((1uL << FXP_FRAC_BITS_DEF) - 1))) \
                        / (1uL << (FXP_FRAC_BITS_DEF)));
float FXP_max_f = \
                ((float) (FXP_MAX >> FXP_FRAC_BITS_DEF) \
                    + ((float) (FXP_MAX & ((1uL << FXP_FRAC_BITS_DEF) - 1))) \
                        / (1uL << (FXP_FRAC_BITS_DEF)));
float FXP_max_fx = \
                (float) (FXP_MAX >> FXP_FRAC_BITS_DEF) \
                    + ((float) (FXP_MAX & ((1uL << FXP_FRAC_BITS_DEF) - 1))) \
                        / (1uL << (FXP_FRAC_BITS_DEF))
                + (((float) 1uL) / (1uL << FXP_FRAC_BITS_DEF)) / 2.0;

double FXP_min_dx = \
                -((double) (FXP_MAX >> FXP_FRAC_BITS_DEF) \
                    + ((double) (FXP_MAX & ((1uL << FXP_FRAC_BITS_DEF) - 1))) \
                        / (1uL << (FXP_FRAC_BITS_DEF))) \
                - (((double) 1uL) / (1uL << FXP_FRAC_BITS_DEF)) / 2.0;
double FXP_min_d = \
                -((double) (FXP_MAX >> FXP_FRAC_BITS_DEF) \
                    + ((double) (FXP_MAX & ((1uL << FXP_FRAC_BITS_DEF) - 1))) \
                        / (1uL << (FXP_FRAC_BITS_DEF)));
double FXP_max_d = \
                ((double) (FXP_MAX >> FXP_FRAC_BITS_DEF) \
                    + ((double) (FXP_MAX & ((1uL << FXP_FRAC_BITS_DEF) - 1))) \
                        / (1uL << (FXP_FRAC_BITS_DEF)));
double FXP_max_dx = \
                (double) (FXP_MAX >> FXP_FRAC_BITS_DEF) \
                    + ((double) (FXP_MAX & ((1uL << FXP_FRAC_BITS_DEF) - 1))) \
                        / (1uL << (FXP_FRAC_BITS_DEF)) \
                + (((double) 1uL) / (1uL << FXP_FRAC_BITS_DEF)) / 2.0;

long double FXP_min_ldx = \
                -((long double) (FXP_MAX >> FXP_FRAC_BITS_DEF) \
                    + ((long double) (FXP_MAX & ((1uL << FXP_FRAC_BITS_DEF) - 1))) \
                        / (1uL << (FXP_FRAC_BITS_DEF))) \
                - (((long double) 1uL) / (1uL << FXP_FRAC_BITS_DEF)) / 2.0L;
long double FXP_min_ld = \
                -((long double) (FXP_MAX >> FXP_FRAC_BITS_DEF) \
                    + ((long double) (FXP_MAX & ((1uL << FXP_FRAC_BITS_DEF) - 1))) \
                        / (1uL << (FXP_FRAC_BITS_DEF)));
long double FXP_max_ld = \
                ((long double) (FXP_MAX >> FXP_FRAC_BITS_DEF) \
                    + ((long double) (FXP_MAX & ((1uL << FXP_FRAC_BITS_DEF) - 1))) \
                        / (1uL << (FXP_FRAC_BITS_DEF)));
long double FXP_max_ldx = \
                ((long double) (FXP_MAX >> FXP_FRAC_BITS_DEF) \
                    + ((long double) (FXP_MAX & ((1uL << FXP_FRAC_BITS_DEF) - 1))) \
                        / (1uL << (FXP_FRAC_BITS_DEF))) \
                + (((long double) 1uL) / (1uL << FXP_FRAC_BITS_DEF)) / 2.0L;
long double FXP_htiniest_ld = \
                (((long double) 1uL) \
                    / (1uL << FXP_FRAC_BITS_DEF)) / 2.0L;

typedef struct tuple {
        int ping;
        ulongy pong;
} tuple;

// transcendental constants (2 whole, 30 and 62 frac bits)
// e = 2.718281828...
const unsigned int FXP_E_I32 = 0xADF85458u;
const unsigned int FXP_E_I32_X = 0xA2BB4A9Bu;

// pi = 3.14159265...
const unsigned int FXP_PI_I32 = 0xC90FDAA2u;
const unsigned int FXP_PI_I32_X = 0x2168C235u;

// transcendental constants (0 whole, 32 and 64 frac bits)
// ln(2) = 0.69314718...
const unsigned int FXP_LN_2_I32 = 0xB17217F7u;
const unsigned int FXP_LN_2_I32_X = 0xD1CF7C72u;
const ulongy FXP_LN_2_ULONGY = {FXP_LN_2_I32, FXP_LN_2_I32_X};

// lg10(2) = 0.30102999...
const unsigned int FXP_LG10_2_I32 = 0x4D104D42u;
const unsigned int FXP_LG10_2_I32_X = 0x7DE7FD01u;
const ulongy FXP_LG10_2_ULONGY = {FXP_LG10_2_I32, FXP_LG10_2_I32_X};

// lg2(e) = 1.44269504... (1 whole, 31 and 63 frac bits)
const unsigned int FXP_LG2_E_I32 = 0xB8AA3B29u;
const unsigned int FXP_LG2_E_I32_X = 0x5C17F19Eu;
const unsigned int FXP_LG2_E_WBITS = 1;
const ulongy FXP_LG2_E_ULONGY = {FXP_LG2_E_I32, FXP_LG2_E_I32_X};

// lg2(10) = 3.32192809... (2 whole, 30 and 62 frac bits)
const unsigned int FXP_LG2_10_I32 = 0xD49A784Bu;
const unsigned int FXP_LG2_10_I32_X = 0xCD1B8B51u;
const unsigned int FXP_LG2_10_WBITS = 2;
const ulongy FXP_LG2_10_ULONGY = {FXP_LG2_10_I32, FXP_LG2_10_I32_X};

static const int FXP_DEF_SHIFT = FXP_INT_BITS - FXP_FRAC_BITS_DEF;
unsigned int FXP_shifted_e = (FXP_E_I32 >> (FXP_DEF_SHIFT - 2));
unsigned int FXP_shifted_pi = (FXP_PI_I32 >> (FXP_DEF_SHIFT - 2));
unsigned int FXP_shifted_ln_2 = (FXP_LN_2_I32 >> FXP_DEF_SHIFT);
unsigned int FXP_shifted_lg10_2 = (FXP_LG10_2_I32 >> FXP_DEF_SHIFT);

// Auxiliary variables used in the lg2 implementations
int FXP_half = 1 << (FXP_FRAC_BITS_DEF - 1);
int FXP_two = 1 << (FXP_FRAC_BITS_DEF + 1);
int FXP_lg2_maxloops = FXP_FRAC_BITS_DEF + 1;
unsigned int FXP_one = 1u << FXP_FRAC_BITS_DEF;
int FXP_pow2_wpos_shift = 2 * FXP_INT_BITS_M1 \
                            - FXP_FRAC_BITS_DEF - 1;
int FXP_pow2_wneg_shift = 2 * FXP_INT_BITS_M1 \
                            - FXP_FRAC_BITS_DEF;

// Default desired max frac decimal value
// (can be changed dynamically calling set_[auto_]frac_max_dec
// TODO: Candidates to convert to ulongy?
// they were ints but made them longs since used only
// for decimal frac conversions, were longs are used.
long fxp_frac_max_dec = 9999;
long fxp_frac_max_dec_p1 = 10000;

// Single array of structs that have two uints (pseudo-ulongs,)
// instead of using two separate arrays of uints as before.
// (in fxp_l.c, this array is of course a simpler array of just
// ulongs.) Values here correspond to:
//                  a[k] = lg2(1 + 1/2^k)
// represented with two uints. The first uint is an fxp with
// 31 frac bits and one magnitude whole bit (notice, not a
// sign but a magnitude single whole bit, needed
// to be able to represent the first value in the table
// == 1.0) (Might optimize this later on removing this value
// and handling it implicitely in the code, while also
// shifting all other values for 1 extra bit of precision.)
// The 2nd uint is the continuation of the frac value
// for extra, long-like precision, 63 frac bits in total
static const ulongy FXP_BKM_LOGS_NEW[] = {
        {   0x80000000, 0x00000000 },   // index [0]
        {   0x4AE00D1C, 0xFDEB4400 },
        {   0x2934F097, 0x9A371600 },
        {   0x15C01A39, 0xFBD68800 },
        {   0xB31FB7D,  0x64898B00 },
        {   0x5AEB4DD,  0x63BF61C0 },
        {   0x2DCF2D0,  0xB85A4540 },
        {   0x16FE50B,  0x6EF08510 },
        {   0xB84E23,   0x6BD563B8 },
        {   0x5C3E0F,   0xFC29D594 },
        {   0x2E24CA,   0x6E87E8A8 },   // [10]
        {   0x1713D6,   0x2F7957C3 },
        {   0xB8A47,    0x6150DFE4 },
        {   0x5C53A,    0xC47E94D8 },
        {   0x2E2A3,    0x2762FA6B },
        {   0x17153,    0x05002E4A },
        {   0xB8A9,     0xDED47C11 },
        {   0x5C55,     0x067F6E58 },
        {   0x2E2A,     0x89050622 },
        {   0x1715,     0x45F3D72B },
        {   0xB8A,      0xA35640A7 },   // [20]
        {   0x5C5,      0x51C23599 },
        {   0x2E2,      0xA8E6E01E },
        {   0x171,      0x5474E163 },
        {   0xB8,       0xAA3ACD06 },
        {   0x5C,       0x551D7D98 },
        {   0x2E,       0x2A8EC491 },
        {   0x17,       0x154763BA },
        {   0xB,        0x8AA3B239 },
        {   0x5,        0xC551D933 },
        {   0x2,        0xE2A8EC9F },   // [30]
        {   0x1,        0x71547651 },
        {   0x0,        0xB8AA3B28 },
        {   0x0,        0x5C551D94 },
        {   0x0,        0x2E2A8ECA },
        {   0x0,        0x17154765 }    // [35]
            //0xB8AA3B2,        // <---- *
            //0x5C551D9,
            //0x2E2A8EC,
            //0x1715476,
            //0xB8AA3B, // [40]
            //0x5C551D,
            //0x2E2A8E,
            //0x171547,
            //0xB8AA3,
            //0x5C551,
            //0x2E2A8,
            //0x17154,
            //0xB8AA,
            //0x5C55,
            //0x2E2A,   // [50]
            //0x1715,
            //0xB8A,
            //0x5C5,
            //0x2E2,
            //0x171,
            //0xB8,
            //0x5C,
            //0x2E,
            //0x17,
            //0xB,      // [60]
            //0x5,
            //0x2,
            //0x1,      // [63]
            //0x0       // [64]
};
// Insteresting that starting with the row marked
// with the *, each entry is exactly a 4-bit
// right-shift of the value 4 positions earlier

// BKM One aligned for Array/Argument: unsigned fxp with 1 whole bit
static const unsigned int FXP_BKM_A_ONE = 1u << (FXP_INT_BITS - 1);
const ulongy FXP_BKM_A_ONE_ULONGY = { FXP_BKM_A_ONE, 0u };
// BKM One aligned for X: unsigned fxp with 2 whole bits
static const unsigned int FXP_BKM_X_ONE = 1u << (FXP_INT_BITS - 2);
const ulongy FXP_BKM_X_ONE_ULONGY = { FXP_BKM_X_ONE, 0u };
// Maxed Unsigned int, all bits set to 1
static const unsigned int FXP_BKM_MAXU = ~(0u);

// Auxiliary variables for the implementations that use longs (fxp_l.c)
const unsigned long FXP_BKM_A_ONE_L = \
                    ((unsigned long) FXP_BKM_A_ONE) << FXP_INT_BITS;

const unsigned long FXP_BKM_X_ONE_L = \
                    ((unsigned long) FXP_BKM_X_ONE) << FXP_INT_BITS;
unsigned long FXP_max_lshifted = (FXP_MAX_L << FXP_FRAC_BITS_DEF) \
                    | (((1 << FXP_FRAC_BITS_DEF) - 1) / 2);
int FXP_lg2_l_mshift    = 2 * FXP_INT_BITS - FXP_FRAC_BITS_DEF - 1;
int FXP_lg2_l_mshift_m1 = 2 * FXP_INT_BITS - FXP_FRAC_BITS_DEF - 2;
int FXP_lg2_l_mshift_p1 = 2 * FXP_INT_BITS - FXP_FRAC_BITS_DEF;

//static const unsigned int UINT_ALL_ONES = ~0u;
//static const unsigned int UINT_ALL_ONES_RS1 = ~0u >> 1;

/*
 * Given an fxp with x number of frac bits, returns
 * the rounded representation using y frac bits.
 * Used internally to adjust the unsigned
 * transcendental constants when changing the number
 * of frac bits to use.
 */
static inline unsigned int fxp_rshift_tconst(unsigned int fxp, int x, int y)
{
        int shift = x - y;
        if (shift <= 0) return (unsigned int) FXP_POS_INF;
        unsigned int rbit = (fxp >> (shift - 1)) & 1u;
        return (fxp >> shift) + rbit;
}

int fxp_get_whole_bits()
{
        return FXP_whole_bits;
}

int fxp_get_frac_bits()
{
        return FXP_frac_bits;
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
        FXP_frac_bits = (nfracbits < FXP_FRAC_BITS_MIN? \
                            FXP_FRAC_BITS_MIN:
                            (nfracbits > FXP_FRAC_BITS_MAX?
                                FXP_FRAC_BITS_MAX: nfracbits));
        FXP_frac_bits_m1 = FXP_frac_bits - 1;
        FXP_whole_bits = FXP_INT_BITS - FXP_frac_bits;
        FXP_whole_bits_m1 = FXP_whole_bits - 1;
        FXP_whole_bits_m2 = FXP_whole_bits - 2;

        // fxp_frac_mask should correspond to 2^FXP_FRAC_BITS - 1
        FXP_frac_mask = (1u << FXP_frac_bits) - 1;

        // When using all bits for frac (except for the sign bit),
        // then our max valid frac cannot be equal to frac_mask
        // (because in that case that value is already the largest
        // positive integer == POS_INF), so we must substract one
        // from the frac mask to get the largest valid frac value
        // in that case
        FXP_frac_max = FXP_frac_mask - (FXP_whole_bits == 1? 1: 0);
        FXP_frac_max_int = (int) FXP_frac_max;
        FXP_frac_max_neg = -FXP_frac_max_int;
        FXP_frac_max_p1 = FXP_frac_max + 1;

        // Auxiliary variables for fxp_l.c implementations
        FXP_max_lshifted = ((FXP_MAX_L) << FXP_frac_bits) \
                                | (FXP_frac_mask / 2);

        // Max and min valid values for the whole part of the fxp's
        FXP_whole_max = FXP_MAX >> FXP_frac_bits;
        FXP_whole_min = (-FXP_whole_max);
        FXP_whole_min_m1 = FXP_whole_min - 1;

        // Max and min floating-point conversion values
        float FXP_htiniest_f = (((float) 1uL) / \
                        (1uL << FXP_frac_bits)) / 2.0;

        FXP_max_f = (float) FXP_whole_max \
                    + ((float) ((unsigned long) FXP_frac_max)) \
                        / ((float) (1uL << (FXP_frac_bits)));

        double FXP_htiniest_d = (((double) 1uL) / \
                        (1uL << FXP_frac_bits)) / 2.0;

        FXP_max_d = (double) FXP_whole_max \
                    + ((double) ((unsigned long) FXP_frac_max)) \
                        / ((double) (1uL << (FXP_frac_bits)));

        FXP_htiniest_ld = (((long double) 1uL) / \
                        (1uL << FXP_frac_bits)) / 2.0L;

        FXP_max_ld = (long double) FXP_whole_max \
                    + ((long double) ((unsigned long) FXP_frac_max)) \
                        / ((long double) (1uL << (FXP_frac_bits)));

        FXP_min_f = -FXP_max_f;
        FXP_min_d = -FXP_max_d;
        FXP_min_ld = -FXP_max_ld;
        FXP_min_fx = FXP_min_f - FXP_htiniest_f;
        FXP_max_fx = FXP_max_f + FXP_htiniest_f;
        FXP_min_dx = FXP_min_d - FXP_htiniest_d;
        FXP_max_dx = FXP_max_d + FXP_htiniest_d;
        FXP_min_ldx = FXP_min_ld - FXP_htiniest_ld;
        FXP_max_ldx = FXP_max_ld + FXP_htiniest_ld;

        // Adjust precision of e, pi, etc. to the frac bits in use
        FXP_shifted_e = fxp_rshift_tconst(FXP_E_I32, \
                        FXP_INT_BITS_M2, FXP_frac_bits);
        FXP_shifted_pi = fxp_rshift_tconst(FXP_PI_I32, \
                        FXP_INT_BITS_M2, FXP_frac_bits);
        FXP_shifted_ln_2 = fxp_rshift_tconst(FXP_LN_2_I32, \
                        FXP_INT_BITS, FXP_frac_bits);
        FXP_shifted_lg10_2 = fxp_rshift_tconst(FXP_LG10_2_I32, \
                        FXP_INT_BITS, FXP_frac_bits);

        // Auxiliary variables used for lg2 and pow2 calculations
        FXP_one = 1 << FXP_frac_bits;
        //FXP_one_l = 1l << FXP_frac_bits;
        FXP_half = FXP_one >> 1;
        FXP_two = FXP_one << 1;
        FXP_lg2_maxloops = FXP_frac_bits + 1;
        FXP_lg2_l_mshift = 2 * FXP_INT_BITS - 1 - FXP_frac_bits;
        FXP_lg2_l_mshift_m1 = FXP_lg2_l_mshift - 1;
        FXP_lg2_l_mshift_p1 = FXP_lg2_l_mshift + 1;
        FXP_pow2_wpos_shift = FXP_INT_BITS_M1 + FXP_whole_bits_m1;
        FXP_pow2_wneg_shift = FXP_INT_BITS_M1 + FXP_whole_bits;

        return FXP_frac_bits;
}

/*
 * Automatically set the fractional max decimal to use.
 * The number of nines in fxp_frac_max_dec will be
 * floor(fxp_frac_bits / 4), e.g. one nine for every 4 bits
 */
int fxp_set_auto_frac_max_dec()
{
        int nnines = FXP_frac_bits / 4;
        fxp_frac_max_dec = 9;
        for (int i = 1; i < nnines; i++) {
                fxp_frac_max_dec = (fxp_frac_max_dec * 10) + 9;
        }
        fxp_frac_max_dec_p1 = fxp_frac_max_dec + 1;
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
    fxp_frac_max_dec_p1 = fxp_frac_max_dec + 1;
    return (int) fxp_frac_max_dec;
}

/*
 * Return the max decimal value 99...9 to be used as
 * max posible decimal fraction value
 */
inline int fxp_get_frac_max_dec()
{
    return (int) fxp_frac_max_dec;
}

/*
 * Create an fxp number given its whole part only
 */
inline int fxp(int whole)
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
        if (whole > FXP_whole_max)
                return FXP_POS_INF;
        if (whole < FXP_whole_min)
                return FXP_NEG_INF;
        if (bin_frac > FXP_frac_max_int) {
                bin_frac = FXP_frac_max;
        } else {
                if (bin_frac < FXP_frac_max_neg) {
                        bin_frac = FXP_frac_max_neg;
                }
        }
        int sign = 0;
        if ((whole == 0) && (bin_frac < 0)) {
                // Special case for negative numbers when whole part is
                // zero, then the fxp gets sign from the frac
                sign = 1;
                bin_frac = -bin_frac;
        } else {
                // All other cases: fxp sign == sign of whole part
                if (whole < 0) {
                        sign = 1;
                        whole = -whole;
                }
                if (bin_frac < 0) {
                        bin_frac = -bin_frac;
                }
        }
        int positive_fxp = (whole << FXP_frac_bits) | bin_frac;
        return sign? -positive_fxp: positive_fxp;
}

/*
 * Create an fxp number given its whole and (decimal) frac parts.
 * dec_frac is expected to be a decimal number between 0 and
 * fxp_frac_max_dec, e.g. between 0 and 999
 * Usage examples:
 *     For fxp=16.001, you would invoke: fxp_dec(16, 1)
 *     For 20.09: fxp_dec(24, 90)
 *     For 24.5:  fxp_dec(24, 500)
 * Note last example frac is not 5 but 500. (5 would correspond to
 * 24.005) This decimal value gets scaled into the binary range
 * available for frac.
 * For negative numbers with whole part=0, the frac value must be
 * negative.
 * If frac is too large it will get truncated to its most
 * significant decimal digits until under the value of
 * FXP_FRAC_MAX_DEC, e.g. for a max of 999, a frac=987654
 * will be truncated to 987
 * (Truncating and not rounding, because the latter would require
 * changing the whole part in some border cases)
 */
int fxp_dec(int whole, int dec_frac)
{
        if ((whole == FXP_UNDEF) || (dec_frac == FXP_UNDEF))
                return FXP_UNDEF;
        if (whole > FXP_whole_max)
                return FXP_POS_INF;
        if (whole < FXP_whole_min)
                return FXP_NEG_INF;
        int neg_frac = 0;
        if (dec_frac < 0) {
                neg_frac = 1;
                dec_frac = -dec_frac;
        }
        long trunc_frac = dec_frac;
        while (trunc_frac > fxp_frac_max_dec) {
                trunc_frac = trunc_frac / 10;
                //printf("   fxp_from_dec_frac: frac trimmed to: %d\n", \
                //    trunc_frac);
        }
        // Watch out this conversion itself can overflow when frac_bits
        // is large, and the frac_max_dec value is also large.
        // Using longs here because of this
        // TODO: we could use ulongy?
        int bin_frac = (int) ((trunc_frac * \
                                (long) FXP_frac_max_p1) \
                                    / fxp_frac_max_dec);
        return fxp_bin(whole, neg_frac? -bin_frac: bin_frac);
}

inline int fxp_get_whole_part(int fxp)
{
        if (fxp < 0)
                if (fxp <= FXP_NEG_INF)
                        // Technically, the whole part of -INF would
                        // still be -INF, and both the whole and frac
                        // parts of UNDEF would still be UNDEF. But here
                        // returning just what our whole-part values for
                        // -INF (and also UNDEF) actually correspond to.
                        // However, see also the comment in
                        // fxp_get_bin_frac()
                        return -FXP_whole_max;
                else
                        return -((-fxp) >> FXP_frac_bits);
        else
                return (fxp >> FXP_frac_bits);
}

/*
 * Get the frac part directly (binary)
 */
inline int fxp_get_bin_frac(int fxp)
{
        if (fxp < 0)
                if (fxp <= FXP_NEG_INF)
                        // In the whole-part function above we are
                        // returning the same whole parts for both
                        // -INF and UNDEF. Here for -INF we return the
                        // value used, but for UNDEF something
                        // different and odd in order to enable
                        // differenciating -INF from UNDEF when
                        // checking their whole and frac constituents.
                        // It must be something that no other fxp
                        // could ever have, not even the +/-INF
                        // values. That can be some (alleged) frac
                        // bits that are actually invalid, e.g.
                        // outside the range of valid frac bits, like
                        // -frac_max - 1. Notice that this value
                        // technically "overflows" the frac bits, but
                        // in a way makes sense to mark and
                        // differenciate precisely only UNDEF this way.
                        // Conveniently, when the # of frac bits is
                        // 31, this scheme returns the full FXP_UNDEF
                        // as the frac part.
                        return -FXP_frac_mask \
                                - (fxp == FXP_UNDEF? 1: 0);
                else
                        return -((-fxp) & FXP_frac_mask);
        else
                return (fxp & FXP_frac_mask);
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
        // TODO: use ulongy?
        long num, ldiv;
        long positive_frac;
        if (fxp < 0) {
                positive_frac = (-fxp) & FXP_frac_mask;
                num = -(((long) positive_frac) \
                            * fxp_frac_max_dec_p1);
        } else {
                positive_frac = fxp & FXP_frac_mask;
                num = ((long) positive_frac) \
                            * fxp_frac_max_dec_p1;
        }
        ldiv = num / FXP_frac_max_p1;
        return (int) ldiv;
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
        return (int) (p >> FXP_frac_bits);
}

int fxp_unsafe_div(int fxp1, int fxp2)
{
        long n1 = ((long) fxp1) << FXP_frac_bits;
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
 * Safe fxp1 - fxp2 just using the safe add functions :)
 */
inline int fxp_sub(int fxp1, int fxp2)
{
        // Conservative check to never attempt to change sign
        // to the true most negative int (aka. FXP_UNDEF)
        if (fxp2 == FXP_UNDEF) return FXP_UNDEF;
        return fxp_add(fxp1, -fxp2);
}

/*
 * Return number of bits used by x, that is,
 * most significant bit in x that is a 1
 * (returned value is between 0 and FXP_INT_BITS)
 */
inline int fxp_nbits(unsigned int x)
{
    if (x == 0) return 0;
    // Replacing with gcc builtin function that counts leading zeros
    return FXP_INT_BITS - __builtin_clz(x);
}

/*
 * Safe implementation of fxp multiplication using only ints.
 *
 * Works for systems in which sizeof(long) is not larger
 * than sizeof(int).
 * Does not use divisions.
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
        if (nba + nbb > FXP_whole_bits) {
                // The product will for sure overflow just by
                // multiplying the whole parts.
                // Return appropriately signed infinity
                //printf("01. Overflowing!!!!!\n");
                return ((fxp1 >= 0 && fxp2 >= 0) \
                            || (fxp1 < 0 && fxp2 < 0))?
                            FXP_POS_INF: FXP_NEG_INF;
        }

        ab = a * b;
        // Frac parts, and nbits in them
        x = fxp_get_bin_frac(v1);
        y = fxp_get_bin_frac(v2);
        nbx = fxp_nbits(x);
        nby = fxp_nbits(y);
        // Compute pf1, pf2, pf3, and pfsum
        ay = a * y;
        bx = b * x;
        //printf("ab: %d\nay: %d\nbx: %d\n", ab, ay, bx);
        int pf1 = fxp_get_bin_frac(ay);
        int pf2 = fxp_get_bin_frac(bx);
        int pf3, rbit;

        // Calculate pf3 using distributive multiplication
        unsigned int ux = ((unsigned int) x) << FXP_whole_bits;
        unsigned int uy = ((unsigned int) y) << FXP_whole_bits;
        unsigned int pf3u = dmul_into_uint(ux, uy);
        rbit = (pf3u >> FXP_whole_bits_m1) & 1u;
        pf3 = (int) ((pf3u >> FXP_whole_bits) + rbit);

        //printf("x:%X, ux:%X\n", x, ux);
        //printf("y:%X, uy:%X\n", y, uy);
        //printf("pf3u: %u, pf3:%X\n", pf3u, pf3);;
        // We must sum safely because depending on the fxp config
        // there might not be enough whole bits to hold the
        // carry on from just the sum of these three frac pieces
        int pfsum = fxp_add(pf1, fxp_add(pf2, pf3));
        //printf("pf1:%X\npf2:%X\npf3:%x\n", pf1, pf2, pf3);
        //printf("pfsum:%X\n", pfsum);
        if (pfsum == FXP_POS_INF) {
                // Overflow. Return appropriately signed infinity
                return ((fxp1 >= 0 && fxp2 >= 0) || (fxp1 < 0 && fxp2 < 0))?
                            FXP_POS_INF: FXP_NEG_INF;
        }
        int pffrac = fxp_get_bin_frac(pfsum);

        // Compute remaining whole parts (pw2, pw3, pw4)
        int pw1 = fxp_get_whole_part(ay);
        int pw2 = fxp_get_whole_part(bx);
        int pw3 = fxp_get_whole_part(pfsum);
        // Sum all whole parts safely
        int pwsum = ab + pw1 + pw2 + pw3;
        if (pwsum > FXP_whole_max) {
                // Overflow.
                // Return appropriately signed infinity
                return ((fxp1 >= 0 && fxp2 >= 0) \
                            || (fxp1 < 0 && fxp2 < 0))?
                            FXP_POS_INF: FXP_NEG_INF;
        }
        //printf("pwsum:%d, pfsum_frac:%d\n", pwsum, fxp_get_bin_frac(pfsum));
        // No overflow, return the appropriately signed product
        int pproduct = fxp_bin(pwsum, pffrac);
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
                return ((fxp1 > 0 && fxp2 > 0) || \
                            (fxp1 < 0 && fxp2 < 0))?
                                FXP_POS_INF: FXP_NEG_INF;

        //printf("x:%d, mul upper lim:%d\n", x, fxp_mul(FXP_MAX, y));
        if (y <= FXP_frac_mask) {
                int mulcheck = fxp_mul(FXP_MAX, y);
                //printf("FXP_MAX : %10d, (x%X)\n", FXP_MAX, FXP_MAX);
                //printf("y       : %10d, (x%X)\n", y, y);
                //printf("mulcheck: %10d, (x%X)\n", mulcheck, mulcheck);
                // Due to rounding and border cases, using > is not
                // enough here, we need >=, otherwise in some cases
                // expecting +INF the division can return UNDEF
                //if (x > mulcheck)
                if (x >= mulcheck)
                        return ((fxp1 > 0 && fxp2 > 0) || \
                                    (fxp1 < 0 && fxp2 < 0))?
                                        FXP_POS_INF: FXP_NEG_INF;
        }

        int bx = fxp_nbits(x);
        int by = fxp_nbits(y);

        int m, next_bit_mask;
        int bmax = bx + FXP_frac_bits;
        if (by == FXP_INT_BITS_M1) {
                // Border case, we need to shrink the divisor so that
                // the minuend can still have 1 more bit for the cases
                // when m < y
                y = (y >> 1);
                by--;
                // if (VERBOSE) printf("\n\t\tNew shrinked divisor %d\n\n", y);
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

inline unsigned int fxp_get_e()
{
        return FXP_shifted_e;
}

inline unsigned int fxp_get_pi()
{
        return FXP_shifted_pi;
}

inline unsigned int fxp_get_ln_2()
{
        return FXP_shifted_ln_2;
}

inline unsigned int fxp_get_lg10_2()
{
        return FXP_shifted_lg10_2;
}

/*
 * Internal auxiliary function that calculates the characteristic
 * and mantissa for lg2(fxp1), returning them separately in a
 * struct tuple, at maximum precision:
 * - The characteristic as int with 0 frac bits
 * - The mantissa now as a ulongy (emulated ulong using two uints)
 *   with the same alignment as the BKM array values:
 *      1st mantissa uint with 31 frac bits
 *      2nd mantissa uint with 32 frac bits
 * Uses the BKM L-Mode algorithm (requires table of pre-calculated
 * values) and only ints.
 *
 * For more details on BKM:
 * https://en.wikipedia.org/wiki/BKM_algorithm
 */
static inline struct tuple fxp_get_lg2_as_tuple(int fxp1, \
                                        const int MAX_LOOPS)
{
        struct tuple result;
        // Here fxp1 for sure > 0
        int clz = __builtin_clz((unsigned int) fxp1);
        // Notice clz > 0 since at least sign bit in fxp1 is 0
        int nbx = FXP_INT_BITS - clz;
        // Assign the characteristic
        result.ping = ((fxp1 < FXP_one) || (fxp1 >= FXP_two))?
                (nbx - FXP_frac_bits - 1): 0;
        // Here we replicate what the lg2_l implementation does,
        // but now using the 'ulongy' struct to emulate an ulong,
        // instead of handling two separate uints for each
        unsigned int aa = ((unsigned int) fxp1) << (clz - 1);
        ulongy argument = ulongy_create(aa, 0u);
        ulongy x = FXP_BKM_X_ONE_ULONGY;
        ulongy y = ULONGY_ZERO;
        ulongy xs, z;
        #ifdef VERBOSE
        printf("\nfxp_get_lg2_as_tuple: clz%d, argument: {x%X,%X}\n", \
                    clz, argument.hi, argument.lo);
        #endif
        for (unsigned int shift = 1; shift <= MAX_LOOPS; shift++) {
                //unsigned long xs = (x >> shift);
                xs = ulongy_rshift(x, shift);

                // z = x + xs;
                z = ulongy_add(x, xs);

                int c = ulongy_compare(z, argument);

                #ifdef VERBOSE
                printf("shift:%u looping\n", shift);
                printf("\txs      : {x%X,%X}\n", xs.hi, xs.lo);
                printf("\targument: {x%X,%X}\n", argument.hi, argument.lo);
                printf("\tz       : {x%X,%X}\n", z.hi, z.lo);
                printf("\tc:%d  ", c);
                #endif

                // if (z <= argument) {
                if (c <= 0) {
                        x = z;
                        // Here we use the precalculated values
                        // y += FXP_BKM_LOGS_L[shift];
                        y = ulongy_add(y, FXP_BKM_LOGS_NEW[shift]);
                        #ifdef VERBOSE
                        printf("\tUpdated x:{x%X,%X}  Updated y:{x%X,%X}", \
                                x.hi, x.lo, y.hi, y.lo);
                        #endif
                }
                #ifdef VERBOSE
                printf("\n");
                #endif
        }
        result.pong = y;
        return result;
}

/*
 * Default log2 using BKM and only ints.
 */
int fxp_lg2(int fxp1)
{
        if (fxp1 < 0) return FXP_UNDEF;
        if (fxp1 == 0) return FXP_NEG_INF;
        if (fxp1 == FXP_POS_INF) return FXP_POS_INF;
        // Get the separate characteristic and full mantissa
        struct tuple tup = fxp_get_lg2_as_tuple(fxp1, FXP_lg2_maxloops);
        // Shift the mantissa (rounding it) for current fxp config
        int lfr = ulongy_hi_to_uint_rounded(tup.pong);
        int rbit = (lfr >> FXP_whole_bits_m2) & 1u;
        int m = (lfr >> FXP_whole_bits_m1) + rbit;
        //printf("Final mantissa:%d  (rounding bit was %d)\n", m, rbit);
        // Return the complete logarithm (c + m) as fxp
        if (tup.ping < FXP_whole_min) {
                if ((tup.ping == FXP_whole_min_m1) \
                        && (m > 0)) {
                        // Special case: in spite of not enough whole
                        // bits available for characteristic, there are
                        // for the characteristic + 1, which is what
                        // gets returned here as whole part
                        return INT_MIN + m;
                }
                // Overflow: the characteristic of the logarithm
                // will not fit in the whole bits part of our
                // current fxp settings
                return FXP_NEG_INF;
        }
        return (tup.ping >= 0)?
                (tup.ping << FXP_frac_bits) + m:
                (-(-tup.ping << FXP_frac_bits)) + m;
}

static inline int neg_lg2_x_factor(tuple tup, const ulongy FACTOR)
{
        if ((tup.ping < FXP_whole_min_m1) \
                || ((tup.ping == FXP_whole_min_m1) \
                    && (ulongy_compare(tup.pong, ULONGY_ZERO) == 0))) {
                return FXP_NEG_INF;
        }
                #ifdef VERBOSE
                if (tup.ping== FXP_whole_min_m1) {
                        printf("\n\tCharacteristic right beyond the -LIMIT!\n\t");
                }
                #endif

        unsigned int posc = (-tup.ping);
        int shift_for_posc = __builtin_clz(posc) - 1;
        unsigned int shifted_posc = posc << shift_for_posc;
        ulongy cxf = dmul_ulongy_x_uint(FACTOR, shifted_posc);

                #ifdef VERBOSE
                printf("lg2 c : %d (x%X,  b",
                        tup.ping, tup.ping);
                print_int_as_bin(tup.ping, 0);
                printf(")\n\t");
                printf("lg2 m : ");
                print_ulongy_as_hex(tup.pong); printf("\n\t");
                print_ulongy_as_bin(tup.pong); printf(")\n\t");
                printf("FACTOR: ");
                print_ulongy_as_hex(FACTOR); printf("\n\t");
                printf("+Characteristic L-shifted (by %d), and factor:\n\t", \
                        shift_for_posc);
                print_uint_as_bin(shifted_posc);
                printf("\n\t");
                print_ulongy_as_bin(FACTOR); printf("\n");
                //printf(" (%LE)\n\t", ((long double) FACTOR) / (~0u));
                printf("\tcxf: ");
                print_ulongy_as_hex(cxf); printf("\n\t");
                print_ulongy_as_bin(cxf); printf("\n");
                printf("\t-cxf: ");
                print_ulongy_as_hex(ulongy_negate(cxf)); printf("\n\t");
                print_ulongy_as_bin(ulongy_negate(cxf)); printf("\n");

                printf("\tshift_for_c is %d\n", shift_for_posc);
                //printf(" (%LE)\n\t", ((long double) cxf) / (1u << shift_for_c));
                #endif

        //unsigned int mxf = mul_distrib_ulongies(tup.f, FACTOR);
        // mxf always non-negative
        ulongy mxf = dmul_ulongys(tup.pong, FACTOR);
        int shift_for_m = FXP_INT_BITS_M1 - shift_for_posc;
        ulongy shifted_mxf = ulongy_rshift_rounding(mxf, shift_for_m);
        ulongy sum = ulongy_sub(cxf, shifted_mxf);

        // Round and then final shift to leave result aligned
        // with current fxp configuration
        int final_rshift = FXP_whole_bits + shift_for_posc;
        int final_lg = -((int) ulongy_get_lo(
                            ulongy_rshift_rounding(sum, final_rshift)));

                #ifdef VERBOSE
                //printf("Mantissa and factor:\n\t");
                //print_ulongy_as_bin(tup.pong); printf("\n\t");
                //printf(" (%LE)\n\t", ((long double) tup.ping) / (~0u >> 1));
                //print_ulongy_as_bin(FACTOR); printf("\n\t");
                //printf(" (%LE)\n\t", ((long double) FACTOR) / (~0u));
                printf("\tmxf: ");
                print_ulongy_as_hex(mxf); printf("\n\t");
                print_ulongy_as_bin(mxf); printf("\n");
                //printf(" (%LE)\n\t", ((long double) mxf) / (~0u >> 1));
                printf("\tshift for m is %d\n\t", shift_for_m);
                printf("R-shifted mxf (by %d, to add it to cxf) is:\n\t", shift_for_m);
                print_ulongy_as_bin(shifted_mxf);
                //printf(" (%LE)\n", ((long double) shifted_mxf) / \
                //                        (~0u >> (1 + shift_for_m)));
                printf("\n\tpsum: ");
                print_ulongy_as_hex(sum); printf("\n\t");
                print_ulongy_as_bin(sum);
                //printf(" (%LE)\n\t", ((long double) sum) / (1u << shift_for_c));
                printf("\n\tfinal_lg: %X", final_lg);
                print_uint_as_bin(final_lg); printf("\n");
                //printf(" (%LE)\n", ((long double) final_lg) / FXP_frac_max_p1);
                #endif

        return final_lg;
}

static inline int pos_lg2_x_factor(tuple tup, const ulongy FACTOR)
{
        int shift_for_c;
        unsigned int shifted_c = 0;
        ulongy cxf;
        if (tup.ping > 0) {
                shift_for_c = __builtin_clz(tup.ping) - 1;
                shifted_c = tup.ping << shift_for_c;
                cxf = dmul_ulongy_x_uint(FACTOR, shifted_c);
        } else {
                // tup.ping == 0
                shift_for_c = FXP_INT_BITS_M1;
                shifted_c = 0;
                cxf = ULONGY_ZERO;
        }

                #ifdef VERBOSE
                printf("\n\tlg2 c : %d (x%X,  b",
                        tup.ping, tup.ping);
                print_int_as_bin(tup.ping, 0);
                printf(")\n\t");
                printf("lg2 m : x");
                print_ulongy_as_hex(tup.pong); printf("\n\t");
                printf("FACTOR: x");
                print_ulongy_as_hex(FACTOR); printf("\n\t");
                printf("Shifted c (by %d): x%X,  b", shift_for_c, shifted_c);
                print_uint_as_bin(shifted_c); printf("\n\t");
                //print_ulongy_as_bin(FACTOR); printf("\n");
                //printf(" (%LE)\n\t", ((long double) FACTOR) / (~0u));
                printf("cxf: x");
                print_ulongy_as_hex(cxf); printf("\n\tb");
                print_ulongy_as_bin(cxf); printf("\n\t");
                //printf(" (%LE)\n\t", ((long double) cxf) / (1u << shift_for_c));
                #endif

        //unsigned int mxf = mul_distrib_ulongies(tup.f, FACTOR);
        // mxf always non-negative
        ulongy mxf = dmul_ulongys(tup.pong, FACTOR);
        int shift_for_m = FXP_INT_BITS_M1 - shift_for_c;
        ulongy shifted_mxf = ulongy_rshift_rounding(mxf, shift_for_m);
        ulongy sum = ulongy_add(cxf, shifted_mxf);

        // Round and then final shift to leave result aligned
        // with current fxp configuration
        int final_rshift = FXP_whole_bits + shift_for_c;
        //int final_lg = ulongy_hi_to_uint_rounded(
        //                    ulongy_rshift_rounding(sum, final_rshift));
        int final_lg = ulongy_get_lo(
                            ulongy_rshift_rounding(sum, final_rshift));

                #ifdef VERBOSE
                //printf("Mantissa and factor:\n\t");
                //print_ulongy_as_bin(tup.pong); printf("\n\t");
                //printf(" (%LE)\n\t", ((long double) tup.ping) / (~0u >> 1));
                //print_ulongy_as_bin(FACTOR); printf("\n\t");
                //printf(" (%LE)\n\t", ((long double) FACTOR) / (~0u));
                printf("mxf: x");
                print_ulongy_as_hex(mxf); printf("\n\tb");
                print_ulongy_as_bin(mxf); printf("\n\t");
                //printf(" (%LE)\n\t", ((long double) mxf) / (~0u >> 1));
                printf("R-shifted mxf (by %d): x", shift_for_m);
                print_ulongy_as_hex(shifted_mxf); printf("\n\tb");
                print_ulongy_as_bin(shifted_mxf);
                //printf(" (%LE)\n", ((long double) shifted_mxf) / \
                //                        (~0u >> (1 + shift_for_m)));
                printf("\n\tSum: x");
                print_ulongy_as_hex(sum); printf("\n\tb");
                print_ulongy_as_bin(sum); printf("\n");
                printf("\tfinal_rshift: %d\n\t", final_rshift);
                printf("Final lg: x%X\n\tb", final_lg);
                print_int_as_bin(final_lg, 0);
                printf(" (%LE)\n", ((long double) final_lg) / FXP_frac_max_p1);
                #endif

        return final_lg;
}

/*
 * Calculates log2 and then multiplies by the given factor
 * Analogous to the one in fxp_l, but here using only ints.
 */
static inline int lg2_x_factor(int fxp1, const ulongy FACTOR)
{
                #ifdef VERBOSE
                printf("\n\tArgument is %d (x%X,  b",
                        fxp1, fxp1);
                print_int_as_bin(fxp1, 0);
                printf(")\n\t");
                printf("Whole: %d,  dec frac: %d (/%ld)\n\t", \
                        fxp_get_whole_part(fxp1),
                        fxp_get_dec_frac(fxp1),
                        fxp_frac_max_dec_p1);
                #endif
        struct tuple tup = fxp_get_lg2_as_tuple(fxp1, FXP_INT_BITS);
        if (tup.ping < 0)
                return neg_lg2_x_factor(tup, FACTOR);
        else
                return pos_lg2_x_factor(tup, FACTOR);
}

/*
 * Default implementation of ln() using lg2 and only ints
 */
int fxp_ln(int fxp1)
{
        if (fxp1 < 0) return FXP_UNDEF;
        if (fxp1 == 0) return FXP_NEG_INF;
        if (fxp1 == FXP_POS_INF) return FXP_POS_INF;
        return lg2_x_factor(fxp1, FXP_LN_2_ULONGY);
}

/*
 * Default implementation of lg10() using lg2 and only ints
 */
int fxp_lg10(int fxp1)
{
        if (fxp1 < 0) return FXP_UNDEF;
        if (fxp1 == 0) return FXP_NEG_INF;
        if (fxp1 == FXP_POS_INF) return FXP_POS_INF;
        return lg2_x_factor(fxp1, FXP_LG10_2_ULONGY);
}

#ifdef VERBOSE
static long double print_ulongy_as_ld(ulongy x, int nfbits)
{
    unsigned long long twopower = 1llu << nfbits;
    long double ldfrac = 0.0L;
    ulongy frac = x;
    while (ulongy_compare(frac, ULONGY_ZERO) > 0) {
        if (frac.lo & 1u) ldfrac += ((long double) 1.0L) / twopower;
        frac = ulongy_rshift(frac, 1);
        twopower = twopower >> 1;
    }
    return ldfrac;
}
#endif


static inline struct ulongy fxp_bkm_emode(ulongy argument)
{
        #ifdef VERBOSE
        printf("\nbkm_emode running with argument {%X,%X}\n", \
                    argument.hi, argument.lo);
        #endif
        ulongy x = FXP_BKM_X_ONE_ULONGY;
        ulongy y = ULONGY_ZERO;
        ulongy z, xs;
        for (unsigned int k = 0; k < FXP_INT_BITS; k++) {
                z = ulongy_add(y, FXP_BKM_LOGS_NEW[k]);
                #ifdef VERBOSE
                printf("k:%d,  z:{%X,%X}\n", k, z.hi, z.lo);
                #endif
                if (ulongy_compare(z, argument) <= 0) {
                        y = z;
                        xs = ulongy_rshift(x, k);
                        x = ulongy_add(x, xs);
                        #ifdef VERBOSE
                        printf("\tUpdating y {%X,%X} und x {%X,%X}\n", \
                                        y.hi, y.lo, x.hi, x.lo);
                        #endif
                }
        }
        #ifdef VERBOSE
        printf("bkm_emode returning x = {%X,%X}  (%.19Lf)\n", \
                        x.hi, x.lo, print_ulongy_as_ld(x, 62));
        #endif
        return x;
}

// pow2 calculation (2^fxp1) of a NON-negative argument,
// using the BKM (E-mode) algorithm and only ints.
// Argument comes as a {w, f} tuple, where the f component
// is a ulongy, and must already have the same alignment
// as the BKM array values (1 whole bit)
static inline int fxp_pow2_wpos_tuple(tuple tup)
{
        if (tup.ping >= FXP_whole_bits_m1) return FXP_POS_INF;
        ulongy argument = tup.pong;
        //printf("pow2_wpos argument is {%X,%X} (clz: %d)\n",\
        //            argument.hi, argument.lo, \
        //            __builtin_clz(argument.hi));
        ulongy x = fxp_bkm_emode(argument);
        //int shift = FXP_INT_BITS_M1 + FXP_whole_bits_m1 - tup.ping;
        int shift = FXP_pow2_wpos_shift - tup.ping;
        //printf("pow2_wpos shift is: %d\n", shift);
        ulongy shifted = ulongy_rshift_rounding(x, shift);
        return shifted.lo;
}

// pow2 calculation (2^fxp1) of a negative argument,
// using the BKM (E-mode) algorithm and ints.
// Argument comes as a {w,f} tuple, where the
// f component must already have the same alignment
// as the BKM array values (1 whole bit)
static inline int fxp_pow2_wneg_tuple(struct tuple tup)
{
        if (tup.ping > FXP_frac_bits) return 0;
        ulongy argument = ulongy_sub(FXP_BKM_A_ONE_ULONGY, tup.pong);
        //printf("pow2_wneg argument is {%X,%X} (clz: %d)\n",\
        //            argument.hi, argument.lo, \
        //            __builtin_clz(argument.hi));
        ulongy x = fxp_bkm_emode(argument);
        //int shift = FXP_INT_BITS_M1 + FXP_whole_bits + tup.ping;
        int shift = FXP_pow2_wneg_shift + tup.ping;
        //printf("pow2_wneg shift is: %d\n", shift);
        ulongy shifted = ulongy_rshift_rounding(x, shift);
        return (int) shifted.lo;
}

// pow2 calculation (2^fxp1) using the BKM (E-mode) algorithm
// Analogous to implementation in bkm.c, but here tailored
// for fxp's using only ints
int fxp_pow2(int fxp1)
{
        if (fxp1 == FXP_UNDEF) return FXP_UNDEF;
        if (fxp1 == FXP_NEG_INF) return 0;
        if (fxp1 == FXP_POS_INF) return FXP_POS_INF;
        #ifdef VERBOSE
        printf("\n\npow2 running for x%X (int = %d)\n", fxp1, fxp1);
        #endif
        if (fxp1 >= 0) {
                ulongy frac = {
                        ((unsigned int) fxp_get_bin_frac(fxp1)) \
                                    << FXP_whole_bits_m1, 0u };
                tuple tup = { fxp_get_whole_part(fxp1), frac };
                return fxp_pow2_wpos_tuple( tup );
        } else {
                int pfxp1 = -fxp1;
                ulongy frac = {
                        ((unsigned int) fxp_get_bin_frac(pfxp1)) \
                                << FXP_whole_bits_m1, 0u };
                tuple tup = { fxp_get_whole_part(pfxp1), frac };
                return fxp_pow2_wneg_tuple( tup );
        }
}


// Internal function to calculate the product of
// x times C and return it as an {int, ulongy} tuple
// with the ulongy frac part having the same alignment
// as the BKM array values (1 whole bit)
// This function is basically a line-to-line equivalent
// of the mirror function in fxp_l.c for longs, but here
// using ulongys
static inline tuple get_xc_as_tuple(int x, \
                                    const ulongy C, \
                                    int c_nwbits)
{
        ulongy fc, wfc, ffc, ffmask;
        unsigned int sf;
        int margin = c_nwbits + 1;

        // First calculating frac(x) * C
        // fraction part of x, l-shifted all the way
        // so we leave no whole bits
        sf = ((unsigned int) fxp_get_bin_frac(x)) \
                                    << FXP_whole_bits_m1;
        // sf times the factor C
        fc = dmul_ulongy_x_uint(C, sf);
        // whole and frac parts of that product f * C
        ffmask = ulongy_rshift(ULONGY_ALL_ONES, margin);
        wfc = ulongy_rshift(fc, FXP_LONG_BITS - margin);
        ffc = ulongy_bitwise_and(fc, ffmask);

                #ifdef VERBOSE
                printf("\n\tf (1 whole bits) and factor C (%d whole bits):\n\t", \
                                c_nwbits);
                print_uint_as_bin(sf);
                printf(" (%LE)\n\t", ((long double) sf / (1ul << FXP_INT_BITS_M1)));
                print_ulongy_as_bin(C);
                //printf(" (%LE)\n\t", ((long double) C) \
                //                        / (1ul << FXP_INT_BITS - c_nwbits));
                printf("\n\tProduct fxC (%d whole bits) is:\n\t", margin);
                print_ulongy_as_bin(fc);
                //printf(" (%LE)\n\t", ((long double) fc) \
                //                        / (1ul << FXP_INT_BITS - c_nwbits));
                printf("\n\tw_fxC: {%X,%X}", wfc.hi, wfc.lo);
                printf("\n\tf_fxC (%d whole bits):\n\t", margin);
                print_ulongy_as_bin(ffc);
                //printf(" (%LE)\n\t", ((long double) ffc) \
                //                        / (1u << FXP_INT_BITS - c_nwbits));
                #endif

        // Now calculating whole(x) * C
        ulongy wc, wfmask, wwc, fwc;
        unsigned int wx, swx;
        int wx_nbits, wx_clz, wx_clz_m1, w_margin;
        // wx whole part of x l-shifted all the way
        wx = fxp_get_whole_part(x);
        wx_clz = (wx == 0)? FXP_INT_BITS: __builtin_clz(wx);
        wx_nbits = FXP_INT_BITS - wx_clz;
        // swx whole(x) l-shifted all the way
        swx = wx << wx_clz;
        // times the factor C
        wc = dmul_ulongy_x_uint(C, swx);

        // Whole and frac parts of that product wc
        w_margin = wx_nbits + c_nwbits;
        wfmask = ulongy_rshift(ULONGY_ALL_ONES, w_margin);
        wwc = ulongy_rshift(wc, FXP_LONG_BITS - w_margin);
        fwc = ulongy_bitwise_and(wc, wfmask);

                #ifdef VERBOSE
                printf("\n\tw (%d whole bits) and factor C (%d whole bits):\n\t", \
                            wx_nbits, c_nwbits);
                print_uint_as_bin(swx);
                printf(" (%LE)\n\t", ((long double) swx / (1u << wx_clz)));
                print_ulongy_as_bin(C);
                //printf(" (%LE)\n\t", ((long double) C) \
                //                        / (1u << FXP_INT_BITS - c_nwbits));
                printf("\n\tProduct wxC (%d whole bits) is:\n\t", w_margin);
                print_ulongy_as_bin(wc);
                //printf(" (%LE)\n\t", ((long double) wc) \
                //                        / (1u << wx_clz_m1));
                printf("\n\tw_wxC: {%u,%u}\n\t", wwc.hi, wwc.lo);
                print_ulongy_as_bin(wwc);
                printf("\n\t");
                printf("f_wxC (%d whole bits):\n\t", w_margin);
                print_ulongy_as_bin(fwc);
                printf("\n");
                //printf(" (%LE)\n\t", ((long double) fwc) \
                //                        / (1ul << FXP_INT_BITS - w_margin));
                #endif

        // Left-align both f_fxC and f_wxC leaving 1 whole bit
        // in order to add them up
        ffc = ulongy_lshift(ffc, margin - 1);
        fwc = ulongy_lshift(fwc, w_margin - 1);

        // Get final whole and frac parts, and return as tuple
        ulongy fplusf = ulongy_add(ffc, fwc);
        ulongy w_fplusf = ulongy_rshift(fplusf, FXP_LONG_BITS_M1);
        ulongy vping_sum = ulongy_add( ulongy_add(wwc, wfc), w_fplusf );
        int vping = vping_sum.lo;
        struct tuple xc = { vping, \
                            ulongy_bitwise_and(fplusf, \
                                        ULONGY_ALL_ONES_RS1) };

                #ifdef VERBOSE
                printf("\n\tvping is: %d", vping);
                printf("\n\tfinal wsum:\n\t");
                print_uint_as_bin(xc.ping);
                printf(" (%u)\n\t", xc.ping);
                printf("final fsum:\n\t");
                print_ulongy_as_bin(xc.pong);
                printf("\n");
                //printf(" (%LE)\n", ((long double) xc.f) \
                //                        / (1u << FXP_INT_BITS_M1));
                #endif

        return xc;
}

#ifdef VERBOSE
static void print_tuple(char * msg, struct tuple tup)
{
        unsigned long pong = \
                ((unsigned long) tup.pong.hi << FXP_INT_BITS) \
                        | tup.pong.lo;
        long double num = (long double) tup.ping \
                + ((long double) pong) / (~0ul >> 1);
        printf("%s: ", msg);
        printf("{x%X, {x%X,x%X}} == %.10Lf\n", \
                tup.ping, tup.pong.hi, tup.pong.lo, num);
}
#endif


// Calculate the pow2 of a NON-negative fxp argument x
// multiplied by a factor C (with C having the
// indicated number of whole bits)
static inline int fxp_pow2_pos_arg_xfactor( \
                        int x, \
                        const ulongy C, \
                        int c_nwbits)
{
        tuple xc = get_xc_as_tuple(x, C, c_nwbits);
                #ifdef VERBOSE
                print_tuple("Pos tuple is: ", xc);
                #endif
        return fxp_pow2_wpos_tuple( xc );
}

// Calculate the pow2 of a negative fxp argument x
// multiplied by a factor C (with C having the
// indicated number of whole bits)
static inline int fxp_pow2_neg_arg_xfactor( \
                        int x, \
                        const ulongy C, \
                        int c_nwbits)
{
        tuple xc = get_xc_as_tuple(x, C, c_nwbits);
                #ifdef VERBOSE
                print_tuple("Neg tuple is: ", xc);
                #endif
        return fxp_pow2_wneg_tuple( xc );
}

// Implementation of exp(x) using pow2() and lg2():
// e^x  == 2^( lg2(e^x) )
//      == 2^( x * lg2(e) )
int fxp_exp(int fxp1)
{
        if (fxp1 == FXP_UNDEF) return FXP_UNDEF;
        if (fxp1 == FXP_NEG_INF) return 0;
        if (fxp1 == FXP_POS_INF) return FXP_POS_INF;
        return (fxp1 >= 0)?
                fxp_pow2_pos_arg_xfactor(fxp1, \
                                    FXP_LG2_E_ULONGY, \
                                    FXP_LG2_E_WBITS):
                fxp_pow2_neg_arg_xfactor( -fxp1, \
                                    FXP_LG2_E_ULONGY, \
                                    FXP_LG2_E_WBITS);
}

// Implementation of pow10(x) using pow2() and lg2():
// 10^x == 2^( lg2(10^x) )
//      == 2^( x * lg2(10) )
//
int fxp_pow10(int fxp1)
{
        if (fxp1 == FXP_UNDEF) return FXP_UNDEF;
        if (fxp1 == FXP_NEG_INF) return 0;
        if (fxp1 == FXP_POS_INF) return FXP_POS_INF;
        return (fxp1 >= 0)?
                fxp_pow2_pos_arg_xfactor(fxp1, \
                                    FXP_LG2_10_ULONGY, \
                                    FXP_LG2_10_WBITS):
                fxp_pow2_neg_arg_xfactor( -fxp1, \
                                    FXP_LG2_10_ULONGY, \
                                    FXP_LG2_10_WBITS);
}

/*
// Implementation of powxy(x) using pow2() and lg2():
// x^y  == 2^( lg2(x^y) )
//      == 2^( y * lg2(x) )
int fxp_powxy(int x, int y)
{
        if ((x < 0) || (x == FXP_UNDEF) || (y == FXP_UNDEF))
                return FXP_UNDEF;
        if (x == 0) {
                // lg2(x) is -INF
                if (y < 0) return FXP_POS_INF;  // 2^(+INF)
                if (y > 0) return 0;            // 2^(-INF)
                return FXP_UNDEF;               // 2^(0 * -INF)
        }
        // x > 0
        // TODO

        return 0;
}

// Square root implementation
int fxp_sqrt(int fxp1)
{
        // TODO
        return 0;
}

*/
