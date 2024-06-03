
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
//#define VERBOSE_BKM 1
#ifdef VERBOSE
#include "print_as_bits.h"
#include "fxp_aux.h"
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
// As a minimum at least 4 frac bits, for at least 1
// decimal digit of precision
const int FXP_FRAC_BITS_MIN = 4;

// # of loops to use for different calculations
// (control the achievable precision)
int FXP_lg2_maxloops = FXP_FRAC_BITS_DEF;           // for lg2
const int FXP_LOGX_LOOPS = FXP_INT_BITS_M1 - 1;     // for ln and lg10
const int FXP_POWX_LOOPS = FXP_INT_BITS_M1;         // for pow2, exp, and pow10
const int FXP_SQRT_LG_LOOPS = FXP_INT_BITS_M2;      // for lg2 in sqrt
int FXP_sqrt_pw_loops = FXP_FRAC_BITS_DEF + 1 \
            + (FXP_INT_BITS - FXP_FRAC_BITS_DEF)/2; // for pow2 in sqrt
int FXP_cordic_loops = FXP_FRAC_BITS_DEF + 2;       // for cordic

int FXP_sqrt_cordic_loops = 11;                      // for the implementation of sqrt using cordic

const int FXP_POWXY_LG_LOOPS = FXP_INT_BITS + 3;    // for lg2 in powxy
const int FXP_POWXY_PW_LOOPS = FXP_INT_BITS + 1;    // for pow2 in powxy

const int FXP_WORD_BITS = FXP_INT_BITS >> 1;
const int FXP_WORD_BITS_M1 = FXP_WORD_BITS - 1;
const unsigned int FXP_RWORD_MASK = ((1u << FXP_WORD_BITS) - 1);
const unsigned int FXP_LWORD_MASK = FXP_RWORD_MASK \
                                    << FXP_WORD_BITS;
const unsigned long FXP_RINT_MASK = ((1uL << FXP_INT_BITS) - 1);
const unsigned long FXP_LINT_MASK = FXP_RINT_MASK \
                                    << FXP_INT_BITS;
const int FXP_LONG_BITS = ((int) sizeof(long)) * 8;
const int FXP_LONG_BITS_M1 = FXP_LONG_BITS - 1;
const int FXP_LONG_BITS_M2 = FXP_LONG_BITS - 2;

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
int FXP_frac_bits_p1 = FXP_FRAC_BITS_DEF + 1;

// Default number of bits for the whole part (includes sign bit)
int FXP_whole_bits = FXP_INT_BITS - FXP_FRAC_BITS_DEF;
int FXP_whole_bits_m1 = FXP_INT_BITS_M1 - FXP_FRAC_BITS_DEF;
int FXP_whole_bits_m2 = FXP_INT_BITS_M1 - FXP_FRAC_BITS_DEF - 1;

int FXP_int_plus_whole_bits_m1 = 2*FXP_INT_BITS - FXP_FRAC_BITS_DEF - 1;
int FXP_int_plus_whole_bits_m2 = 2*FXP_INT_BITS - FXP_FRAC_BITS_DEF - 2;

int FXP_long_mfrac_bits = FXP_LONG_BITS - FXP_FRAC_BITS_DEF;

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

//int FXP_rsqrt_1p5 = (11u << (FXP_FRAC_BITS_DEF - 1));

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

const unsigned int UINT_ALL_ONES = ~0u;
const unsigned int UINT_ALL_ONES_RS1 = UINT_ALL_ONES >> 1;

// A super_fxp is not just a "bigger FXP," it also has its own
// fxp configuration (it's own number of whole vs. frac bits,)
// independent of the global fxp setting, always carrying this
// individual configuration within.
typedef struct super_fxp {
        int sign;
        int nwbits;     // <- Number of whole bits
        ulongy number;  // <- The super fxp number
} super_fxp;

// Aux struct used internally for lg2
typedef struct lg2tuple {
        int characteristic;
        ulongy mantissa;
} lg2tuple;

// A tuple of ulongy's
typedef struct fxptuple_ul {
        ulongy a;
        ulongy b;
} tuple_ul;

const fxptuple FXP_TUPLE_UNDEFS = {FXP_UNDEF, FXP_UNDEF};

const unsigned int SFXP_MAX_WBITS = FXP_LONG_BITS - FXP_FRAC_BITS_MIN;

// transcendental constants
// e = 2.718281828...
const unsigned int FXP_E_I32 = 0x56FC2A2Cu;     // <- 29 frac bits
const unsigned int FXP_E_I32_X = 0x515DA54Du;

// pi = 3.14159265...
// ----*----1----*----2----*----3-- --*----4----*----5----*----6----
// 01100100100001111110110101010001 00010000101101000110000100011010011000100110001100110001010001011100000
const unsigned int FXP_PI_I32 = 0x6487ED51u;    // <- 29 frac bits
const unsigned int FXP_PI_I32_X = 0x10B4611Au;
const unsigned long FXP_PI_L = (((unsigned long) FXP_PI_I32) << FXP_INT_BITS) \
                                | FXP_PI_I32_X;

// transcendental constants
// lg2(10) = 3.32192809...
const unsigned int FXP_LG2_10_I32 = 0x6A4D3C25u;
const unsigned int FXP_LG2_10_I32_X = 0xE68DC57Fu;
const ulongy FXP_LG2_10_ULONGY = {FXP_LG2_10_I32, FXP_LG2_10_I32_X};
const super_fxp SFXP_LG2_10_FACTOR = {0, 3, FXP_LG2_10_ULONGY};

// lg2(e) = 1.44269504... (1 whole, 31 and 63 frac bits)
const unsigned int FXP_LG2_E_I32 = 0x5C551D94u;
const unsigned int FXP_LG2_E_I32_X = 0xAE0BF85Eu;
const ulongy FXP_LG2_E_ULONGY = {FXP_LG2_E_I32, FXP_LG2_E_I32_X};
const super_fxp SFXP_LG2_E_FACTOR = {0, 2, FXP_LG2_E_ULONGY};

// ln(2) = 0.69314718...
const unsigned int FXP_LN_2_I32 = 0x58B90BFBu;
const unsigned int FXP_LN_2_I32_X = 0xE8E7BCD6u;
const ulongy FXP_LN_2_ULONGY = {FXP_LN_2_I32, FXP_LN_2_I32_X};
const super_fxp SFXP_LN_2_FACTOR = {0, 1, FXP_LN_2_ULONGY};

// lg10(2) = 0.30102999...
const unsigned int FXP_LG10_2_I32 = 0x268826A1u;
const unsigned int FXP_LG10_2_I32_X = 0x3EF3FDE6u;
const ulongy FXP_LG10_2_ULONGY = {FXP_LG10_2_I32, FXP_LG10_2_I32_X};
const super_fxp SFXP_LG10_2_FACTOR = {0, 1, FXP_LG10_2_ULONGY};

// Auxiliary variables used in the lg2 and cordic implementations
int FXP_one = 1 << FXP_FRAC_BITS_DEF;
int FXP_minus_one = -(1 << FXP_FRAC_BITS_DEF);
int FXP_half = 1 << (FXP_FRAC_BITS_DEF - 1);
int FXP_quarter = 1 << (FXP_FRAC_BITS_DEF - 2);
int FXP_two = 1 << (FXP_FRAC_BITS_DEF + 1);
int FXP_almost1 = (1 << FXP_FRAC_BITS_DEF) - 1;

long FXP_one_l = 1L << (FXP_INT_BITS + FXP_FRAC_BITS_DEF);
long FXP_half_l = 1L << (FXP_INT_BITS_M1 + FXP_FRAC_BITS_DEF);
long FXP_quarter_l = 1L << (FXP_INT_BITS_M2 + FXP_FRAC_BITS_DEF);
long FXP_two_l = 1L << (FXP_INT_BITS + 1 + FXP_FRAC_BITS_DEF);

/*
// Auxiliary magic number for inverse square root
const unsigned int FXP_RSQRT_MAGIC_I32 = 0x5FE6EB50;
const unsigned int FXP_RSQRT_MAGIC_I32_X = 0xC7B537A9;
//const unsigned int FXP_RSQRT_MAGIC_I32 = 0x5F375A86;
//const unsigned int FXP_RSQRT_MAGIC_I32_X = 0x0;
const ulongy FXP_RSQRT_MAGIC_ULONGY = {FXP_RSQRT_MAGIC_I32, FXP_RSQRT_MAGIC_I32_X};
const super_fxp SFXP_RSQRT_MAGIC = {0, 1, FXP_RSQRT_MAGIC_ULONGY};
*/

// Default desired max frac decimal value
// (can be changed dynamically calling set_[auto_]frac_max_dec
// TODO: Candidates to convert to ulongy?
// they were ints but made them longs since used only
// for decimal frac conversions, were longs are used.
long fxp_frac_max_dec = 9999L;
long fxp_frac_max_dec_p1 = 10000L;

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
        { 0x80000000u, 0x00000000u },     // [0]
        { 0x4AE00D1Cu, 0xFDEB43CFu },
        { 0x2934F097u, 0x9A3715FCu },
        { 0x15C01A39u, 0xFBD6879Fu },
        { 0x0B31FB7Du, 0x64898B3Eu },
        { 0x05AEB4DDu, 0x63BF61CCu },
        { 0x02DCF2D0u, 0xB85A4531u },
        { 0x016FE50Bu, 0x6EF08517u },
        { 0x00B84E23u, 0x6BD563BAu },
        { 0x005C3E0Fu, 0xFC29D593u },
        { 0x002E24CAu, 0x6E87E8A8u },     // [10]
        { 0x001713D6u, 0x2F7957C3u },
        { 0x000B8A47u, 0x6150DFE4u },
        { 0x0005C53Au, 0xC47E94D8u },
        { 0x0002E2A3u, 0x2762FA6Bu },
        { 0x00017153u, 0x05002E4Au },
        { 0x0000B8A9u, 0xDED47C11u },
        { 0x00005C55u, 0x067F6E58u },
        { 0x00002E2Au, 0x89050622u },
        { 0x00001715u, 0x45F3D72Bu },
        { 0x00000B8Au, 0xA35640A7u },     // [20]
        { 0x000005C5u, 0x51C23599u },
        { 0x000002E2u, 0xA8E6E01Eu },
        { 0x00000171u, 0x5474E163u },
        { 0x000000B8u, 0xAA3ACD06u },
        { 0x0000005Cu, 0x551D7D98u },
        { 0x0000002Eu, 0x2A8EC491u },
        { 0x00000017u, 0x154763BAu },
        { 0x0000000Bu, 0x8AA3B239u },
        { 0x00000005u, 0xC551D933u },
        { 0x00000002u, 0xE2A8EC9Fu },     // [30]
        { 0x00000001u, 0x71547651u },
        { 0x00000000u, 0xB8AA3B28u },
        { 0x00000000u, 0x5C551D94u },
        { 0x00000000u, 0x2E2A8ECAu },
        { 0x00000000u, 0x17154765u },
        { 0x00000000u, 0x0B8AA3B2u },
        { 0x00000000u, 0x05C551D9u },     // *
        { 0x00000000u, 0x02E2A8ECu },
        { 0x00000000u, 0x01715476u },
        { 0x00000000u, 0x00B8AA3Bu }      // [40]
        /*
        { 0x00000000u, 0x005C551Du },
        { 0x00000000u, 0x002E2A8Eu },
        { 0x00000000u, 0x00171547u },
        { 0x00000000u, 0x000B8AA3u },
        { 0x00000000u, 0x0005C551u },
        { 0x00000000u, 0x0002E2A8u },
        { 0x00000000u, 0x00017154u },
        { 0x00000000u, 0x0000B8AAu },
        { 0x00000000u, 0x00005C55u },
        { 0x00000000u, 0x00002E2Au },     // [50]
        { 0x00000000u, 0x00001715u },
        { 0x00000000u, 0x00000B8Au },
        { 0x00000000u, 0x000005C5u },
        { 0x00000000u, 0x000002E2u },
        { 0x00000000u, 0x00000171u },
        { 0x00000000u, 0x000000B8u },
        { 0x00000000u, 0x0000005Cu },
        { 0x00000000u, 0x0000002Eu },
        { 0x00000000u, 0x00000017u },
        { 0x00000000u, 0x0000000Bu },     // [60]
        { 0x00000000u, 0x00000005u },
        { 0x00000000u, 0x00000002u },
        { 0x00000000u, 0x00000001u },
        { 0x00000000u, 0x00000000u },
        */
};
// Interesting that starting with the row marked with the *
// ( corresponds to [37]), each entry is exactly a 4-bit
// right-shift of the value 4 positions earlier

// BKM One aligned for X: unsigned fxp with 2 whole bits
const unsigned int FXP_BKM_X_ONE = 1u << FXP_INT_BITS_M2;
// BKM One aligned for Array/Argument: unsigned fxp with 1 whole bit
const unsigned int FXP_BKM_A_ONE = 1u << FXP_INT_BITS_M1;
const unsigned int FXP_BKM_A_POINT5 = FXP_BKM_A_ONE >> 1;

const ulongy FXP_BKM_X_ONE_ULONGY = { FXP_BKM_X_ONE, 0u };
const ulongy FXP_BKM_A_ONE_ULONGY = { FXP_BKM_A_ONE, 0u };
const ulongy FXP_BKM_A_POINT5_ULONGY = { FXP_BKM_A_POINT5, 0u };

// Auxiliary variables for the implementations that use longs (fxp_l.c)
const unsigned long FXP_BKM_X_ONE_L = 1ul << FXP_LONG_BITS_M2;
const unsigned long FXP_BKM_A_ONE_L = 1ul << FXP_LONG_BITS_M1;
const unsigned long FXP_BKM_A_POINT5_L = FXP_BKM_A_ONE_L >> 1;

unsigned long FXP_max_lshifted = (FXP_MAX_L << FXP_FRAC_BITS_DEF) \
                    | (((1 << FXP_FRAC_BITS_DEF) - 1) / 2);

static const int FXP_DEF_SHIFT = FXP_INT_BITS - FXP_FRAC_BITS_DEF;

// Initializations for FXP_FRAC_BITS_DEF = 16 bits
int FXP_shifted_e = 0x2B7E1;
int FXP_shifted_pi = 0x3243F;
int FXP_shifted_phalfpi = 0x19220;
int FXP_shifted_nhalfpi = -0x19220;
long FXP_shifted_pi_l = 0x3243F6A8885A3L;
long FXP_shifted_phalfpi_l = 0x1921FB54442D2L;
long FXP_shifted_nhalfpi_l = -0x1921FB54442D2L;

static const ulongy FXP_CORDIC_ANGLES[] = {
{ 0x3243F6A8u, 0x885A308Du },
{ 0x1DAC6705u, 0x61BB4F68u },
{ 0x0FADBAFCu, 0x96406EB1u },
{ 0x07F56EA6u, 0xAB0BDB71u }, // 4 entries
{ 0x03FEAB76u, 0xE59FBD39u },
{ 0x01FFD55Bu, 0xBA97624Au },
{ 0x00FFFAAAu, 0xDDDB94D5u },
{ 0x007FFF55u, 0x56EEEA5Cu }, // 8
{ 0x003FFFEAu, 0xAAB7776Eu },
{ 0x001FFFFDu, 0x5555BBBBu },
{ 0x000FFFFFu, 0xAAAAADDEu },
{ 0x0007FFFFu, 0xF555556Fu },
{ 0x0003FFFFu, 0xFEAAAAABu },
{ 0x0001FFFFu, 0xFFD55555u },
{ 0x0000FFFFu, 0xFFFAAAAAu },
{ 0x00007FFFu, 0xFFFF5555u }, // 16
{ 0x00003FFFu, 0xFFFFEAAAu },
{ 0x00001FFFu, 0xFFFFFD55u },
{ 0x00000FFFu, 0xFFFFFFAAu },
{ 0x000007FFu, 0xFFFFFFF5u },
{ 0x000003FFu, 0xFFFFFFFEu },
{ 0x00000200u, 0x00000000u },
{ 0x00000100u, 0x00000000u },
{ 0x00000080u, 0x00000000u }, // 24
{ 0x00000040u, 0x00000000u },
{ 0x00000020u, 0x00000000u },
{ 0x00000010u, 0x00000000u },
{ 0x00000008u, 0x00000000u },
{ 0x00000004u, 0x00000000u },
{ 0x00000002u, 0x00000000u },
{ 0x00000001u, 0x00000000u },
{ 0x00000000u, 0x80000000u }  // 32 entries
/*
{ 0x00000000u, 0x40000000u },
{ 0x00000000u, 0x20000000u },
{ 0x00000000u, 0x10000000u },
{ 0x00000000u, 0x08000000u },
{ 0x00000000u, 0x04000000u },
{ 0x00000000u, 0x02000000u },
{ 0x00000000u, 0x01000000u },
{ 0x00000000u, 0x00800000u },
{ 0x00000000u, 0x00400000u },
{ 0x00000000u, 0x00200000u },
{ 0x00000000u, 0x00100000u },
{ 0x00000000u, 0x00080000u },
{ 0x00000000u, 0x00040000u },
{ 0x00000000u, 0x00020000u },
{ 0x00000000u, 0x00010000u },
{ 0x00000000u, 0x00008000u },   // 48
{ 0x00000000u, 0x00004000u },
{ 0x00000000u, 0x00002000u },
{ 0x00000000u, 0x00001000u },
{ 0x00000000u, 0x00000800u },
{ 0x00000000u, 0x00000400u },
{ 0x00000000u, 0x00000200u },
{ 0x00000000u, 0x00000100u },
{ 0x00000000u, 0x00000080u },
{ 0x00000000u, 0x00000040u },
{ 0x00000000u, 0x00000020u },
{ 0x00000000u, 0x00000010u },
{ 0x00000000u, 0x00000008u },
{ 0x00000000u, 0x00000004u },
{ 0x00000000u, 0x00000002u },
{ 0x00000000u, 0x00000001u },
{ 0x00000000u, 0x00000000u }    // 64
*/
};

static const ulongy FXP_CORDIC_KFACTOR = { 0x26DD3B6Au, 0x10D7969Au };

#ifdef VERBOSE
void print_sfxp(char * msg, super_fxp x)
{
        printf("\n%s", msg);
        printf("Super_fxp  : {x%X, ", x.nwbits);
        print_ulongy_as_hex(x.number);
        printf("}\n");
}

static long double ulongy_as_ld(ulongy x, int nfbits)
{
    unsigned long twopower = 1ul << nfbits;
    long double ldfrac = 0.0L;
    ulongy frac = x;
    while (ulongy_compare(frac, ULONGY_ZERO) > 0) {
        if (frac.lo & 1u) ldfrac += ((long double) 1.0L) / twopower;
        frac = rshift_ulongy(frac, 1);
        twopower = twopower >> 1;
    }
    return ldfrac;
}

void print_lg2tuple(char * msg, lg2tuple tup)
{
        long double num = (long double) tup.characteristic \
                + ulongy_as_ld(tup.mantissa, FXP_LONG_BITS_M1);
        printf("%s: ", msg);
        printf("{x%X, ", tup.characteristic);
        print_ulongy_as_hex(tup.mantissa);
        printf("} == %.10Lf\n", num);
}
#endif

static inline super_fxp sfxp_create(unsigned int sign, \
                                    unsigned int nwbits, \
                                    ulongy num)
{
        nwbits = (nwbits > SFXP_MAX_WBITS)? SFXP_MAX_WBITS: nwbits;
        super_fxp x = { (sign != 0), nwbits, num };
        return x;
}

static inline int sfxp_2_fxp(super_fxp x)
{
        int sfb = FXP_LONG_BITS - x.nwbits;
        int fbdiff = sfb - FXP_frac_bits;
        if (x.sign) {
                if (fbdiff >= 0) {
                        // As many or more frac bits than in fxp's
                        unsigned int rbit = rshift_ulongy(x.number, fbdiff - 1).lo & 1u;
                        ulongy num = ulongy_add_uint(rshift_ulongy(x.number, fbdiff), rbit);
                        ulongy whole = rshift_ulongy(num, FXP_frac_bits);
                        if ((whole.lo > FXP_whole_max) || (whole.hi > 0))
                                return FXP_NEG_INF;
                        return -((int) num.lo);
                } else {
                        // Fewer frac bits than in fxp's
                        ulongy whole = rshift_ulongy(x.number, sfb);
                        if ((whole.lo > FXP_whole_max) || (whole.hi > 0))
                                return FXP_NEG_INF;
                        int frac = (int) (x.number.lo & (UINT_ALL_ONES >> x.nwbits - FXP_INT_BITS));
                        frac <<= -fbdiff;
                        return -fxp_bin((int) whole.lo, frac);
                }
        } else {
                if (fbdiff >= 0) {
                        // As many or more frac bits than in fxp's
                        unsigned int rbit = rshift_ulongy(x.number, fbdiff - 1).lo & 1u;
                        ulongy num = ulongy_add_uint(rshift_ulongy(x.number, fbdiff), rbit);
                        ulongy whole = rshift_ulongy(num, FXP_frac_bits);
                        if ((whole.lo > FXP_whole_max) || (whole.hi > 0))
                                return FXP_POS_INF;
                        return (int) num.lo;
                } else {
                        // Fewer frac bits than in fxp's
                        ulongy whole = rshift_ulongy(x.number, sfb);
                        if ((whole.lo > FXP_whole_max) || (whole.hi > 0))
                                return FXP_POS_INF;
                        int frac = (int) (x.number.lo & (UINT_ALL_ONES >> x.nwbits - FXP_INT_BITS));
                        frac <<= -fbdiff;
                        return fxp_bin((int) whole.lo, frac);
                }
        }
}

static inline unsigned int sfxp_get_poswhole(super_fxp x)
{
        return (x.nwbits)? rshift_ulongy(x.number, FXP_LONG_BITS - x.nwbits).lo: 0;
}

static inline ulongy get_sfxp_frac_for_bkme(super_fxp x)
{
        return rshift_ulongy(lshift_ulongy(x.number, x.nwbits), 1);
}

static inline super_fxp posfxp_x_possfxp(int x, super_fxp c)
{
        int clz_m1 = __builtin_clz(x) - 1;
        unsigned int shifted_x = x << clz_m1;
        ulongy product = dmul_ulongy_x_uint(c.number, x);
        int nwb = c.nwbits + FXP_whole_bits;
        super_fxp super_product = sfxp_create(0, nwb, product);
        return super_product;
}

static inline super_fxp get_fxp_x_sfxp(int x, super_fxp c)
{
        super_fxp product;
        if (x >= 0) {
                product = posfxp_x_possfxp(x, c);
                product.sign = c.sign;
        } else {
                product = posfxp_x_possfxp(-x, c);
                product.sign = !c.sign;
        }
        return product;
}

/*
 * R-shifts an unsigned int rounding the last bit
 */
inline unsigned int rshift_uint_rounding(unsigned int x, \
                                         unsigned int shift)
{
        unsigned int rbit = shift? (x >> (shift - 1)) & 1u: 0;
        return (x >> shift) + rbit;
}

/*
 * Given an fxp with x number of frac bits, returns
 * the rounded representation using y frac bits.
 * Used internally to adjust the unsigned
 * transcendental constants when changing the number
 * of frac bits to use.
 */
static inline unsigned int fxp_rshift_tconst(unsigned int fxp, int x, int y)
{
        if (x < y) return (unsigned int) FXP_POS_INF;
        return rshift_uint_rounding(fxp, x - y);
}

static inline unsigned long fxp_rshift_tconst_l(unsigned long fxp_l, int x, int y)
{
        if (x < y) return ((unsigned long) FXP_POS_INF) << FXP_INT_BITS;
        return rshift_ulong_rounding(fxp_l, x - y);
}

inline unsigned long rshift_ulong_rounding( \
                                    unsigned long n,
                                    unsigned int shift)
{
        unsigned long rbit = (!shift)? 0: (n >> (shift - 1)) & 1ul;
        return (n >> shift) + rbit;
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
        FXP_frac_bits_p1 = FXP_frac_bits + 1;
        FXP_whole_bits = FXP_INT_BITS - FXP_frac_bits;
        FXP_whole_bits_m1 = FXP_whole_bits - 1;
        FXP_whole_bits_m2 = FXP_whole_bits - 2;
        FXP_int_plus_whole_bits_m1 = FXP_INT_BITS + FXP_whole_bits_m1;
        FXP_int_plus_whole_bits_m2 = FXP_INT_BITS + FXP_whole_bits_m2;
        FXP_long_mfrac_bits = FXP_LONG_BITS - FXP_frac_bits;

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

        // Adjust shifted e, pi, etc. to the frac bits in use
        FXP_shifted_e = fxp_rshift_tconst(FXP_E_I32, \
                        FXP_INT_BITS - 3, FXP_frac_bits);
        FXP_shifted_pi = fxp_rshift_tconst(FXP_PI_I32, \
                        FXP_INT_BITS - 3, FXP_frac_bits);
        FXP_shifted_phalfpi = fxp_rshift_tconst(FXP_PI_I32, \
                        FXP_INT_BITS_M2, FXP_frac_bits);
        FXP_shifted_nhalfpi = -FXP_shifted_phalfpi;
        FXP_shifted_pi_l = fxp_rshift_tconst_l(FXP_PI_L, \
                        FXP_LONG_BITS - 3, FXP_INT_BITS + FXP_frac_bits);
        FXP_shifted_phalfpi_l = fxp_rshift_tconst_l(FXP_PI_L, \
                        FXP_LONG_BITS_M2, FXP_INT_BITS + FXP_frac_bits);
        FXP_shifted_nhalfpi_l = -FXP_shifted_phalfpi_l;

        // Auxiliary variables used for lg2 and pow2 calculations
        FXP_half = 1 << (FXP_frac_bits - 1);
        FXP_quarter = FXP_half >> 1;
        FXP_almost1 = (int) ((1u << FXP_frac_bits) - 1);
        if (FXP_whole_bits == 1) {
                FXP_one = FXP_POS_INF;
                FXP_two = FXP_POS_INF;
                FXP_minus_one = FXP_NEG_INF;
                //FXP_rsqrt_1p5 = FXP_POS_INF;
        } else {
                FXP_one = 1 << FXP_frac_bits;
                FXP_two = (FXP_whole_bits >= 3)? FXP_one << 1: FXP_POS_INF;
                FXP_minus_one = -FXP_one;
                //FXP_rsqrt_1p5 = (11u << (FXP_frac_bits - 1));
        }
        FXP_one_l = ((long) FXP_one) << FXP_INT_BITS;
        FXP_half_l = ((long) FXP_half) << FXP_INT_BITS;
        FXP_quarter_l = ((long) FXP_quarter) << FXP_INT_BITS;
        FXP_two_l = ((long) FXP_two) << FXP_INT_BITS;
        FXP_lg2_maxloops = FXP_frac_bits;
        FXP_sqrt_pw_loops = FXP_frac_bits + FXP_whole_bits/2 + 1;
        FXP_cordic_loops = FXP_frac_bits + 2;
        // Minimum number of sqrt cordic loops needed to avoid
        // asserts in the tester program (fxp_tester) with its
        // default WDELTA_MAX = 3.0L.
        // Note that at least 3 whole bits required for this implementation
        // fracbits     Loops needed
        // 29           15
        // 28           14
        // 27           14
        // 26           14
        // 25           14
        // 24           13
        // 23           13
        // 22           13
        // 21           13
        // 20           13
        // 19           12
        // 18           12
        // 17           12
        // 16           12
        // 15           11
        // 14           11
        // 13           11
        // 12           10
        // 11           10
        // 10           10
        //  9           10
        //  8            9
        //  7            9
        //  6            9
        //  5            8
        //  4            8
        // The following formula closely matches that table
        FXP_sqrt_cordic_loops = (FXP_frac_bits - 1)/4 + \
                (FXP_frac_bits==20 || FXP_frac_bits==16? 9: 8);

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
                        // fxp_get_frac_part_bin()
                        return -FXP_whole_max;
                else
                        return -((-fxp) >> FXP_frac_bits);
        else
                return (fxp >> FXP_frac_bits);
}

/*
 * Get the frac part directly (binary)
 */
inline int fxp_get_frac_part_bin(int fxp)
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
int fxp_get_frac_part_dec(int fxp)
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

inline unsigned int fxp_get_lshifted_frac(unsigned int fxp1)
{
        return (unsigned int) ((fxp1 & FXP_frac_mask) << FXP_whole_bits_m1);
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
        x = fxp_get_frac_part_bin(v1);
        y = fxp_get_frac_part_bin(v2);
        nbx = fxp_nbits(x);
        nby = fxp_nbits(y);
        // Compute pf1, pf2, pf3, and pfsum
        ay = a * y;
        bx = b * x;
        //printf("ab: %d\nay: %d\nbx: %d\n", ab, ay, bx);
        int pf1 = fxp_get_frac_part_bin(ay);
        int pf2 = fxp_get_frac_part_bin(bx);
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
        int pffrac = fxp_get_frac_part_bin(pfsum);

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
        //printf("pwsum:%d, pfsum_frac:%d\n", pwsum, fxp_get_frac_part_bin(pfsum));
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

        #ifdef VERBOSE
        printf("x:%d, mul upper lim:%d\n", x, fxp_mul(FXP_MAX, y));
        #endif
        if (y <= FXP_frac_mask) {
                int mulcheck = fxp_mul(FXP_MAX, y);
                #ifdef VERBOSE
                printf("FXP_MAX : %10d, (x%X)\n", FXP_MAX, FXP_MAX);
                printf("y       : %10d, (x%X)\n", y, y);
                printf("mulcheck: %10d, (x%X)\n", mulcheck, mulcheck);
                #endif
                // Due to rounding and border cases, using > is not
                // enough here, we need >=, otherwise in some cases
                // expecting +INF the division can return UNDEF
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
                #ifdef VERBOSE
                printf("\n\t\tNew shrinked divisor %d\n\n", y);
                #endif
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
        #ifdef VERBOSE
        printf("\td: starting m: %d (%d bits)\n", m, bidx);
        int loops = 0;
        int oldm = 0, ba = 0;
        #endif
        int q = 0;
        // bidx is (from left to right) the highest bit # we are
        // currently processing, so we loop till bidx exceeds the
        // right-most bit in the full (left-shifted) dividend
        while (bidx <= bmax) {
                #ifdef VERBOSE
                printf("\tm: %d (x%X)  bidx:%d\n", m, m, bidx);
                oldm = m;
                ba = (m >= y);
                #endif
                if (m >= y) {
                        // Append a 1 to the right of the quotient
                        q = (q << 1) | 1;
                        m = m - y;
                } else {
                        // Append a 0 to the right of the quotient
                        q = (q << 1);
                }
                #ifdef VERBOSE
                trace_fxp_div("div:", loops, FXP_frac_bits,
                                bidx, x, y, q, ba, oldm, (oldm - ba * y));
                loops ++;
                #endif
                bidx++;
                if (next_bit_mask > 0) {
                        // Pull down next bit from the dividend,
                        // and append it to the right of the minuend
                        m = (m << 1) | \
                                    ((x & next_bit_mask) >> shift);
                        next_bit_mask = (next_bit_mask >> 1);
                        shift--;
                } else {
                        // Pull down a 0 bit and append it to minuend
                        m = (m << 1);
                }
        }
        // Return properly signed final quotient
        int finalq = ((fxp1 >= 0 && fxp2 > 0) || (fxp1 < 0 && fxp2 < 0))? q: -q;
        #ifdef VERBOSE
        printf("Final quotient is: %d (x%X)\n", finalq, finalq);
        #endif
        return finalq;
}

inline unsigned int fxp_get_e()
{
        return FXP_shifted_e;
}

inline unsigned int fxp_get_pi()
{
        return FXP_shifted_pi;
}

inline unsigned int fxp_get_halfpi()
{
        return FXP_shifted_phalfpi;
}

inline unsigned long fxp_get_pi_l()
{
        return FXP_shifted_pi_l;
}

inline unsigned long fxp_get_halfpi_l()
{
        return FXP_shifted_phalfpi_l;
}

static inline ulongy fxp_bkm_lmode(unsigned int argument, \
                                   const int MAX_LOOPS)
{
        ulongy x = FXP_BKM_X_ONE_ULONGY;
        ulongy y = ULONGY_ZERO;
        ulongy xs, z;
        #ifdef VERBOSE_BKM
        printf("bkm_lmode: argument: {x%X}\n", argument);
        #endif
        for (int shift = 1; shift <= MAX_LOOPS; shift++) {
                //unsigned long xs = (x >> shift);
                xs = rshift_ulongy(x, shift);

                // z = x + xs;
                z = ulongy_add(x, xs);

                int c = ulongy_compare_to_uint(z, argument);

                #ifdef VERBOSE_BKM
                printf("bkm_lmode shift:%u looping\n", shift);
                printf("\txs      : {x%X,%X}\n", xs.hi, xs.lo);
                printf("\targument: {x%X}\n", argument);
                printf("\tz       : {x%X,%X}\n", z.hi, z.lo);
                printf("\tc:%d  ", c);
                #endif

                // if (z <= argument) {
                if (c <= 0) {
                        x = z;
                        // Here we use the precalculated values
                        // y += FXP_BKM_LOGS_L[shift];
                        y = ulongy_add(y, FXP_BKM_LOGS_NEW[shift]);
                        #ifdef VERBOSE_BKM
                        printf("\tUpdated x:{x%X,%X}  Updated y:{x%X,%X}", \
                                x.hi, x.lo, y.hi, y.lo);
                        #endif
                }
                #ifdef VERBOSE_BKM
                printf("\n");
                #endif
        }
        #ifdef VERBOSE
        printf("bkm_lmode final mantissa is: ");
        print_ulongy_as_hex(y); printf("\n");
        #endif
        return y;
}

/*
 * Internal auxiliary function that calculates the characteristic
 * and mantissa for lg2(fxp1), returning them separately in a
 * tuple struct, at maximum precision:
 * - The characteristic as int with 0 frac bits
 * - The mantissa as a ulongy
 *   with the same alignment as the BKM array values: 63 frac bits
 * Uses the BKM L-Mode algorithm (requires table of pre-calculated
 * values) and only ints.
 *
 * For more details on BKM:
 * https://en.wikipedia.org/wiki/BKM_algorithm
 */
static inline lg2tuple fxp_get_lg2tuple(int fxp1, \
                                const int MAX_LOOPS)
{
        lg2tuple result;
        // Here fxp1 for sure > 0
        int clz = __builtin_clz((unsigned int) fxp1);
        // Notice clz > 0 since at least sign bit in fxp1 is 0
        int nbx = FXP_INT_BITS - clz;
        // Assign the characteristic
        result.characteristic = ((fxp1 <= FXP_almost1) || (fxp1 >= FXP_two))?
                (nbx - FXP_frac_bits_p1): 0;
        #ifdef VERBOSE
        printf("\nlg2_as_tuple: fxp1: x%X  clz: %d,  nbx: %d\n", fxp1, clz, nbx);
        printf("lg2_as_tuple: characteristic is : %d\n", result.characteristic);
        #endif
        // Here we replicate what the lg2_l implementation does,
        // just not using longs but ulongys
        unsigned int argument = ((unsigned int) fxp1) << (clz - 1);
        result.mantissa = fxp_bkm_lmode(argument, MAX_LOOPS);
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
        lg2tuple tup = fxp_get_lg2tuple(fxp1, FXP_lg2_maxloops);

        // Shift the mantissa (rounding it) for current fxp config
        ulongy r = rshift_ulongy(tup.mantissa, FXP_int_plus_whole_bits_m2);
        unsigned int rbit = r.lo & 1u;
        r = rshift_ulongy(r, 1);
        int m = (int) ulongy_add_uint(r, rbit).lo;

        // Return the complete logarithm (c + m) as fxp
        if (tup.characteristic < FXP_whole_min) {
                if ((tup.characteristic == FXP_whole_min_m1) \
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
        return (tup.characteristic >= 0)?
                (tup.characteristic << FXP_frac_bits) + m:
                (-(-tup.characteristic << FXP_frac_bits)) + m;
}

static inline super_fxp sfxp_from_cm(int c, \
                                     ulongy m)
{
        #ifdef VERBOSE
        printf("c: %d, m: {%X,%X}\n", c, m.hi, m.lo);
        printf("characteristic:  "); print_int_as_bin(c, 0); printf("\n");
        printf("mantissa      :  "); print_ulongy_as_bin(m); printf("\n");
        #endif
        if (c >= 0) {
                int nb = fxp_nbits(c);
                ulongy num = ulongy_create(0u, c);
                num = lshift_ulongy(num, FXP_LONG_BITS_M1 - nb);
                unsigned int rbit = rshift_ulongy(m, (nb - 1)).lo & 1u;
                m = ulongy_add_uint(rshift_ulongy(m, nb), rbit);
                num = ulongy_add(num, m);
                super_fxp x = { 0, nb + 1, num };
                #ifdef VERBOSE
                printf("sfxp num      :  "); print_ulongy_as_bin(num); printf("\n");
                #endif
                return x;
        } else {
                int pc = -c;
                int nb = fxp_nbits(pc);
                ulongy num = ulongy_create(0u, pc);
                num = lshift_ulongy(num, FXP_LONG_BITS_M1 - nb);
                unsigned int rbit = rshift_ulongy(m, (nb - 1)).lo & 1u;
                m = ulongy_add_uint(rshift_ulongy(m, nb), rbit);
                num = ulongy_sub(num, m);
                super_fxp x = { 1, nb + 1, num };
                #ifdef VERBOSE
                printf("sfxp num      : -"); print_ulongy_as_bin(num); printf("\n");
                #endif
                return x;
        }
}

static inline super_fxp fxp_get_lg2_as_sfxp(int fxp1, \
                                        const int MAX_LOOPS)
{
        int clz = __builtin_clz((unsigned int) fxp1);
        int nbx = FXP_INT_BITS - clz;
        // Assign the characteristic
        int c = ((fxp1 <= FXP_almost1) || (fxp1 >= FXP_two))? \
                        (nbx - FXP_frac_bits_p1): 0;
        #ifdef VERBOSE
        printf("\nlg2_as_sfxp: fxp1: x%X  clz: %d,  nbx: %d\n", fxp1, clz, nbx);
        printf("lg2_as_sfxp: characteristic is : %d\n", c);
        #endif
        // Prepare argument for bkm_lmode
        unsigned int argument = ((unsigned int) fxp1) << (clz - 1);
        // Get the mantissa with bkm_lmode
        ulongy m = fxp_bkm_lmode(argument, MAX_LOOPS);
        // Return result as a super_fxp
        return sfxp_from_cm(c, m);
}

static inline super_fxp get_sfxp_x_sfxp(super_fxp x, super_fxp y)
{
        // Sign of the product assuming arguments are never both negative
        super_fxp prod = { x.sign | y.sign, \
                           x.nwbits + y.nwbits,
                           dmul_ulongys(x.number, y.number) };
        return prod;
}

/*
 * Calculates log2 and then multiplies by the given factor
 * Analogous to the one in fxp_l, but here using only ints.
 */
static inline int lg2_x_factor(int fxp1, super_fxp factor)
{
        super_fxp lg2x = fxp_get_lg2_as_sfxp(fxp1, FXP_LOGX_LOOPS);
        super_fxp prod = get_sfxp_x_sfxp(factor, lg2x);
        return sfxp_2_fxp(prod);
}

/*
 * Default implementation of ln() using lg2 and only ints
 */
int fxp_ln(int fxp1)
{
        if (fxp1 < 0) return FXP_UNDEF;
        if (fxp1 == 0) return FXP_NEG_INF;
        if (fxp1 == FXP_POS_INF) return FXP_POS_INF;
        return lg2_x_factor(fxp1, SFXP_LN_2_FACTOR);
}

/*
 * Default implementation of lg10() using lg2 and only ints
 */
int fxp_lg10(int fxp1)
{
        if (fxp1 < 0) return FXP_UNDEF;
        if (fxp1 == 0) return FXP_NEG_INF;
        if (fxp1 == FXP_POS_INF) return FXP_POS_INF;
        return lg2_x_factor(fxp1, SFXP_LG10_2_FACTOR);
}

static inline ulongy fxp_bkm_emode(ulongy argument, \
                                   const int MAX_LOOPS)
{
        #ifdef VERBOSE_BKM
        printf("\nbkm_emode running with argument ");
        print_ulongy_as_hex(argument); printf("\n");
        #endif
        ulongy x = FXP_BKM_X_ONE_ULONGY;
        ulongy y = ULONGY_ZERO;
        ulongy z, xs;
        // TODO: check optimizing here possibly starting at k = 1
        for (unsigned int k = 0; k < MAX_LOOPS; k++) {
                z = ulongy_add(y, FXP_BKM_LOGS_NEW[k]);
                #ifdef VERBOSE_BKM
                printf("k:%d,  z:{%X,%X}\n", k, z.hi, z.lo);
                #endif
                if (ulongy_compare(z, argument) <= 0) {
                        y = z;
                        xs = rshift_ulongy(x, k);
                        x = ulongy_add(x, xs);
                        #ifdef VERBOSE_BKM
                        printf("\tUpdating y {%X,%X} und x {%X,%X}\n", \
                                        y.hi, y.lo, x.hi, x.lo);
                        #endif
                }
        }
        #ifdef VERBOSE_BKM
        printf("bkm_emode returning x = {%X,%X}  (%.19Lf)\n", \
                        x.hi, x.lo, ulongy_as_ld(x, FXP_LONG_BITS_M2));
        #endif
        return x;
}

/*
 * pow2 calculation (2^fxp1) of a NON-negative argument,
 * using the BKM (E-mode) algorithm and only ints.
 * The frac ulongy argument must already have the same alignment
 * as the BKM array values (1 whole bit)
 * Input value for 2^x here is x = whole + frac
 */
static inline int fxp_pow2_wpos(int whole, \
                                ulongy frac, \
                                const int MAX_LOOPS)
{
        #ifdef VERBOSE
        printf("pow2_wpos for x%X.{%X,%X}\n", whole, frac.hi, frac.lo);
        #endif
        if (whole >= FXP_whole_bits_m1) return FXP_POS_INF;
        if ((whole == 0) && (ulongy_compare(frac, ULONGY_ZERO) == 0))
                return FXP_one;
        // Argument frac will be in [0, 1)
        ulongy x = fxp_bkm_emode(frac, MAX_LOOPS);
        // When one of the operands is a shifted 1, we don't really
        // need to call the (expensive) dmul operation, since the
        // result will be identical to the other operand anyway,
        // just shifted appropriately
        int shift = FXP_int_plus_whole_bits_m2 - whole;
        int shifted = (int) rshift_ulongy_rounding(x, shift).lo;
        #ifdef VERBOSE
        printf("x from bkm_e: {x%X,%X}\n", x.hi, x.lo);
        print_ulongy_as_bin(x); printf("\n");
        printf("shift is %d\n", shift);
        printf("pow2_wpos result is ");
        print_uint_as_bin(shifted); printf("\n");
        #endif
        return (shifted != FXP_UNDEF)? shifted: FXP_POS_INF;

}

// pow2 calculation (2^fxp1) of a negative argument,
// using the BKM (E-mode) algorithm and ints.
// The frac argument must already have the same alignment
// as the BKM array values (1 whole bit)
static inline int fxp_pow2_wneg(int whole, \
                                ulongy frac, \
                                const int MAX_LOOPS)
{
        #ifdef VERBOSE
        printf("pow2_wneg for -x%X.{%X,%X}\n", whole, frac.hi, frac.lo);
        #endif
        if (whole > FXP_frac_bits) return 0;
        if ((whole == 0) && (ulongy_compare(frac, ULONGY_ZERO) == 0))
                return FXP_one;
        // Notice argument a for bkm_emode will be in (0, 1]
        ulongy a = ulongy_sub(FXP_BKM_A_ONE_ULONGY, frac);
        ulongy x = fxp_bkm_emode(a, MAX_LOOPS);
        int shift = FXP_int_plus_whole_bits_m1 + whole;
        int shifted = (int) rshift_ulongy_rounding(x, shift).lo;
        #ifdef VERBOSE
        printf("x from bkm_e: {x%X,%X}\n", x.hi, x.lo);
        print_ulongy_as_bin(x); printf("\n");
        printf("shift is %d\n", shift);
        printf("pow2_wneg result is ");
        print_uint_as_bin(shifted); printf("\n");
        #endif
        return (shifted != FXP_UNDEF)? shifted: FXP_POS_INF;
}

static inline ulongy get_fxp_frac_for_bkme(int frac)
{
        #ifdef VERBOSE
        printf("get_fxp_frac_for_bkme: Fracbits is %d\n", FXP_frac_bits);
        printf("Frac is x%X (d: %d)\n", frac, frac);
        printf("Shift for frac is: %d\n", FXP_whole_bits_m1);
        #endif
        ulongy frac4bkme = { ((unsigned int) frac) << FXP_whole_bits_m1, 0ul };
        return frac4bkme;
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
                int w = fxp_get_whole_part(fxp1);
                ulongy bkmearg = get_fxp_frac_for_bkme( \
                                        fxp_get_frac_part_bin(fxp1));
                return fxp_pow2_wpos( w, bkmearg, FXP_POWX_LOOPS );
        } else {
                int pfxp1 = -fxp1;
                int w = fxp_get_whole_part(pfxp1);
                ulongy bkmearg = get_fxp_frac_for_bkme( \
                                        fxp_get_frac_part_bin(pfxp1));
                return fxp_pow2_wneg( w, bkmearg, FXP_POWX_LOOPS );
        }
}

// Calculate the pow2 of a NON-negative fxp argument x
// multiplied by a factor c
static inline int fxp_pow2_pos_arg_xfactor(int x, \
                                        super_fxp factorc, \
                                        const int MAX_LOOPS)
{
        super_fxp xc = get_fxp_x_sfxp(x, factorc);
        #ifdef VERBOSE
        print_sfxp("Pos xc (sfxp) is: ", xc);
        #endif
        int w = (int) sfxp_get_poswhole(xc);
        ulongy bkmearg = get_sfxp_frac_for_bkme(xc);
        return fxp_pow2_wpos( w, bkmearg, MAX_LOOPS );
}

// Calculate the pow2 of a negative fxp argument x
// multiplied by a factor C
static inline int fxp_pow2_neg_arg_xfactor(int x, \
                                        super_fxp factorc, \
                                        const int MAX_LOOPS)
{
        super_fxp xc = get_fxp_x_sfxp(x, factorc);
        #ifdef VERBOSE
        print_sfxp("Neg xc (sfxp) is: ", xc);
        #endif
        int w = (int) sfxp_get_poswhole(xc);
        ulongy bkmearg = get_sfxp_frac_for_bkme(xc);
        return fxp_pow2_wneg( w, bkmearg, MAX_LOOPS );
}

// Implementation of exp(x) using pow2() and lg2():
// e^x  == 2^( x * lg2(e) )
int fxp_exp(int fxp1)
{
        if (fxp1 == FXP_UNDEF) return FXP_UNDEF;
        if (fxp1 == FXP_NEG_INF) return 0;
        if (fxp1 == FXP_POS_INF) return FXP_POS_INF;
        return (fxp1 >= 0)?
                fxp_pow2_pos_arg_xfactor(fxp1, \
                                    SFXP_LG2_E_FACTOR, \
                                    FXP_POWX_LOOPS):
                fxp_pow2_neg_arg_xfactor( -fxp1, \
                                    SFXP_LG2_E_FACTOR, \
                                    FXP_POWX_LOOPS);
}

// Implementation of pow10(x) using pow2() and lg2():
// 10^x == 2^( x * lg2(10) )
int fxp_pow10(int fxp1)
{
        if (fxp1 == FXP_UNDEF) return FXP_UNDEF;
        if (fxp1 == FXP_NEG_INF) return 0;
        if (fxp1 == FXP_POS_INF) return FXP_POS_INF;
        return (fxp1 >= 0)?
                fxp_pow2_pos_arg_xfactor(fxp1, \
                                    SFXP_LG2_10_FACTOR, \
                                    FXP_POWX_LOOPS):
                fxp_pow2_neg_arg_xfactor( -fxp1, \
                                    SFXP_LG2_10_FACTOR, \
                                    FXP_POWX_LOOPS);
}

// Implementation of square root as:
// sqrt(x) = 2^( 0.5 * lg2(x) )
int fxp_sqrt(int fxp1)
{
        if (fxp1 < 0) return FXP_UNDEF;
        if (fxp1 == 0) return 0;
        if (fxp1 == FXP_POS_INF) return FXP_POS_INF;
        // First get the lg2 of the argument
        super_fxp slg = fxp_get_lg2_as_sfxp(fxp1, FXP_SQRT_LG_LOOPS);
        // Halving that value
        slg.nwbits--;
        // parameters for pow2
        int w = (int) sfxp_get_poswhole(slg);
        ulongy bkmearg = get_sfxp_frac_for_bkme(slg);
        return (slg.sign)?
                fxp_pow2_wneg(w, bkmearg, FXP_sqrt_pw_loops):
                fxp_pow2_wpos(w, bkmearg, FXP_sqrt_pw_loops);
}


// Implementation of fast inverse square root with
// the same fast approach used in Quake III
int fxp_rsqrt(int fxp1)
{
    return 0;
}


// Implementation of powxy_l(x) using pow2() and lg2():
// x^y == 2^( y * lg2(x) )
int fxp_powxy(int x, int y)
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
        super_fxp slg2x = fxp_get_lg2_as_sfxp(x, FXP_POWXY_LG_LOOPS);
        #ifdef VERBOSE
        printf("powxy:\n");
        printf("      x: "); print_fxp(x); printf("\n");
        printf("      y: "); print_fxp(y);
        print_sfxp(" lg2(x): ", slg2x);
        #endif
        // Return appropriately signed 2^( y * lg2x )
        if (slg2x.sign) {
                slg2x.sign = 0; // negating in place
                if (y >= 0) {
                        #ifdef VERBOSE
                        int pnwbits = slg2x.nwbits + fxp_nbits(y) - FXP_frac_bits;
                        printf("\ncase 2: -lgx +y, product would have %d whole bits\n", pnwbits);
                        printf("Whole bit counts:  factor: %d,  y: %d\n", \
                                        slg2x.nwbits, fxp_nbits(y) - FXP_frac_bits);
                        #endif
                        return fxp_pow2_neg_arg_xfactor(y, slg2x, FXP_POWXY_PW_LOOPS);
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
                                        fxp_pow2_pos_arg_xfactor(posy, slg2x, \
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
                                        fxp_pow2_pos_arg_xfactor(y, slg2x, \
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
                                        fxp_pow2_neg_arg_xfactor(posy, slg2x, \
                                                                FXP_POWXY_PW_LOOPS);
                }
        }
}




// Trigonometric functions

/*
*/



int fxp_tan(int fxp1)
{
        return 0;
}

int fxp_asin(int fxp1)
{
        return 0;
}

int fxp_acos(int fxp1)
{
        return 0;
}

int fxp_atan(int fxp1){
        return 0;
}
