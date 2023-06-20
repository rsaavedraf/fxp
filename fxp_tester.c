/* SPDX-License-Identifier: MIT */
/*
 * fxp_tester.c
 * Tests the implementation of binary fixed point numbers
 * (fxp.c and fxp_l.c)
 *
 * By Raul Saavedra, Bonn Germany
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <float.h>
#include <time.h>
#include <math.h>
#include "fxp.h"
#include "fxp_extern.h"
#include "fxp_l.h"
#include "fxp_aux.h"
#include "fxp_conv.h"
#include "print_as_bits.h"
#include "ulongy.h"
#define DASHES "============================================================\n"

// Set to 0 in order to be able to replicate runs
#define SET_RAND_SEED 0

#define TEST_WITH_RANDS 1
#define MAX_RAND_LOOPS 5
//#define MAX_RAND_LOOPS 5000
#define TEST_BASICS 0
#define TEST_SUPER_FXP_L 0
#define TEST_LOGARITHMS 0
#define TEST_LG2_MUL_L 0
#define TEST_POWERS 0
#define TEST_SQRT 0
#define TEST_POWXY 0
#define AVOID_EXTREME_INPUTS_FOR_POWXY 1
#define TEST_TRIGONOM 1

// Max allowed error will be WDELTA_MAX * tiniest (least signif frac bit)
#define WDELTA_MAX 3.0

// Warnings will start appearing when error >= MIN_DELTA = WDELTA_MAX*tiniest / WDELTA_DIV
#define WDELTA_DIV 1.2

static int fracbit_configs[] = {8, 11, 16, 24, 28, 31};
//static int fracbit_configs[] = {28};
/*
static int fracbit_configs[] = {
31, 30,
29, 28, 27, 26, 25, 24, 23, 22, 21, 20,
19, 18, 17, 16, 15, 14, 13, 12, 11, 10,
9, 8, 7, 6, 5, 4
};
*/

#define NINTBITS (sizeof(int) * 8)
static unsigned int fracbit_nwarnings[NINTBITS];
static unsigned int fracbit_threshold_cases[NINTBITS];
static unsigned int fracbit_extr_powxy_cases[NINTBITS];
static long double fracbit_maxdelta_allowed[NINTBITS];
static long double fracbit_maxdelta_observed[NINTBITS];

static long double min_warn_delta = 0.0;
static long double max_warn_delta = 0.0;
static long double larger_delta = 0.0;
static long double largest_delta = 0.0;
static long double largest_madelta = 0.0;
static int largest_delta_fbits = 0;
static int twarnings = 0;
static int nwarnings = 0;

static int frac_bits;
static int whole_bits;
static int frac_mask;
static int frac_max;
static int frac_max_dec;
static int whole_max;
static int whole_min;
static int fxp_largest;

static void test_fxp(char *s, long double d_assert_val, int fxp1)
{
        printf("%s\n", s);
        printf("    exp: ");
        if (d_assert_val <= FXP_UNDEF_LD)
                printf("UNDEF");
        else if (d_assert_val <= FXP_NINF_LD)
                printf("-INF");
        else if (d_assert_val >= FXP_PINF_LD)
                printf("+INF");
        else {
                printf("%17.10Le", d_assert_val);
        }
        printf("\n    act: ");
        print_fxp(fxp1);
        printf("\n");
        int b1 = ((fxp1 == FXP_UNDEF) && (d_assert_val == FXP_UNDEF_LD));
        int b2 = ((fxp1 == FXP_NEG_INF) && (d_assert_val == FXP_NINF_LD));
        int b3 = ((fxp1 == FXP_POS_INF) && (d_assert_val == FXP_PINF_LD));
        if (b1 || b2 || b3) {
                // Expected and actual are both same-signed infinities,
                // or both undef. Either way, they're the same
                //printf(" (~same)\n");
                return;
        }
        fflush(stdout);
        // Make sure it's not an UNDEF vs. NON-UNDEF case
        assert((fxp1 != FXP_UNDEF) && (d_assert_val != FXP_UNDEF_LD));
        long double df = fxp2ld(fxp1);
        long double delta;
        char * msg;
        if ((fxp1 == FXP_POS_INF) || (d_assert_val == FXP_PINF_LD)) {
                // One of them +infinity, but not the other ->
                // Calculate delta of the latter with respect to
                // the +INF threshold (FXP_max_ldx)
                delta = FXP_max_ldx - ((fxp1 == FXP_POS_INF)? d_assert_val: df);
                msg = ", rare case of extremely close values on either side of +INF threshold";
        } else if ((fxp1 == FXP_NEG_INF) || (d_assert_val == FXP_NINF_LD)) {
                // One of them -infinity, but not the other ->
                // Calculate delta of the latter with respect to
                // the -INF threshold (FXP_min_ldx)
                delta = ((fxp1 == FXP_NEG_INF)? d_assert_val: df) - FXP_min_ldx;
                msg = ", rare case of extremely close values on either side of -INF threshold";
        } else {
                // Values are either both within infinity thresholds
                // (normal case,) or both are inf's with different signs
                delta = (df >= d_assert_val)?
                                df - d_assert_val:
                                d_assert_val - df;
                msg = "";
        }

        if (delta <= min_warn_delta) {
                // No warning up to the min_warn_delta value
                if (msg[0] != '\0') {
                        printf("    (~same%s)\n", msg);
                        fracbit_threshold_cases[frac_bits]++;
                }
        } else {
                nwarnings++;
                fracbit_nwarnings[frac_bits]++;
                printf("\n***** Warning %d: d=%1.2LE for %1.2LE (from %1.2LE to MAX %1.2LE allowed for %d f.bits)\n",
                            nwarnings, delta, df,
                            min_warn_delta,
                            max_warn_delta,
                            frac_bits);
                if (delta > larger_delta) {
                        larger_delta = delta;
                        printf("***** For %d f.bits, that's the largest delta so far: %LE\n", \
                                frac_bits, larger_delta);
                        fracbit_maxdelta_observed[frac_bits] = larger_delta;
                }
                printf("\n");
        }
        fflush(stdout);

        // Here assert that we never exceed the max_warn_delta
        assert( (delta <= max_warn_delta) &&
            ((df >= 0 && d_assert_val >= 0) ||
                (df <= 0 && d_assert_val <= 0)));
}

static long double get_target(long double x)
{
        if (x <= FXP_UNDEF_LD) return FXP_UNDEF_LD;
        if (x <= FXP_min_ldx) return FXP_NINF_LD;
        if (x >= FXP_max_ldx) return FXP_PINF_LD;
        return x;
}

static long double get_lg2_target(int x)
{
        if (x < 0) return FXP_UNDEF_LD;
        if (x == 0) return FXP_NINF_LD;
        if (x == FXP_POS_INF) return FXP_PINF_LD;
        return get_target(log2l(fxp2ld(x)));
}

static long double get_ln_target(int x)
{
        if (x < 0) return FXP_UNDEF_LD;
        if (x == 0) return FXP_NINF_LD;
        if (x == FXP_POS_INF) return FXP_PINF_LD;
        return get_target(logl(fxp2ld(x)));
}

static long double get_lg10_target(int x)
{
        if (x < 0) return FXP_UNDEF_LD;
        if (x == 0) return FXP_NINF_LD;
        if (x == FXP_POS_INF) return FXP_PINF_LD;
        return get_target(log10l(fxp2ld(x)));
}

static long double get_pow2_target(int x)
{
        if (x == FXP_UNDEF) return FXP_UNDEF_LD;
        if (x == FXP_NEG_INF) return 0.0L;
        if (x == FXP_POS_INF) return FXP_PINF_LD;
        long double ldvalue = powl(2.0L, fxp2ld(x));
        long double target = get_target(ldvalue);
        return target;
}

static long double get_exp_target(int x)
{
        if (x == FXP_UNDEF) return FXP_UNDEF_LD;
        if (x == FXP_NEG_INF) return 0.0L;
        if (x == FXP_POS_INF) return FXP_PINF_LD;
        return get_target(expl(fxp2ld(x)));
}

static long double get_pow10_target(int x)
{
        if (x == FXP_UNDEF) return FXP_UNDEF_LD;
        if (x == FXP_NEG_INF) return 0.0L;
        if (x == FXP_POS_INF) return FXP_PINF_LD;
        return get_target(powl(10.0L, fxp2ld(x)));
}

static long double get_sqrt_target(int x)
{
        if ((x == FXP_UNDEF) || (x < 0)) return FXP_UNDEF_LD;
        if (x == FXP_POS_INF) return FXP_PINF_LD;
        return get_target(sqrtl(fxp2ld(x)));
}


/*
Expected output from 2^( y * lg2(x) ) given the combination of inputs x and y:

  x   y: Und. -INF [MIN,-1)  -1  (-1,0)   0   (0,1)   1  (1,MAX] +INF
< 0      Und.  Und.   Und.   Und.  Und.  Und.  Und.  Und.  Und.   Und.
  0      Und.  +INF  +INF   +INF  +INF   Und.   0     0     0      0
(0,1)    Und.   0      *      *     *     1     *     *     *      0
  1      Und.  Und.    1      1     1     1     1     1     1     Und.
(1,MAX]  Und.   0      *      *     *     1     *     *     *     +INF
 +INF    Und.   0      0      0     0    Und.  +INF  +INF  +INF   +INF
*/
static long double get_powxy_target(int x, int y)
{
        if ((x < 0) || (y == FXP_UNDEF)) {
                return FXP_UNDEF_LD;
        }
        if (x == 0) {
                return (y < 0)? FXP_PINF_LD: \
                            (y == 0)? FXP_UNDEF_LD: 0.0L;
        }
        if (x == FXP_POS_INF) {
                return (y < 0)? 0.0L: \
                            (y == 0)? FXP_UNDEF_LD: FXP_PINF_LD;
        }
        if (y == FXP_NEG_INF) {
                return (x == FXP_one)? FXP_UNDEF_LD: 0.0L;
        }
        if (y == FXP_POS_INF) {
                return (x < FXP_one)? 0.0L: \
                        (x == FXP_one)? FXP_UNDEF_LD: FXP_PINF_LD;
        }
        return ((x == FXP_one) || (y == 0))? 1.0L: \
                    get_target(powl(fxp2ld(x), fxp2ld(y)));
}


static long double get_div_target(int x, int y)
{
        long double target;
        if ((x == FXP_UNDEF) || (y == FXP_UNDEF) \
            || ((x == 0) && (y == 0)) \
            || ((x == FXP_POS_INF) && (y == FXP_POS_INF)) \
            || ((x == FXP_POS_INF) && (y == FXP_NEG_INF)) \
            || ((x == FXP_NEG_INF) && (y == FXP_POS_INF)) \
            || ((x == FXP_NEG_INF) && (y == FXP_NEG_INF))) {
            target = FXP_UNDEF_LD;
        } else {
            if (((x > 0) && (y == 0)) \
                || ((x == FXP_POS_INF) && (y > 0)) \
                || ((x == FXP_NEG_INF) && (y < 0))) {
                target = FXP_PINF_LD;
            } else {
                if (((x < 0) && (y == 0)) \
                    || ((x == FXP_NEG_INF) && (y > 0)) \
                    || ((x == FXP_POS_INF) && (y < 0))) {
                    target = FXP_NINF_LD;
                } else {
                    if ((y == FXP_POS_INF) || (y == FXP_NEG_INF)) {
                        target = 0;
                    } else {
                        long double ldx = fxp2ld(x);
                        long double ldy = fxp2ld(y);
                        target = ldx / ldy;
                        if (target >= FXP_max_ldx) {
                            target = FXP_PINF_LD;
                        } else {
                            if (target <= FXP_min_ldx) {
                                target = FXP_NINF_LD;
                            }
                        }
                    }
                }
            }
        }
        return target;
}

static long double get_mul_target(int x, int y)
{
        long double target;
        if (((x == FXP_UNDEF) || (y == FXP_UNDEF)) \
            || ((x == 0) && ((y == FXP_POS_INF) || (y == FXP_NEG_INF))) \
            || ((y == 0) && ((x == FXP_POS_INF) || (x == FXP_NEG_INF)))) {
            target = FXP_UNDEF_LD;
        } else {
            if ((x == FXP_POS_INF) || (x == FXP_NEG_INF)) {
                if (y < 0)
                    target = (x == FXP_POS_INF)? FXP_NINF_LD: FXP_PINF_LD;
                else
                    target = (x == FXP_POS_INF)? FXP_PINF_LD: FXP_NINF_LD;
            } else {
                if ((y == FXP_POS_INF) || (y == FXP_NEG_INF)) {
                    if (x < 0)
                        target = (y == FXP_POS_INF)? FXP_NINF_LD: FXP_PINF_LD;
                    else
                        target = (y == FXP_POS_INF)? FXP_PINF_LD: FXP_NINF_LD;
                } else {
                    target = fxp2ld(x) * fxp2ld(y);
                    if (target <= FXP_min_ldx)
                        target = FXP_NINF_LD;
                    else if (target >= FXP_max_ldx)
                        target = FXP_PINF_LD;
                }
            }
        }
        return target;
}

static long double get_f_target(int x)
{
        if (x == FXP_UNDEF) return FXP_UNDEF_LD;
        if (x == FXP_NEG_INF) return FXP_NINF_LD;
        if (x == FXP_POS_INF) return FXP_PINF_LD;
        float f = fxp2f(x);
        if (f <= FXP_min_fx) return FXP_NINF_LD;
        if (f >= FXP_max_fx) return FXP_PINF_LD;
        return (long double) f;
}

static long double get_sin_target(int x)
{
        if ((x == FXP_UNDEF) || (x == FXP_NEG_INF) || (x == FXP_POS_INF))
                return FXP_UNDEF_LD;
        long double ldvalue = sinl(fxp2ld(x));
        long double target = get_target(ldvalue);
        return target;
}

static long double get_cos_target(int x)
{
        if ((x == FXP_UNDEF) || (x == FXP_NEG_INF) || (x == FXP_POS_INF))
                return FXP_UNDEF_LD;
        long double ldvalue = cosl(fxp2ld(x));
        long double target = get_target(ldvalue);
        return target;
}

static long double get_tan_target(int x)
{
        if ((x == FXP_UNDEF) || (x == FXP_NEG_INF) || (x == FXP_POS_INF))
                return FXP_UNDEF_LD;
        long double ldvalue = tanl(fxp2ld(x));
        long double target = get_target(ldvalue);
        return target;
}

static long double get_asin_target(int x)
{
        if ((x > FXP_one) || (x < FXP_minus_one)) return FXP_UNDEF_LD;
        long double ldvalue = asinl(fxp2ld(x));
        long double target = get_target(ldvalue);
        return target;
}

static long double get_acos_target(int x)
{
        if ((x > FXP_one) || (x < FXP_minus_one)) return FXP_UNDEF_LD;
        long double ldvalue = acosl(fxp2ld(x));
        long double target = get_target(ldvalue);
        return target;
}

static long double get_atan_target(int x)
{
        if (x == FXP_UNDEF) return FXP_UNDEF_LD;
        long double ldvalue = atanl(fxp2ld(x));
        long double target = get_target(ldvalue);
        return target;
}

static int FXP_TINIEST = 1;
static long double dfxp_tiniest;

static int zero;
static int fxp_ten;
static int fxp_two;
static int fxp_one;
static int fxp_point5;
static int fxp_halfmax;
static int fxp_halfp2;

void tests_01()
{
        printf("\nChecking extreme int values, part I, ");
        printf("frac bits: %d\n", frac_bits);
        test_fxp("+Inf",
                    FXP_PINF_LD,
                    FXP_POS_INF);
        test_fxp("Largest",
                    fxp2ld(FXP_MAX),
                    fxp_largest);
        test_fxp("HalfMax",
                    fxp2ld(FXP_MAX)/2,
                    fxp_halfmax);
        test_fxp("0.5",
                    ((long double) 0.5),
                    fxp_point5);
        test_fxp("tiniest",
                    dfxp_tiniest,
                    FXP_TINIEST);
        test_fxp("zero",
                    0.0,
                    fxp_bin(0, 0));
        test_fxp("-tiniest",
                    fxp2ld(-FXP_TINIEST),
                    -FXP_TINIEST);
        test_fxp("Most negative",
                    fxp2ld(FXP_MIN),
                    -fxp_largest);
        test_fxp("-Inf",
                    FXP_NINF_LD,
                    FXP_NEG_INF);
        test_fxp("Undefined",
                    FXP_UNDEF_LD,
                    FXP_UNDEF);
}

void tests_02()
{
        printf("\nChecking extreme int values, part II, ");
        printf("frac bits: %d\n", frac_bits);

        test_fxp("Almost most negative",
                    fxp2ld(-fxp_largest) + dfxp_tiniest,
                    fxp_add(-fxp_largest, FXP_TINIEST));
        test_fxp(" Largest -Largest",
                    0.0,
                    fxp_add(fxp_largest, -fxp_largest));
        test_fxp("-Largest +Largest",
                    0.0,
                    fxp_add(-fxp_largest, fxp_largest));
        test_fxp("Largest + 0",
                    fxp2ld(fxp_largest),
                    fxp_add(fxp_largest, zero));
        test_fxp("-Largest - 0",
                    fxp2ld(-fxp_largest),
                    fxp_add(-fxp_largest, -zero));
        test_fxp("Largest - tiniest",
                    fxp2ld(fxp_largest) - dfxp_tiniest,
                    fxp_sub(fxp_largest, FXP_TINIEST));
        test_fxp("Largest + tiniest safe",
                    FXP_PINF_LD,
                    fxp_add(fxp_largest, FXP_TINIEST));
        test_fxp("Largest + tiniest unsafe",
                    FXP_PINF_LD,
                    fxp_unsafe_add(fxp_largest, FXP_TINIEST));
        test_fxp("-(+inf)",
                    FXP_NINF_LD,
                    -FXP_POS_INF);
        test_fxp("-(-inf)",
                    FXP_PINF_LD,
                    -FXP_NEG_INF);
        test_fxp("+inf + +inf",
                    FXP_PINF_LD,
                    fxp_add(FXP_POS_INF, FXP_POS_INF));
        test_fxp("-inf - +inf",
                    FXP_NINF_LD,
                    fxp_sub(FXP_NEG_INF, FXP_POS_INF));
        test_fxp("+inf + -inf",
                    FXP_UNDEF_LD,
                    fxp_add(FXP_POS_INF, FXP_NEG_INF));
        test_fxp("-inf + -inf",
                    FXP_NINF_LD,
                    fxp_add(FXP_NEG_INF, FXP_NEG_INF));
        test_fxp("-inf - -inf",
                    FXP_UNDEF_LD,
                    fxp_sub(FXP_NEG_INF, FXP_NEG_INF));
        test_fxp("+inf * -inf",
                    FXP_NINF_LD,
                    fxp_mul(FXP_POS_INF, FXP_NEG_INF));
        test_fxp("+inf - 0.5",
                    FXP_PINF_LD,
                    fxp_sub(FXP_POS_INF, fxp_point5));
        test_fxp("-inf + 0.5",
                    FXP_NINF_LD,
                    fxp_add(FXP_NEG_INF, fxp_point5));
        test_fxp("+num / zero",
                    FXP_PINF_LD,
                    fxp_div(fxp_largest, zero));
        test_fxp("zero / zero",
                    FXP_UNDEF_LD,
                    fxp_div(zero, zero));
        test_fxp("zero * zero",
                    0.0,
                    fxp_mul(zero, zero));
        test_fxp("zero + zero",
                    0.0,
                    fxp_add(zero, zero));
        test_fxp("zero - zero",
                    0.0,
                    fxp_sub(zero, zero));
        test_fxp("zero - undef",
                    FXP_UNDEF_LD,
                    fxp_sub(zero, FXP_UNDEF));
        test_fxp("-num / zero",
                    FXP_NINF_LD,
                    fxp_div(-fxp_largest, zero));
        test_fxp("zero * +inf",
                    FXP_UNDEF_LD,
                    fxp_mul(zero, FXP_POS_INF));
        test_fxp("zero * -inf",
                    FXP_UNDEF_LD,
                    fxp_mul(zero, FXP_NEG_INF));
        test_fxp("zero * undef",
                    FXP_UNDEF_LD,
                    fxp_mul(zero, FXP_UNDEF));
        test_fxp("-inf * undef",
                    FXP_UNDEF_LD,
                    fxp_mul(FXP_NEG_INF, FXP_UNDEF));
        test_fxp("+inf * undef",
                    FXP_UNDEF_LD,
                    fxp_mul(FXP_POS_INF, FXP_UNDEF));
        test_fxp("undef * undef",
                    FXP_UNDEF_LD,
                    fxp_mul(FXP_UNDEF, FXP_UNDEF));
        test_fxp("tiniest * inf",
                    FXP_PINF_LD,
                    fxp_mul(FXP_TINIEST, FXP_POS_INF));
}

void tests_03()
{
        printf("\nChecking extreme int values, part III, ");
        printf("frac bits: %d\n", frac_bits);
        test_fxp("Way Too Large whole part!",
                    FXP_PINF_LD,
                    fxp_add(fxp(whole_max), fxp(5)));
        test_fxp("Largest * 1",
                    fxp2ld(fxp_mul(fxp_largest, fxp(1))),
                    fxp_mul(fxp_largest,  fxp_one));
        test_fxp("Largest * -1",
                    fxp2ld(fxp_mul(fxp_largest, fxp(-1))),
                    fxp_mul(fxp_largest,  -fxp_one));
        test_fxp("Largest + two safe",
                    FXP_PINF_LD,
                    fxp_add(fxp_largest, fxp_two));
        test_fxp("Largest + two unsafe",
                    fxp2ld(fxp_largest + fxp_two),
                    fxp_unsafe_add(fxp_largest, fxp_two));
        test_fxp("Safe Too neg substraction",
                    FXP_NINF_LD,
                    fxp_sub( -fxp_largest, fxp_point5));
        test_fxp("Unsafe Too neg substraction",
                    fxp2ld(-fxp_largest - fxp_point5),
                    fxp_unsafe_sub( -fxp_largest, fxp_point5));
        test_fxp("Largest + 0.5",
                    FXP_PINF_LD,
                    fxp_add(fxp_largest, fxp_point5));
        test_fxp("-Largest - 0.5",
                    FXP_NINF_LD,
                    fxp_add(-fxp_largest, -fxp_point5));
        test_fxp("+HalfMax + HMaxp2",
                    FXP_PINF_LD,
                    fxp_add(fxp_halfmax, fxp_halfp2));
        test_fxp("-HalfMax - HMaxp2",
                    FXP_NINF_LD,
                    fxp_add(-fxp_halfmax, -fxp_halfp2));
        test_fxp("HalfMax + HalfMax",
                    fxp2ld(fxp_halfmax + fxp_halfmax),
                    fxp_add(fxp_halfmax, fxp_halfmax));
        test_fxp("FXP_MAX - HalfMax",
                    fxp2ld(FXP_MAX - fxp_halfmax),
                    fxp_sub(FXP_MAX, fxp_halfmax));
        test_fxp("HalfMax + FXP_MAX",
                    FXP_PINF_LD,
                    fxp_add(fxp_halfmax, FXP_MAX));
        test_fxp("-FXP_MAX - HalfMax",
                    FXP_NINF_LD,
                    fxp_sub(-FXP_MAX, fxp_halfmax));
        test_fxp("HalfMax * 2",
                    fxp2ld(fxp_mul(fxp_halfmax, fxp(2))),
                    fxp_mul(fxp_halfmax, fxp(2)));
        test_fxp("HalfMax * 2 (long)",
                    fxp2ld(fxp_mul_l(fxp_halfmax, fxp(2))),
                    fxp_mul_l(fxp_halfmax, fxp(2)));
        test_fxp("HalfMax * 3",
                    FXP_PINF_LD,
                    fxp_mul(fxp_halfmax, fxp(3)));
        test_fxp("-HalfMax * 3",
                    FXP_NINF_LD,
                    fxp_mul(-fxp_halfmax, fxp(3)));
        test_fxp("(HalfMax+0.5)*2",
                    FXP_PINF_LD,
                    fxp_mul(fxp_add(fxp_halfmax, fxp_point5), fxp(2)));
        test_fxp("(HalfMax+0.5)*2 (long)",
                    FXP_PINF_LD,
                    fxp_mul_l(fxp_add(fxp_halfmax, fxp_point5), fxp(2)));
}

void test_decbin_mappings()
{
        printf("\nChecking decimal <=> bin mappings of frac ranges, ");
        printf("frac bits: %d\n", frac_bits);
        int fmbin = frac_max;
        int fmdec = frac_max_dec;
        printf("Max frac dec: %d (bin %d)", fmdec, fmbin);
        for (int i = 0; i <= 5; i++) {
                printf("\nShowing fxp for 0.%7d: ", i);
                int vf = fxp_dec(0, i);
                print_fxp(vf);
        }
        printf("\n:");
        int m = (fmdec + 1) / 2;
        for (int i = (m >= 2 ? m - 2: 0); i <= m + 2; i++) {
                printf("\nShowing fxp for 0.%7d: ", i);
                int vf = fxp_dec(0, i);
                print_fxp(vf);
        }
        printf("\n:");
        for (int i = (fmdec >= 5? fmdec - 5: 0);
                    i <= fmdec; i++) {
                printf("\nShowing fxp for 0.%7d: ", i);
                int vf = fxp_dec(0, i);
                print_fxp(vf);
        }
        printf("\n");
}

void test_fracs()
{
        printf("\nChecking sign taken from frac when whole == 0, ");
        printf("frac bits: %d\n", frac_bits);
        // Temporary switching to a fixed frac_max_dec for these tests,
        // and relaxing the max allowed delta
        long double bkp_minwd = min_warn_delta;
        long double bkp_ld = larger_delta;
        long double bkp_maxwd = max_warn_delta;
        int fmdec = fxp_get_frac_max_dec();

        fxp_set_frac_max_dec(999);
        min_warn_delta = fxp2ld(fxp_dec(0,100));
        max_warn_delta = min_warn_delta * 2;

        int a = fxp_dec(-0, 500);
        int b = fxp_dec(-0, -500);
        printf("-0. 500: "); print_fxp(a); printf("\n");
        printf("-0.-500: "); print_fxp(b); printf("\n");
        test_fxp("-0.(+)500",  0.5, fxp_dec(-0, 500));
        test_fxp("-0.(-)500", -0.5, fxp_dec(-0, -500));

        printf("\nChecking truncation of longer frac decimal arguments, ");
        printf("frac bits: %d\n", frac_bits);
        test_fxp("0.22222",  0.222, fxp_dec(0, 22222));
        test_fxp("0.444444", 0.444, fxp_dec(0, 444444));
        test_fxp("0.771999", 0.771, fxp_dec(0, 771999));
        test_fxp("0.999999", 0.999, fxp_dec(0, 999999));
        // Restoring original frac_max_dec and delta vars
        fxp_set_frac_max_dec(fmdec);
        min_warn_delta = bkp_minwd;
        larger_delta = bkp_ld;
        max_warn_delta = bkp_maxwd;

        printf("\nChecking extreme frac values, ");
        printf("frac bits: %d\n", frac_bits);
        test_fxp("Largest frac",
                    fxp2ld(frac_max),
                    frac_max);
        test_fxp("-Largest frac",
                    fxp2ld(-frac_max),
                    -frac_max);
        test_fxp("0.5 + fracmax",
                    ((whole_bits > 1)?
                        0.5 + fxp2ld(frac_max):
                        fxp2ld(FXP_POS_INF)),
                    fxp_add(fxp_point5, frac_max));
        test_fxp("fracmax +tiniest",
                    fxp2ld(fxp_one),
                    fxp_add(frac_max, FXP_TINIEST));
        test_fxp("-fracmax -tiniest",
                    fxp2ld(-fxp_one),
                    fxp_add(-frac_max, -FXP_TINIEST));
        test_fxp("fracmax - fracmax",
                    0.0,
                    fxp_add(frac_max, -frac_max));
}

void test_ops_with_whole_bits()
{
        if (whole_bits < 3) return;
        printf("\nChecking simple ops when using 3+ bits for whole part:\n");
        int whole = 0;
        int bin_frac = 1u << (frac_bits - 1);  // == 0.500
        int num = fxp_bin(whole, bin_frac);
        int dec_frac = fxp_get_frac_part_dec(num);
        int fxp1 = num;
        int fxp2 = fxp(2);
        long double dnum = fxp2ld(num);
        test_fxp(" 1 + 1",
                    fxp2ld(fxp2),
                    fxp_add(fxp_one, fxp_one));
        test_fxp("-1 - 1",
                    fxp2ld(-fxp2),
                    fxp_add(-fxp_one, -fxp_one));
        test_fxp("Ok sum == 2",
                    fxp2ld(fxp_add(-fxp_halfmax, fxp_halfp2)),
                    fxp_add(-fxp_halfmax, fxp_halfp2));
        test_fxp("Ok sum == -2",
                    fxp2ld(fxp_add(fxp_halfmax, -fxp_halfp2)),
                    fxp_add(fxp_halfmax, -fxp_halfp2));
        test_fxp(" num", dnum, fxp1);
        test_fxp(" num +  2",
                    dnum + 2.0,
                    fxp_add(fxp1, fxp2));
        test_fxp(" num + -2",
                    dnum - 2.0,
                    fxp_add(fxp1, -fxp2));
        test_fxp("-num +  2",
                    -dnum + 2.0,
                    fxp_add(-fxp1, fxp2));
        test_fxp("-num + -2",
                    -dnum - 2.0,
                    fxp_add(-fxp1, -fxp2));
        test_fxp(" num -  2",
                    dnum - 2.0,
                    fxp_sub(fxp1, fxp2));
        test_fxp(" num - -2",
                    dnum + 2.0,
                    fxp_sub(fxp1, -fxp2));
        test_fxp("-num -  2",
                    -dnum - 2.0,
                    fxp_sub(-fxp1, fxp2));
        test_fxp("-num - -2",
                    -dnum + 2.0,
                    fxp_sub(-fxp1, -fxp2));
        test_fxp(" num *  2",
                    dnum * 2.0,
                    fxp_mul(fxp1, fxp2));
        test_fxp(" num * -2",
                    dnum * -2.0,
                    fxp_mul(fxp1, -fxp2));
        test_fxp("-num *  2",
                    -dnum * 2.0,
                    fxp_mul(-fxp1, fxp2));
        test_fxp("-num * -2",
                    -dnum * -2.0,
                    fxp_mul(-fxp1, -fxp2));
        test_fxp(" num *  2 (long)",
                    dnum * 2.0,
                    fxp_mul_l( fxp1, fxp2));
        test_fxp(" num * -2 (long)",
                    dnum * -2.0,
                    fxp_mul_l( fxp1, -fxp2));
        test_fxp("-num *  2 (long)",
                    -dnum * 2.0,
                    fxp_mul_l(-fxp1, fxp2));
        test_fxp("-num * -2 (long)",
                    -dnum * -2.0,
                    fxp_mul_l(-fxp1, -fxp2));
        test_fxp(" num /  2",
                    dnum / 2.0,
                    fxp_div( fxp1, fxp2));
        test_fxp(" num / -2",
                    dnum / -2.0,
                    fxp_div(fxp1, -fxp2));
        test_fxp("-num /  2",
                    -dnum / 2.0,
                    fxp_div(-fxp1, fxp2));
        test_fxp("-num / -2",
                    -dnum / -2.0,
                    fxp_div(-fxp1, -fxp2));
}

void test_ops_with_values_of_interest()
{
        int sign1, sign2, sign3, n1, n2, n3, n4, fxp1, fxp2;
        long double ldx, ldy, ldz, tgt1, tgt2;

        //printf("\nVerifying multiplication with values of interest:\n");
        printf("\nVerifying ops with values of interest, ");
        printf("frac bits: %d\n", frac_bits);
        int ax[] = {
                    FXP_UNDEF, FXP_NEG_INF, FXP_POS_INF, FXP_MAX,
                    fxp_dec(33333, 333), fxp_dec(33, 33),
                    fxp(2),
                    fxp(1),
                    fxp_bin(0, FXP_frac_max),
                    1, 0, fxp_bin(0, -1),
                    fxp_bin(0, -FXP_frac_max),
                    fxp(-1),
                    fxp(-2),
                    fxp_dec(-33, 33), fxp_dec(-33333, 333),
                    -1046690937,
                    FXP_MIN
                    };
        int ay[] = {
                    FXP_UNDEF, FXP_NEG_INF, FXP_POS_INF, FXP_MAX,
                    fxp(2), fxp(1),
                    fxp_bin(0, FXP_frac_max),
                    2, 1, 0, -1, -2,
                    fxp(-1), fxp(-2),
                    -1094861345,
                    FXP_MIN
                    };
        int x, y, posx, posy;
        int ndd = (int) (sizeof(ax) / sizeof(int));
        int ndr = (int) (sizeof(ay) / sizeof(int));
        for (int i = 0; i < ndd; i++) {
                x = ax[i];
                fflush(stdout);
                // Test that round trip conversions fxp -> double or
                // long double -> fxp result in exactly the same fxp
                assert( d2fxp(fxp2d(x)) == x );
                assert( ld2fxp(fxp2ld(x)) == x );

                // With float conversions we cannot expect perfect
                // accurary, but we can check that the
                // conversion keeps values very close
                if ((i > 2) && ((x == FXP_POS_INF) \
                        || (x == FXP_NEG_INF))) {
                        continue;
                }
                float fx = fxp2f(x);
                int x2 = f2fxp(fx);
                printf("\nRoundtrip fxp->float->fxp conversion of x:\n");
                test_fxp("x: ", get_f_target(x), x2);

                for (int j = 0; j<ndr; j++) {
                        y = ay[j];
                        if ((j > 2) && ((y == FXP_POS_INF) \
                                || (y == FXP_NEG_INF))) {
                                continue;
                        }
                        printf("\nx      : "); print_fxp(x); printf("\n");
                        printf("y      : "); print_fxp(y); printf("\n");
                        ldx = fxp2ld(x);
                        ldy = fxp2ld(y);

                        //For multiplication
                        tgt1 = get_mul_target(x, y);
                        //n1 = fxp_mul_l(x, y);
                        //test_fxp("\nmul_l  (x*y)", tgt1, n1);

                        n1 = fxp_mul(x, y);
                        test_fxp("\nmul  (x*y)", tgt1, n1);

                        //For division
                        tgt1 = get_div_target(x, y);
                        n1 = fxp_div(x, y);
                        test_fxp("div  (x/y)", tgt1, n1);
                }
        }
}

void test_super_fxp_l()
{
        printf("\nTesting super_fxp_l\n");
        int pos = fxp_bin(5, 5);
        super_fxp_l spos = sfxp_l_from_fxp(pos);
        int vpos = sfxp_l_2_fxp(spos);
        test_fxp("5.5<->sfxp_l", fxp2ld(pos), vpos);
        int v0 = sfxp_l_2_fxp(sfxp_l_from_fxp(0));
        test_fxp("0<->sfxp_l", 0.0L, v0);
        int neg = -pos;
        super_fxp_l sneg = sfxp_l_from_fxp(neg);
        int vneg = sfxp_l_2_fxp(sneg);
        test_fxp("-5.5<->sfxp_l", fxp2ld(neg), vneg);
}

void test_lg2(char * msg, int x)
{
        long double tgt = get_lg2_target(x);
        printf("lg2_l("); test_fxp(msg, tgt, fxp_lg2_l(x));
        printf("lg2(");   test_fxp(msg, tgt, fxp_lg2(x));
}

void test_lg2_mul_l(char * msg, int x)
{
        printf("lg2_mul_l(");
        test_fxp(msg, get_lg2_target(x), fxp_lg2_mul_l(x));
}

void test_ln(char * msg, int x)
{
        long double tgt = get_ln_target(x);
        printf("ln_l("); test_fxp(msg, tgt, fxp_ln_l(x));
        printf("ln(");   test_fxp(msg, tgt, fxp_ln(x));
}

void test_lg10(char * msg, int x)
{
        long double tgt = get_lg10_target(x);
        printf("lg10_l("); test_fxp(msg, tgt, fxp_lg10_l(x));
        printf("lg10(");   test_fxp(msg, tgt, fxp_lg10(x));
}

void test_log(char base, char * msg, int x)
{
        switch(base) {
        case '2':
                test_lg2(msg, x);
                break;
        case 'M':
                test_lg2_mul_l(msg, x);
                break;
        case 'A':
                test_lg10(msg, x);
                break;
        default:    // 'e'
                test_ln(msg, x);
                break;
        }
}

void test_logarithms()
{
        printf("\nShowing Transcendental constants as fxp's: ");
        printf("frac bits: %d\n", frac_bits);
        printf("e      : "); print_fxp(fxp_get_e()); printf("\n");
        printf("pi     : "); print_fxp(fxp_get_pi()); printf("\n");
        printf("ln(2)  : "); print_fxp(fxp_get_ln_2()); printf("\n");
        printf("lg10(2): "); print_fxp(fxp_get_lg10_2()); printf("\n");

        printf("\nTesting logarithms for %d frac bits:", frac_bits);
        // M is also base 2 but using the multiplication algorithm
        char pbases[] = {'2', 'M', 'e', 'A'};
        int npbases = sizeof(pbases) / sizeof(pbases[0]);
        for (int i = 0; i < npbases; i++) {
                char base = pbases[i];
                if (base == 'M') {
                        if (!TEST_LG2_MUL_L) continue;
                        if (whole_bits < 3) {
                                printf("\nOnly %d whole bit(s); 3 or more needed for lg2_mul_l; skipping.", \
                                            whole_bits);
                                 continue;
                        }
                }
                //printf("\nBase %c:\n", base);
                printf("\n");
                test_log(base, "+INF)",       FXP_POS_INF);
                test_log(base, "largest)",    FXP_MAX);
                test_log(base, "100)",        fxp(100));
                test_log(base, "e+1)",        fxp_add(fxp_get_e(), fxp(1)));
                test_log(base, "e)",          fxp_get_e());
                test_log(base, "2.2)",        d2fxp(2.2));
                test_log(base, "2)",          fxp(2));
                test_log(base, "1.99...9)",   fxp_sub(fxp(2), 1));
                test_log(base, "1.5)",        d2fxp(1.5));
                test_log(base, "1.00..01)",   fxp_add(fxp(1), 1));
                test_log(base, "1)",          fxp(1));
                test_log(base, "0.99...9)",   fxp_sub(fxp(1), 1));
                test_log(base, "0.9)",        d2fxp(0.9));
                test_log(base, "0.50..01)",   fxp_add(FXP_half, 1));
                test_log(base, "0.5)",        FXP_half);
                test_log(base, "2*tiniest)",  2);
                test_log(base, "tiniest)",    1);
                test_log(base, "0)",          0);
                test_log(base, "-INF)",       FXP_NEG_INF);
                test_log(base, "UNDEF)",      FXP_UNDEF);
                if (base == '2') {
                        int exponent = (1u << (whole_bits - 1));
                        int k = d2fxp(powl(2.0, -exponent ));
                        if ((k > 0) && (exponent < FXP_INT_BITS)) {
                                printf("\nWith only %d whole bit(s) (%d frac bits,)",
                                            whole_bits, frac_bits);
                                printf(" range for whole part is +/-%d,\n", FXP_whole_max);
                                printf("lg2(x) must necessarily overflow returning -INF for ");
                                printf("all x, 0 < x <= k = \n         ");
                                print_fxp(k); printf("\n");
                                test_lg2("k + tiniest)", fxp_add(k, 1));
                                test_lg2("k)", k);
                                test_lg2("k - tiniest)", fxp_sub(k, 1));
                        }
                }
        }
}

void test_pow2(char * msg, int x)
{
        long double tgt = get_pow2_target(x);
        printf("pow2_l("); test_fxp(msg, tgt, fxp_pow2_l(x));
        printf("pow2("); test_fxp(msg, tgt, fxp_pow2(x));
}

void test_exp(char * msg, int x)
{
        long double tgt = get_exp_target(x);
        printf("exp_l("); test_fxp(msg, tgt, fxp_exp_l(x));
        printf("exp("); test_fxp(msg, tgt, fxp_exp(x));
}

void test_pow10(char * msg, int x)
{
        long double tgt = get_pow10_target(x);
        printf("pow10_l("); test_fxp(msg, tgt, fxp_pow10_l(x));
        printf("pow10("); test_fxp(msg, tgt, fxp_pow10(x));
}

void test_pow(char base, char * msg, int x)
{
        switch(base) {
        case '2':
                test_pow2(msg, x);
                break;
        case 'A':
                test_pow10(msg, x);
                break;
        default:    // 'e'
                test_exp(msg, x);
                break;
        }
}

void test_sqrt(char * msg, int x)
{
        long double tgt = get_sqrt_target(x);
        printf("sqrt_l("); test_fxp(msg, tgt, fxp_sqrt_l(x));
        printf("sqrt(");   test_fxp(msg, tgt, fxp_sqrt(x));
}

void test_powxy(char * msg, int x, int y)
{
        long double tgt = get_powxy_target(x, y);
        printf("powxy_l("); test_fxp(msg, tgt, fxp_powxy_l(x, y));
        printf("powxy(");   test_fxp(msg, tgt, fxp_powxy(x, y));
}

void test_powers()
{
        printf("\nTesting Powers of 2, e, and 10 for %d frac bits:\n", frac_bits);

        char pbases[] = {'2', 'e', 'A'};
        int npbases = sizeof(pbases) / sizeof(pbases[0]);
        for (int i = 0; i < npbases; i++) {
                char base = pbases[i];
                printf("\n");
                // Back when still not using ulongys in fxp.c, for
                // certain number of frac bits these "kcrit" values
                // used to elicit a relatively large frac error when
                // using them as exponent for the exp or pow10
                // functions, so then a more permissive max delta
                // was required. Once using ulongys, these notion of
                // kcrit arguments no longer applies whatsoever,
                // since exp and pow10 are exactly as accurate as
                // exp_l and pow10_l. But leaving these tests in
                // case of checking other fixed point libraries or
                // implementations later on.
                int kcrit;
                switch (base) {
                case '2':
                        kcrit = (int) truncl(log2l(whole_max + 0.1));
                        printf("kcrit argument for pow2() using %d frac bits: %d\n", \
                                frac_bits, kcrit);
                        break;
                case 'A':
                        kcrit = (int) truncl(log10l(whole_max + 0.1));
                        printf("kcrit argument for pow10() using %d frac bits: %d\n", \
                                frac_bits, kcrit);
                        break;
                default:
                        kcrit = (int) truncl(logl(whole_max + 0.1));
                        printf("kcrit argument for exp() using %d frac bits: %d\n", \
                                frac_bits, kcrit);
                        break;
                }
                kcrit = fxp(kcrit);
                test_pow(base, "kcrit+tiniest)", fxp_add(kcrit, 1));
                test_pow(base, "kcrit)",         kcrit);
                test_pow(base, "kcrit-tiniest)", fxp_sub(kcrit, 1));
                test_pow(base, "largest)",    FXP_MAX);
                test_pow(base, "warn29fb)",    316801043);
                test_pow(base, "3.5)",        d2fxp(3.5));
                test_pow(base, "2)",          fxp(2));
                test_pow(base, "2-tiniest)",  fxp_sub(fxp(2), 1));
                test_pow(base, "1.5)",        d2fxp(1.5));
                test_pow(base, "1+tiniest)",  fxp_add(fxp(1), 1));
                test_pow(base, "1)",          fxp(1));
                test_pow(base, "1-tiniest)",  fxp_sub(fxp(1), 1));
                test_pow(base, "1-2*tiniest)", fxp_sub(fxp(1), 2));
                test_pow(base, "0.5)",        d2fxp(0.5));
                test_pow(base, "4*tiniest)",  4);
                test_pow(base, "2*tiniest)",  2);
                test_pow(base, "tiniest)",    1);
                test_pow(base, "0)",          0);
                test_pow(base, "-tiniest)",   -1);
                test_pow(base, "-2*tiniest)", -2);
                test_pow(base, "-3*tiniest)", -3);
                test_pow(base, "-0.5)",       d2fxp(-0.5));
                test_pow(base, "-1+tiniest)", fxp_add(fxp(-1), 1));
                test_pow(base, "-1)",         fxp(-1));
                test_pow(base, "-1-tiniest)", fxp_sub(fxp(-1), 1));
                test_pow(base, "-1.5)",       d2fxp(-1.5));
                test_pow(base, "-2)",         fxp(-2));
                test_pow(base, "-7.9)",       d2fxp(-7.9));
                test_pow(base, "-8)",         fxp(-8));
                test_pow(base, "-8.1)",       d2fxp(-8.1));
                test_pow(base, "-10.5)",      d2fxp(-10.5));
                test_pow(base, "-16.5)",      d2fxp(-16.5));
                test_pow(base, "-30)",        fxp(-30));
                test_pow(base, "-31)",        fxp(-31));
                test_pow(base, "-32)",        fxp(-32));
                test_pow(base, "-42)",        fxp(-42));
                test_pow(base, "most neg)",   -fxp_largest);
                test_pow(base, "-INF)",       FXP_NEG_INF);
        }
}

void test_sqroots()
{
        printf("\nTesting square roots for %d frac bits:\n", FXP_frac_bits);
        test_sqrt("largest)",   FXP_MAX);
        test_sqrt("16)",        fxp(16));
        test_sqrt("13.7)",      d2fxp(13.7));
        test_sqrt("12.1)",      d2fxp(12.1));
        test_sqrt("4)",         fxp(4));
        test_sqrt("2)",         fxp(2));
        test_sqrt("1.5)",       d2fxp(1.5));
        test_sqrt("1+tiniest)", fxp_add(fxp(1), 1));
        test_sqrt("1)",         fxp(1));
        test_sqrt("1-tiniest)", fxp_sub(fxp(1), 1));
        test_sqrt("0.25)",      d2fxp(0.25));
        test_sqrt("0.0625)",    d2fxp(0.0625));
        test_sqrt("0.5)",       d2fxp(0.5));
        test_sqrt("0.125)",     d2fxp(0.125));
        test_sqrt("0.333)",     d2fxp(0.33333));
        test_sqrt("0.07)",      d2fxp(0.07));
        test_sqrt("0.75)",      d2fxp(0.75));
        test_sqrt("0.00364)",   d2fxp(0.00364));
        test_sqrt("tiniest)",   1);
}

void test_powersxy()
{
        printf("\nTesting powxy (x^y) for %d frac bits:\n", FXP_frac_bits);

        test_powxy("0^1)", 0, fxp(1));
        test_powxy("1^0)", fxp(1), 0);
        test_powxy("1^1)", fxp(1), fxp(1));
        test_powxy("1^2)", fxp(1), fxp(2));
        test_powxy("1^-2", fxp(1), fxp(-2));
        test_powxy("3^3)", fxp(3), fxp(3));
        int num1 = d2fxp(3.2);
        int num2 = d2fxp(1.5);
        test_powxy("3.2^1.5)", num1, num2);
        test_powxy("1.5^3.2)", num2, num1);
        test_powxy("3.2^-1.5)", num1, -num2);
        test_powxy("1.5^-3.2)", num2, -num1);
        test_powxy("0.5^4)", d2fxp(0.5), fxp(4));
        test_powxy("0.5^-4)", d2fxp(0.5), fxp(-4));
        test_powxy("+INF^-3)", FXP_POS_INF, fxp(-3));
        test_powxy("largest^1+tiniest)", FXP_MAX, fxp_add(FXP_one, 1));
        test_powxy("largest^1)", FXP_MAX, FXP_one);
        test_powxy("largest^almost1)", FXP_MAX, FXP_almost1);
        test_powxy("largest^tiniest)", FXP_MAX, 1);
        test_powxy("largest^-tiniest)", FXP_MAX, -1);
        test_powxy("almost1^almost1)", FXP_almost1, FXP_almost1);
        test_powxy("almost1^largest)", FXP_almost1, FXP_MAX);
        test_powxy("almost1^tiniest)", FXP_almost1, 1);
        test_powxy("tiniest^almost1)", 1, FXP_almost1);
        test_powxy("tiniest^-almost1)", 1, FXP_almost1);
        test_powxy("almost1^-tiniest)", FXP_almost1, -1);
        test_powxy("almost1^-almost1)", FXP_almost1, -FXP_almost1);
        test_powxy("tiniest^-largest)", 1, FXP_MIN);
        test_powxy("tiniest^4)", 1, fxp(4));
        test_powxy("tiniest^-4)", 1, fxp(-4));
        // Indeterminate powers
        test_powxy("0^0)", 0, 0);
        test_powxy("1^+INF)", fxp(1), FXP_POS_INF);
        test_powxy("1^-INF)", fxp(1), FXP_NEG_INF);
        test_powxy("+INF^0)", FXP_POS_INF, 0);

        /*
        Notice in the binary expansions included here below, the most error-prone
        combinations of arguments for powxy follow a distinct pattern:
             1) One value is very close to but slightly larger than 1, and
             2) The other value is quite a large number.

        Hence, the most critical combination of x and y for powxy ought to be something
        like x = 1 + tiniest, and y close to or == Largest. And that should be the case for
        all frac bit configurations.

        Such cases are allowed in the generation of random inputs for the powxy tests
        only if AVOID_EXTREME_INPUTS_FOR_POWXY is 0. But such combinations might not
        even be likely arguments for the powxy functions within a given application
        domain. Or for such extreme combinations, ultimate precision just to cope
        with them might be impractically expensive, while not of general interest.

        In any case, increasing the overall error tolerance of this tester (WDELTA_MAX), or
        increasing the number of loops used in the powxy functions, are always available
        tuning options to make the tester pass through such cases if desired.
        */

        // Formerly problematic cases (when using just tuple instead of tuple_l in fxp_l.c,
        // and having just uint instead of ulongy in fxp.c's tuple).
        // Only positive values are used for the powxy tests, so ignore the signs from these:
        // For 31 frac bits even with longs
        // n1 = -9.9673205148e-01 (fxp: x(-)7F94EA76 = -2140465782 =b(-).1111111100101001110101001110110)
        // n2 =  7.9005645076e-01 (fxp:    x652091DD =  1696633309 =b   .1100101001000001001000111011101)
        test_powxy("31fb n1,n2)", 2140465782, 1696633309);

        // For 30 frac bits even with longs, requires WDELTA_MAX = 2.8
        // n1 = 1.4423908889e+00 (fxp: x5C5021E0 =1548755424 =b1.011100010100000010000111100000)
        // n2 = 1.7836997155e+00 (fxp: x722822DA =1915232986 =b1.110010001010000010001011011010)
        // n2 =  4.1533096321e-01 (fxp:   445958226,  x1A94C852,  b00.011010100101001100100001010010)
        // n3 = -6.0732139554e-01 (fxp:  -652106383, -x26DE5A8F, -b00.100110110111100101101010001111)
        test_powxy("30fb1 n1,n2)", 1548755424, 1915232986);
        test_powxy("30fb2 n2,n3)", 445958226, 652106383);

        // For 29 frac bits
        // n1 = 1.4221374150e+00 (fxp: x2D822653 =763504211 =b1.01101100000100010011001010011)
        // n2 = 3.9325484559e+00 (fxp: x7DD76FDC =2111270876 =b11.11101110101110110111111011100)
        test_powxy("29fb n1,n2)", 763504211, 2111270876);

        // For 26 frac bits:
        // n1 =  1.1161174327e+00 (fxp:    x 476E77D =    74901373 =b   1.00011101101110011101111101)
        // n2 =  2.8956378385e+01 (fxp:    x73D354DB =  1943229659 =b   11100.11110100110101010011011011)
        // n1 =  1.1112943739e+00 (fxp:    x 471F727 =    74577703 =b   1.00011100011111011100100111)
        // n2 =  3.1679264620e+01 (fxp:    x7EB79125 =  2125959461 =b   11111.10101101111001000100100101)
        // n1 =  1.8891474858e+01 (fxp:    x4B90DEC9 =  1267785417 =b   10010.11100100001101111011001001)
        // n2 = -1.1775795519e+00 (fxp: x(-) 4B5D76A =   -79026026 =b(-)1.00101101011101011101101010)
        test_powxy("26fb1 n1,n2)", 74901373, 1943229659);
        test_powxy("26fb2 n1,n2)", 74577703, 2125959461);
        test_powxy("26fb3 n1,n2)", 1267785417, 79026026);

        // For 25 frac bits:
        // n1 = 1.1232554913e+00 (fxp: x23F1B58 =37690200 =b1.0001111110001101101011000)
        // n2 = 3.5341912180e+01 (fxp: x46AF0F1D =1185877789 =b100011.0101011110000111100011101)
        test_powxy("25fb n1,n2)", 37690200, 1185877789);

        // For 24 frac bits:
        // The following quite troublesome case for 24 frac bits requires either an
        // error tolerance of WDELTA_MAX >= 5.1, or else, if keeping WDELTA_MAX at let's say 3.0,
        // then more precision loops are needed: FXP_POWXY_LG_LOOPS in fxp.c at least
        // FXP_INT_BITS + 5. Otherwise an assert gets triggered:
        // n1 = -1.0454773903e+00 (fxp:   -17540200, -x010BA468, -b00000001.000010111010010001101000)
        // n2 =  1.0314526862e+02 (fxp:  1730490451,  x67253053,  b01100111.001001010011000001010011)
        // test_powxy("24fb3 n1,n2)", 17540200, 1730490451);
        //
        // Here some still tough cases for 24 frac bits, but they do pass under the defaults of
        // WDELTA_MAX = 3.0, and FXP_POWXY_LG_LOOPS = FXP_INT_BITS + 3:
        // n1 = -1.0481667519e+00 (fxp: x(-) 10C54A8 =   -17585320 =b(-)1.000011000101010010101000)
        // n2 = -9.0499958158e+01 (fxp: x(-)5A7FFD42 = -1518337346 =b(-)1011010.011111111111110101000010)
        // n1 =  1.0613912344e+00 (fxp:    17807190,  x010FB756,  b00000001.000011111011011101010110)
        // n2 = -7.4797068596e+01 (fxp: -1254886576, -x4ACC0CB0, -b01001010.110011000000110010110000)
        test_powxy("24fb1 n1,n2)", 17585320, 1518337346);
        test_powxy("24fb2 n1,n2)", 17807190, 1254886576);

        // For 23 frac bits:
        // n1 = 1.0374912024e+00 (fxp: x84CC83 =8703107 =b1.00001001100110010000011)
        // n2 = 1.2279908419e+02 (fxp: x3D664864 =1030113380 =b1111010.11001100100100001100100)
        // n1 =  1.0255248547e+00 (fxp:    x  834466 =     8602726 =b   1.00000110100010001100110)
        // n2 = -2.1888473272e+02 (fxp: x(-)6D713EEC = -1836138220 =b(-)11011010.11100010011111011101100)
        test_powxy("23fb1 n1,n2)", 8703107, 1030113380);
        test_powxy("23fb2 n1,n2)", 8602726, 1836138220);
        // Other things being equal, the following quite troublesome case for 23 frac bits
        // requires a WDELTA_MAX >= 6.4
        //n1 =  1.0186194181e+00 (fxp:     8544799,  x0082621F,  b000000001.00000100110001000011111)
        //n2 =  2.4920319748e+02 (fxp:  2090467936,  x7C9A0260,  b011111001.00110100000001001100000)
        // test_powxy("23fb3 n1,n2)", 8544799, 2090467936);
}

void test_sin(char * msg, int x)
{
        long double tgt = get_sin_target(x);
        printf("sin_l("); test_fxp(msg, tgt, fxp_sin_l(x));
        printf("sin("); test_fxp(msg, tgt, fxp_sin(x));
}

void test_cos(char * msg, int x)
{
        long double tgt = get_cos_target(x);
        printf("cos_l("); test_fxp(msg, tgt, fxp_cos_l(x));
        printf("cos("); test_fxp(msg, tgt, fxp_cos(x));
}

void test_tan(char * msg, int x)
{
        long double tgt = get_tan_target(x);
        printf("tan_l("); test_fxp(msg, tgt, fxp_tan_l(x));
        printf("tan("); test_fxp(msg, tgt, fxp_tan(x));
}


void test_asin(char * msg, int x)
{
        long double tgt = get_asin_target(x);
        printf("asin_l("); test_fxp(msg, tgt, fxp_asin_l(x));
        printf("asin("); test_fxp(msg, tgt, fxp_asin(x));
}


void test_acos(char * msg, int x)
{
        long double tgt = get_acos_target(x);
        printf("acos_l("); test_fxp(msg, tgt, fxp_acos_l(x));
        printf("acos("); test_fxp(msg, tgt, fxp_acos(x));
}


void test_atan(char * msg, int x)
{
        long double tgt = get_atan_target(x);
        printf("atan_l("); test_fxp(msg, tgt, fxp_atan_l(x));
        printf("atan("); test_fxp(msg, tgt, fxp_atan(x));
}


void test_trigonometrics()
{
        printf("\nTesting Trigonometric functions for %d frac bits:\n", FXP_frac_bits);
        int angle = 0;
        test_sin("0)", angle);
        test_cos("0)", angle);
        test_tan("0)", angle);
        int num = FXP_half;
        test_asin("0.5)", num);
        test_acos("0.5)", num);
        test_atan("0.5)", num);

}

void test_ops_with_rand_nums()
{
        int sign1, sign2, sign3, n1, n2, n3, n4, fxp1, fxp2;
        int vn1, vn2, vn3;
        long double ldx, ldy, ldz, tgt1, tgt2, tgt3, pldx, pldy;
        const long double QUASI_ONE = 1.25;
        const long double QUITE_BIG = fxp2ld(FXP_MAX) * 0.25;

        for (int i = 0; i < MAX_RAND_LOOPS; i++) {
                sign1 = rand() % 2 == 1? -1: 1;
                sign2 = rand() % 2 == 1? -1: 1;
                sign3 = rand() % 2 == 1? -1: 1;
                n1 = sign1 * rand();
                n2 = sign2 * rand();
                n3 = sign3 * rand();
                if (frac_bits < FXP_INT_BITS_M1) {
                        // n3 always in (-1, 1)
                        n3 %= frac_max;
                }

                ldx = fxp2ld(n1);
                ldy = fxp2ld(n2);
                ldz = fxp2ld(n3);
                pldx = (ldx >= 0.0)? ldx: -ldx;
                pldy = (ldy >= 0.0)? ldy: -ldy;
                if (((pldx > 1.0) && (pldx <= QUASI_ONE) && (pldy >= QUITE_BIG)) \
                        ||
                     ((pldy > 1.0) && (pldy <= QUASI_ONE) && (pldx >= QUITE_BIG))) {
                        // n1 and n2 are an extreme and particularly error-prone
                        // combination of inputs for powxy:
                        // n1 greater but quite close to 1, and
                        // n2 a rather large number close enough to Largest
                        // (or the other way around.)
                        // Count the number of these cases encountered
                        printf("\nExtreme combination of inputs for powxy encountered.\n");
                        fracbit_extr_powxy_cases[FXP_frac_bits]++;
                        if (AVOID_EXTREME_INPUTS_FOR_POWXY) {
                                // Replace n1 with their average, and n2 with the
                                // average between that new n1 and n2
                                ldx = (ldx + ldy) / 2.0;
                                n1 = ld2fxp(ldx);
                                ldy = (ldx + ldy) / 2.0;
                                n2 = ld2fxp(ldy);
                                ldx = fxp2ld(n1);
                                ldy = fxp2ld(n2);
                        }
                }
                printf("\nRandom Test #%d, ", i);
                printf("frac bits: %d", frac_bits);
                printf("\n    n1 = "); print_fxp(n1);
                printf("\n    n2 = "); print_fxp(n2);
                printf("\n    n3 = "); print_fxp(n3); printf("\n");

                if (TEST_BASICS) {
                        // Test that roundtrip conversions from fxp
                        // to double or long double, then back to fxp
                        // results in exactly the same original fxp
                        fflush( stdout );
                        assert( d2fxp(fxp2d(n1)) == n1 );
                        assert( d2fxp(fxp2d(n2)) == n2 );
                        assert( d2fxp(fxp2d(n3)) == n3 );

                        assert( ld2fxp(fxp2ld(n1)) == n1 );
                        assert( ld2fxp(fxp2ld(n2)) == n2 );
                        assert( ld2fxp(fxp2ld(n3)) == n3 );

                        // Test that roundtrip conversions from fxp
                        // to float and back to fxp remains very close
                        // to the original
                        test_fxp("n1 from float roundtrip: ", get_f_target(n1), f2fxp(fxp2f(n1)));
                        test_fxp("n2 from float roundtrip: ", get_f_target(n2), f2fxp(fxp2f(n2)));
                        test_fxp("n3 from float roundtrip: ", get_f_target(n3), f2fxp(fxp2f(n3)));

                        fxp1 = fxp_add(n1, n2);
                        tgt1 = get_target(ldx + ldy);
                        test_fxp("add   (n1+n2)", tgt1, fxp1);
                        fxp1 = fxp_add(n1, n3);
                        tgt2 = get_target(ldx + ldz);
                        test_fxp("add   (n1+n3)", tgt2, fxp1);

                        fxp1 = fxp_mul(n1, n2);
                        fxp2 = fxp_mul_l(n1, n2);
                        tgt1 = get_target(ldx * ldy);
                        test_fxp("mul_l (n1*n2)", tgt1, fxp2);
                        test_fxp("mul   (n1*n2)", tgt1, fxp1);
                        fxp1 = fxp_mul(n1, n3);
                        fxp2 = fxp_mul_l(n1, n3);
                        tgt2 = get_target(ldx * ldz);
                        test_fxp("mul_l (n1*n3)", tgt2, fxp2);
                        test_fxp("mul   (n1*n3)", tgt2, fxp1);

                        fxp1 = fxp_div_l(n1, n2);
                        fxp2 = fxp_div(n1, n2);
                        tgt1 = get_div_target(n1, n2);
                        test_fxp("div_l (n1/n2)", tgt1, fxp1);
                        test_fxp("div   (n1/n2)", tgt1, fxp2);
                        fxp1 = fxp_div_l(n1, n3);
                        fxp2 = fxp_div(n1, n3);
                        tgt2 = get_div_target(n1, n3);
                        test_fxp("div_l (n1/n3)", tgt2, fxp1);
                        test_fxp("div   (n1/n3)", tgt2, fxp2);
                }
                if (TEST_SUPER_FXP_L) {
                        super_fxp_l snum = sfxp_l_from_fxp(n1);
                        int vn1 = sfxp_l_2_fxp(snum);
                        test_fxp("n1<->sfxp_l", fxp2ld(n1), vn1);
                }

                // Make arguments positive for the log tests
                n1 = n1 < 0? -n1: n1;
                n2 = n2 < 0? -n2: n2;
                n3 = n3 < 0? -n3: n3;
                printf("n1,n2,n3 made positive for logs and powers\n");
                if (TEST_LOGARITHMS) {
                        printf("Testing lg2 and pow2 implementations, ");
                        printf("frac bits: %d\n", frac_bits);

                        test_lg2("n1)", n1);
                        test_lg2("n2)", n2);
                        test_lg2("n3)", n3);
                        if ((whole_bits >= 3) && (TEST_LG2_MUL_L)) {
                                test_lg2_mul_l("n1)", n1);
                                test_lg2_mul_l("n2)", n2);
                                test_lg2_mul_l("n3)", n3);
                        }
                        test_ln("n1)", n1);
                        test_ln("n2)", n2);
                        test_ln("n3)", n3);
                        test_lg10("n1)", n1);
                        test_lg10("n2)", n2);
                        test_lg10("n3)", n3);
                }

                // Test powers with rand nums
                if (TEST_POWERS) {
                        n4 = -fxp_get_frac_part_bin(n1);
                        test_pow2("-frac(n1))", n4);
                        test_pow2("n3)", n3);
                        test_pow2("-n3)", -n3);
                        test_exp("-frac(n1))", n4);
                        test_exp("n3)", n3);
                        test_exp("-n3)", -n3);
                        test_pow10("-frac(n1))", n4);
                        test_pow10("n3)", n3);
                        test_pow10("-n3)", -n3);
                }
                if (TEST_SQRT) {
                        test_sqrt("n1)", n1);
                        test_sqrt("n2)", n2);
                        test_sqrt("n3)", n3);
                }
                if (TEST_POWXY) {
                        test_powxy("n1, n2)", n1, n2);
                        test_powxy("n1, -n3)", n1, -n3);
                        test_powxy("n2, -n3)", n2, -n3);
                }
        }
}

int main(void)
{
        printf("\n%sFXP Tester run\n%s", DASHES, DASHES);
        print_sys_info();

        if (SET_RAND_SEED) srand((unsigned int) time(0));

        // Make sure starting trascendental constants have exactly
        // same values after resetting the # frac bits to default
        int nfb = FXP_frac_bits;
        unsigned int e = fxp_get_e();
        unsigned int pi = fxp_get_pi();
        unsigned int ln2 = fxp_get_ln_2();
        unsigned int lg10_2 = fxp_get_lg10_2();

        fxp_set_frac_bits( nfb );
        if ((e != fxp_get_e()) || (pi != fxp_get_pi()) || \
            (ln2 != fxp_get_ln_2()) || (lg10_2 != fxp_get_lg10_2())) {
                // These differences might happen if rounding for least-significant
                // bit of the constants is done differently between the startup vs.
                // what fxp_set_frac_bits() does when called
                printf("Constant difference detected after resetting frac bits to default:\n");
                printf("  \tstarting\tafter\n");
                printf("e:\t%x\t\t%x\n", e, fxp_get_e());
                printf("pi:\t%x\t\t%x\n", pi, fxp_get_pi());
                printf("ln2:\t%x\t\t%x\n", ln2, fxp_get_ln_2());
                printf("lg10_2:\t%x\t\t%x\n", lg10_2, fxp_get_lg10_2());
                assert( 0 );
        }

        printf("+INF_LD : %.10LE\n", FXP_PINF_LD);
        printf("+INF_D  : %.10lE\n", FXP_PINF_D);
        printf("+INF_F  : %.10E\n", FXP_PINF_F);
        printf("-INF_F  : %.10E\n", FXP_NINF_F);
        printf("-INF_D  : %.10lE\n", FXP_NINF_D);
        printf("-INF_LD : %.10LE\n", FXP_NINF_LD);
        printf("UNDEF_F : %.10E\n", FXP_UNDEF_F);
        printf("UNDEF_D : %.10lE\n", FXP_UNDEF_D);
        printf("UNDEF_LD: %.10LE\n", FXP_UNDEF_LD);
        long double x = FXP_PINF_LD;

        // Any of these tests should trigger the assert
        //test_fxp("\nMust fail 01: +inf vs. -inf", FXP_PINF_LD, FXP_NEG_INF);
        //test_fxp("\nMust fail 02: -inf vs. +inf", FXP_NINF_LD, FXP_POS_INF);
        //test_fxp("\nMust fail 03: +inf vs. undef", FXP_PINF_LD, FXP_UNDEF);
        //test_fxp("\nMust fail 04: -inf vs. undef", FXP_NINF_LD, FXP_UNDEF);
        //test_fxp("\nMust fail 05: undef vs. +inf", FXP_UNDEF_LD, FXP_POS_INF);
        //test_fxp("\nMust fail 06: undef vs. -inf", FXP_UNDEF_LD, FXP_NEG_INF);
        //test_fxp("\nMust fail 07: +inf vs. finite far from inf threshold", \
        //            FXP_PINF_LD, d2fxp(-0.5));
        //test_fxp("\nMust fail 08: -inf vs. finite far from inf threshold", \
        //            FXP_NINF_LD, d2fxp(0.5));
        //test_fxp("\nMust fail 09: finite far from inf threshold vs. +inf", \
        //            -0.5L, FXP_POS_INF);
        //test_fxp("\nMust fail 10: finite far from inf threshold vs. -inf", \
        //            0.5L, FXP_NEG_INF);

        // Testing match between tester's long double vs. FXP's infinities
        test_fxp("Testing +INF from Tester vs. FXP", FXP_PINF_LD, FXP_POS_INF);
        test_fxp("Testing -INF from Tester vs. FXP", FXP_NINF_LD, FXP_NEG_INF);
        test_fxp("Testing UNDEF from Tester vs. FXP", FXP_UNDEF_LD, FXP_UNDEF);

        for (int i = 0; i < FXP_INT_BITS; i++) {
                fracbit_nwarnings[i] = 0;
                fracbit_threshold_cases[i] = 0;
                fracbit_extr_powxy_cases[i] = 0;
                fracbit_maxdelta_allowed[i] = 0.0;
                fracbit_maxdelta_observed[i] = 0.0;
        }
        int nconfigs = sizeof(fracbit_configs) / sizeof(fracbit_configs[0]);
        printf("\nRunning tests for frac-bit sizes: ");
        for (nfb = 0; nfb < nconfigs; nfb++) {
                printf("%d", fracbit_configs[nfb]);
                if ((nfb + 1) < nconfigs) printf(", ");
        }
        printf("\n");

        for (int nfb = 0; nfb < nconfigs; nfb++) {
                fxp_set_frac_bits(fracbit_configs[nfb]);
                fxp_set_auto_frac_max_dec();

                frac_bits = FXP_frac_bits;
                whole_bits = FXP_whole_bits;
                frac_mask = FXP_frac_mask;
                frac_max = FXP_frac_max;
                frac_max_dec = fxp_frac_max_dec;
                whole_max = FXP_whole_max;
                whole_min = FXP_whole_min;
                fxp_largest = FXP_MAX;

                zero = fxp(0);
                fxp_ten = fxp_bin(10, 0);
                fxp_two = fxp_bin(2, 0);
                fxp_one = fxp_bin(1, 0);
                fxp_point5 = fxp_one >> 1;
                fxp_halfmax = FXP_MAX >> 1;
                fxp_halfp2 = fxp_add(fxp_halfmax, fxp_two);

                max_warn_delta = ((long double) WDELTA_MAX) / frac_max;
                min_warn_delta = max_warn_delta / WDELTA_DIV;
                fracbit_maxdelta_allowed[frac_bits] = max_warn_delta;

                nwarnings = 0;
                larger_delta = 0;
                printf("\n%sFXP configuration parameters:\n", DASHES);
                printf("\tfrac bits : %d (requested was %d)\n",
                                frac_bits, fracbit_configs[nfb]);
                printf("\twhole bits: %d\n", whole_bits);
                printf("\t+INF      : x%X (wh,fr: %d, %d,  Ld: %.10LE)\n", \
                                FXP_POS_INF, \
                                fxp_get_whole_part(FXP_POS_INF), \
                                fxp_get_frac_part_bin(FXP_POS_INF), \
                                FXP_PINF_LD);
                printf("\tlargest   : x%X (Ld: %.10LE)\n", \
                                (FXP_MAX),
                                FXP_max_ld);
                //printf("largest <<: %ld (x%lX, largest l-shifted and rounded, for mul_l)\n", \
                //                FXP_max_lshifted, FXP_max_lshifted);
                printf("\twhole max : x%X\n", whole_max);
                printf("\tfrac mask : x%X\n", frac_mask);
                printf("\tfrac max  : x%X (->decimals: .%d)\n",
                                frac_max,
                                frac_max_dec);
                printf("\ttiniest   : x%X (wh,fr: %d, %d,  Ld: %.10LE)\n", \
                                1, 0, 1, fxp2ld(1));

                printf("\twhole min : x%X\n", whole_min);
                printf("\t-largest  : x%X (Ld: %1.3LE)\n", \
                                FXP_MIN,
                                FXP_min_ld);
                printf("\t-INF      : x%X (wh,fr: %d, %d,  Ld: %1.10LE)\n", \
                                FXP_NEG_INF, \
                                fxp_get_whole_part(FXP_NEG_INF), \
                                fxp_get_frac_part_bin(FXP_NEG_INF), \
                                FXP_NINF_LD);
                printf("\tUNDEF     : x%X (wh,fr: %d, %d,  Ld: %1.10LE)\n", \
                                FXP_UNDEF, \
                                fxp_get_whole_part(FXP_UNDEF), \
                                fxp_get_frac_part_bin(FXP_UNDEF), \
                                FXP_UNDEF_LD);

                printf("\nFXP Tester-related floating point values:\n");
                dfxp_tiniest = fxp2ld(FXP_TINIEST);
                printf("\tmin_warn_delta: %.10LE (== %1.1Lfx tiniest)\n", \
                                min_warn_delta, min_warn_delta/dfxp_tiniest);
                printf("\tmax_warn_delta: %.10LE (== %1.1Lfx tiniest)\n", \
                                max_warn_delta, max_warn_delta/dfxp_tiniest);

                printf("\tLower & upper bound values per type:\n");
                printf("\tlong double:\n");
                printf("\t\tFXP_max_ldx: %.10LE\n", FXP_max_ldx);
                printf("\t\tFXP_max_ld : %.10LE\n", FXP_max_ld);
                printf("\t\tFXP_min_ld : %.10LE\n", FXP_min_ld);
                printf("\t\tFXP_min_ldx: %.10LE\n", FXP_min_ldx);
                printf("\tdouble:\n");
                printf("\t\tFXP_max_dx : %.10E\n", FXP_max_dx);
                printf("\t\tFXP_min_dx : %.10E\n", FXP_min_dx);
                printf("\tfloat:\n");
                printf("\t\tFXP_max_fx : %.10E\n", FXP_max_fx);
                printf("\t\tFXP_min_fx : %.10E\n", FXP_min_fx);
                printf("%s\n", DASHES);

                // Tests to run =============================
                if (TEST_BASICS) {
                        tests_01();
                        tests_02();
                        tests_03();
                        test_decbin_mappings();
                        test_fracs();
                        test_ops_with_whole_bits();
                        test_ops_with_values_of_interest();
                }
                if (TEST_SUPER_FXP_L) {
                        test_super_fxp_l();
                }
                if (TEST_LOGARITHMS) {
                        test_logarithms();
                }
                if (TEST_POWERS) {
                        test_powers();
                }
                if (TEST_SQRT) {
                        test_sqroots();
                }
                if (TEST_POWXY) {
                        test_powersxy();
                }
                if (TEST_TRIGONOM) {
                        test_trigonometrics();
                }
                if (TEST_WITH_RANDS) {
                        test_ops_with_rand_nums();
                }
                // ==========================================

                printf("\n%d Warnings for %d frac bits.\n", \
                            nwarnings, frac_bits);
                if (nwarnings > 0) {
                        printf("Largest delta was: %1.2LE (max allowed: %1.2LE)\n", \
                            larger_delta, max_warn_delta);
                }
                printf("All tests passed using %d-bit fracs, ", frac_bits);
                printf("and '%d' as max decimal frac.\n\n", frac_max_dec);
                twarnings += nwarnings;
                if (larger_delta > largest_delta) {
                        largest_delta = larger_delta;
                        largest_madelta = max_warn_delta;
                        largest_delta_fbits = frac_bits;
                }
        }

        printf("\nSummary of FXP Tests:\n");
        printf("Grand total of %d warnings checking %d configurations.\n", \
                    twarnings, nconfigs);
        printf("F.bits  #Warnings  #Thr.cases  Max Delta   Max Delta Allowed   Extr.Powxy Cases\n");
        for (int nfb = 0; nfb < nconfigs; nfb++) {
                int fb = fracbit_configs[nfb];
                printf("%d\t %5d\t    %5d      %1.5Le   %1.5Le      %5d\n", \
                        fb,
                        fracbit_nwarnings[fb],
                        fracbit_threshold_cases[fb],
                        fracbit_maxdelta_observed[fb],
                        fracbit_maxdelta_allowed[fb],
                        fracbit_extr_powxy_cases[fb]);
        }

        return 0;
}
