/* SPDX-License-Identifier: MIT */
/*
 * fxp_tester.c
 * Tests the implementation of binary fixed point numbers
 * (fxp.c)
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

#define DASHES "========================\n"

// Set to 0 in order to be able to replicate runs
#define SET_RAND_SEED 1

#define TEST_WITH_RANDS 1
#define MAX_RAND_NUMS 5
//#define MAX_RAND_NUMS 5000
#define TEST_LG2_MUL_L 0

#define WDELTA 2.0
#define WDELTA_MAX 2.0

static int fracbit_configs[] = {8, 11, 13, 16, 24, 31};
//static int fracbit_configs[] = {16};
/*
static int fracbit_configs[] = {
4, 8, 9, 10, 11, 12,
13, 14, 15, 16, 17, 18, 19, 20,
21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31};
*/

static long double max_warn_delta = 0.0;
static long double larger_delta = 0.0;
static long double largest_delta = 0.0;
static long double largest_madelta = 0.0;
static int largest_delta_fbits = 0;
static long double warn_delta = 0.0;
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

void test_fxp(char *s, long double d_assert_val, int fxp1)
{
        printf("%s\n", s);
        printf("    exp: ");

        long double avplim = lim_frac(d_assert_val, frac_bits);
        if (d_assert_val == FXP_UNDEF_LD)
            printf("UNDEF");
        else if (d_assert_val == FXP_NINF_LD)
            printf("-INF");
        else if (d_assert_val == FXP_PINF_LD)
            printf("+INF");
        else {
            printf("%.10LE", avplim);
            //printf("\n\t     %.20LE", d_assert_val);
            //printf("\n\t     %.20LE", FXP_PINF_LD);
            //printf("\n\tdiff 1: %LE", (FXP_PINF_LD - d_assert_val));
            //printf("\n\tdiff 2: %LE\n", (d_assert_val - FXP_PINF_LD));
        }

        printf("\n    act: "); print_fxp(fxp1);
        int b1 = ((fxp1 == FXP_UNDEF) && (d_assert_val == FXP_UNDEF_LD));
        int b2 = ((fxp1 == FXP_NEG_INF) && (d_assert_val == FXP_NINF_LD));
        int b3 = ((fxp1 == FXP_POS_INF) && (d_assert_val == FXP_PINF_LD));
        if (b1 || b2 || b3) {
            printf(" (~same)\n"); //b1:%d, b2:%d, b3:%d\n", b1, b2, b3);
            return;
        }

        long avlimw = (long) avplim;
        long double df = fxp2ld(fxp1);
        long double delta = (df >= avplim)?
                        df - avplim:
                        avplim - df;

        if (delta <= warn_delta) {
            // No warning up to the warn_delta value
            printf(" (~same)\n"); //delta:%.120LE\n", delta);
        } else {
            nwarnings++;
            printf("\n***** Warning %d: d=%1.2LE for %1.2LE (from %1.2LE to MAX %1.2LE allowed for %d f.bits)\n",
                        nwarnings, delta, df,
                        warn_delta,
                        max_warn_delta,
                        frac_bits);
            if (delta > larger_delta) {
                larger_delta = delta;
                printf("***** It's the largest delta so far: %LE\n", \
                        larger_delta);
            }
            //assert(0);
            printf("\n");
        }
        fflush( stdout );

        // Here assert that we never exceed the max_warn_delta
        assert( (delta <= max_warn_delta) &&
            ((df >= 0 && avplim >= 0) ||
                (df <= 0 && avplim <= 0)));
}

static long double get_target(long double x)
{
        long double y = lim_frac(x, frac_bits);
        if (y <= FXP_UNDEF_LD) return FXP_UNDEF_LD;
        if (y <= FXP_min_ld) return FXP_NINF_LD;
        if (y >= FXP_max_ld) return FXP_PINF_LD;
        return y;
}

static long double get_lg2_target(int x)
{
        if (x < 0) return FXP_UNDEF_LD;
        if (x == 0) return FXP_NINF_LD;
        if (x == FXP_POS_INF) return FXP_PINF_LD;
        return get_target(log2l(fxp2ld(x)));
}

static long double get_lg10_target(int x)
{
        long double l2 = get_lg2_target(x);
        if ((l2 == FXP_UNDEF_LD) \
                || (l2 == FXP_NINF_LD) \
                || (l2 == FXP_PINF_LD))
                return l2;
        return get_target(log10l(fxp2ld(x)));
}

static long double get_ln_target(int x)
{
        long double l2 = get_lg2_target(x);
        if ((l2 == FXP_UNDEF_LD) \
                || (l2 == FXP_NINF_LD) \
                || (l2 == FXP_PINF_LD))
                return l2;
        return get_target(logl(fxp2ld(x)));
}

static long double get_pow2_target(int x)
{
        if (x == FXP_UNDEF) return FXP_UNDEF_LD;
        if (x == FXP_NEG_INF) return 0.0L;
        if (x == FXP_POS_INF) return FXP_PINF_LD;
        return get_target(powl(2.0L, fxp2ld(x)));
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
                        target = lim_frac(ldx, frac_bits) \
                                    / lim_frac(ldy, frac_bits);
                        target = lim_frac(target, frac_bits);
                        if (target > FXP_max_ld) {
                            target = FXP_PINF_LD;
                        } else {
                            if (target < FXP_min_ld) {
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
                    if (target < FXP_min_ld)
                        target = FXP_NINF_LD;
                    else if (target > FXP_max_ld)
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
        if (f < FXP_min_f) return FXP_NINF_LD;
        if (f > FXP_max_f) return FXP_PINF_LD;
        return (long double) f;
}

static int FXP_TINIEST = 1;
static long double dfxp_tiniest;

static int zero;
static int fxp_ten;
static int fxp_two;
static int fxp_one;
static int fxp_p5;
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
        test_fxp("Largest frac",
                    fxp2ld(fxp_bin(0, frac_max)),
                    frac_max);
        test_fxp("tiniest",
                    dfxp_tiniest,
                    FXP_TINIEST);
        test_fxp("0.5",
                    ((long double) 0.5),
                    fxp_p5);
        test_fxp("zero",
                    0.0,
                    fxp_bin(0, 0));
        test_fxp("-tiniest",
                    fxp2ld(-FXP_TINIEST),
                    -FXP_TINIEST);
        test_fxp("-Largest frac",
                    fxp2ld(-frac_max),
                    -frac_max);
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
                    fxp_sub(FXP_POS_INF, fxp_p5));
        test_fxp("-inf + 0.5",
                    FXP_NINF_LD,
                    fxp_add(FXP_NEG_INF, fxp_p5));
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
                    fxp_sub( -fxp_largest, fxp_p5));
        test_fxp("Unsafe Too neg substraction",
                    fxp2ld(-fxp_largest - fxp_p5),
                    fxp_unsafe_sub( -fxp_largest, fxp_p5));
        test_fxp("Largest + 0.5",
                    FXP_PINF_LD,
                    fxp_add(fxp_largest, fxp_p5));
        test_fxp("-Largest - 0.5",
                    FXP_NINF_LD,
                    fxp_add(-fxp_largest, -fxp_p5));
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
                    fxp_mul(fxp_add(fxp_halfmax, fxp_p5), fxp(2)));
        test_fxp("(HalfMax+0.5)*2 (long)",
                    FXP_PINF_LD,
                    fxp_mul_l(fxp_add(fxp_halfmax, fxp_p5), fxp(2)));
}

void test_decbin_mappings()
{
        printf("\nChecking decimal <=> bin mappings of frac ranges, ");
        printf("frac bits: %d\n", frac_bits);
        int fmbin = frac_max;
        int fmdec = frac_max_dec;
        printf("Max frac dec: %d (bin %d)", fmdec, fmbin);
        for (int i = 0; i <= 5; i++) {
                printf("\nShowing fxp for 0.%3d: ", i);
                int vf = fxp_dec(0, i);
                print_fxp(vf);
        }
        printf("\n:");
        int m = (fmdec + 1) / 2;
        for (int i = (m >= 2 ? m - 2: 0); i <= m + 2; i++) {
                printf("\nShowing fxp for 0.%3d: ", i);
                int vf = fxp_dec(0, i);
                print_fxp(vf);
        }
        printf("\n:");
        for (int i = (fmdec >= 5? fmdec - 5: 0);
                    i <= fmdec; i++) {
                printf("\nShowing fxp for 0.%3d: ", i);
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
        long double bkp_wd = warn_delta;
        long double bkp_ld = larger_delta;
        long double bkp_mwd = max_warn_delta;
        int fmdec = fxp_get_frac_max_dec();

        fxp_set_frac_max_dec(999);
        warn_delta = lim_frac(fxp2ld(fxp_dec(0,100)), frac_bits);
        max_warn_delta = lim_frac(warn_delta * 2, frac_bits);

        test_fxp("-0.(+)500",  0.5, fxp_dec(-0, 500));
        test_fxp("-0.(-)500", -0.5, fxp_dec(-0, -500));

        printf("\nChecking truncation of longer frac decimal arguments, ");
        printf("frac bits: %d\n", frac_bits);
        test_fxp("0.22222",     0.222, fxp_dec(0, 22222));
        test_fxp("0.4444444",   0.444, fxp_dec(0, 4444444));
        test_fxp("0.771999",    0.771, fxp_dec(0, 771999));
        test_fxp("0.9999999",   0.999, fxp_dec(0, 9999999));
        // Restoring original frac_max_dec and delta vars
        fxp_set_frac_max_dec(fmdec);
        warn_delta = bkp_wd;
        larger_delta = bkp_ld;
        max_warn_delta = bkp_mwd;

        printf("\nChecking extreme frac values, ");
        printf("frac bits: %d\n", frac_bits);
        test_fxp("fracmax",
                    fxp2ld(frac_max),
                    frac_max);
        test_fxp("0.5 + fracmax",
                    ((whole_bits > 1)?
                        0.5 + fxp2ld(frac_max):
                        fxp2ld(FXP_POS_INF)),
                    fxp_add(fxp_p5, frac_max));
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
        int bin_frac = frac_max / 2;  // == 0.500
        int num = fxp_bin(whole, bin_frac);
        int dec_frac = fxp_get_dec_frac(num);
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
                    FXP_UNDEF, FXP_NEG_INF,
                    FXP_POS_INF,
                    FXP_MAX,
                    fxp_dec(33333, 333), fxp_dec(33, 33),
                    fxp(2),
                    fxp(1),
                    fxp_bin(0, FXP_frac_max),
                    1, 0, fxp_bin(0, -1),
                    fxp_bin(0, -FXP_frac_max),
                    fxp(-1),
                    fxp(-2),
                    fxp_dec(-33, 33), fxp_dec(-33333, 333),
                    FXP_MIN
                    };
        int ay[] = {
                    FXP_UNDEF, FXP_NEG_INF, FXP_POS_INF, FXP_MAX,
                    fxp(2), fxp(1),
                    fxp_bin(0, FXP_frac_max),
                    2, 1, 0, -1, -2,
                    fxp(-1), fxp(-2), FXP_MIN
                    };
        int x, y, posx, posy;
        int ndd = (int) (sizeof(ax) / sizeof(int));
        int ndr = (int) (sizeof(ay) / sizeof(int));
        for (int i = 0; i < ndd; i++) {
                x = ax[i];
                fflush( stdout );
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
                printf("x    : "); print_fxp(x); printf("\n");
                printf("float: %fE\n", fx);
                printf("x2   : "); print_fxp(x2); printf("\n");
                test_fxp("\nx: ", get_f_target(x), x2);

                for (int j=0; j<ndr; j++) {
                        y = ay[j];
                        if ((j > 2) && ((y == FXP_POS_INF) \
                                || (y == FXP_NEG_INF))) {
                                continue;
                        }
                        ldx = fxp2ld(x);
                        ldy = fxp2ld(y);
                        printf("y: "); print_fxp(y);

                        //For multiplication
                        tgt1 = get_mul_target(x, y);
                        //n1 = fxp_mul_l(x, y);
                        //test_fxp("\nmul_l  (x*y)", n1, tgt1);

                        n1 = fxp_mul(x, y);
                        test_fxp("\nmul  (x*y)", tgt1, n1);

                        //For division
                        tgt1 = get_div_target(x, y);
                        n1 = fxp_div(x, y);
                        test_fxp("div  (x/y)", tgt1, n1);
                }
        }
}

void test_lg2_mul_l(char * msg, int x)
{
        printf("lg2_mul_l(");
        test_fxp(msg, get_lg2_target(x), fxp_lg2_mul_l(x));
}

void test_lg2(char * msg, int x)
{
        long double tgt = get_lg2_target(x);
        printf("lg2_l("); test_fxp(msg, tgt, fxp_lg2_l(x));
        printf("lg2(");   test_fxp(msg, tgt, fxp_lg2(x));
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

void test_logarithms()
{
        printf("\nShowing Transcendental constants as fxp's: ");
        printf("frac bits: %d\n", frac_bits);
        printf("e      : "); print_fxp(fxp_get_e()); printf("\n");
        printf("pi     : "); print_fxp(fxp_get_pi()); printf("\n");
        printf("ln(2)  : "); print_fxp(fxp_get_ln_2()); printf("\n");
        printf("lg10(2): "); print_fxp(fxp_get_lg10_2()); printf("\n");

        if ((whole_bits >= 3) && (TEST_LG2_MUL_L)) {
                // To test logarithm calculation with fxp_log2_mul_l
                // we need an fxp configuration with at least 3 whole bits
                test_lg2_mul_l("\n+INF):",        FXP_POS_INF);
                test_lg2_mul_l("largest):",     FXP_MAX);
                test_lg2_mul_l("2.2):",         fxp_bin(2, FXP_frac_mask/5));
                test_lg2_mul_l("2):",           fxp(2));
                test_lg2_mul_l("1.99...9):",    fxp_bin(1, FXP_frac_max));
                test_lg2_mul_l("1.00..01):",    fxp_bin(1, 1));
                test_lg2_mul_l("1):",           fxp(1));
                test_lg2_mul_l("0.99...9):",    fxp(1) - 1);
                test_lg2_mul_l("0.50..01):",    FXP_half + 1);
                test_lg2_mul_l("0.5):",         FXP_half);
                test_lg2_mul_l("tiniest):",     1);
                test_lg2_mul_l("0):",           0);
                test_lg2_mul_l("-INF):",        FXP_NEG_INF);
                test_lg2_mul_l("UNDEF):",       FXP_UNDEF);
        } else {
                if (TEST_LG2_MUL_L) {
                        printf("Only %d whole bit(s); 3 or more needed for lg2_mul_l; skipping.", \
                                    whole_bits);
                }
                printf("\n");
        }

        test_lg2("+INF):",      FXP_POS_INF);
        test_lg2("largest):",   FXP_MAX);
        test_lg2("100):",       fxp(100));
        test_lg2("2.2):",       fxp_bin(2, FXP_frac_mask/5));
        test_lg2("2):",         fxp(2));
        test_lg2("1.99...9):",  fxp_bin(1, FXP_frac_max));
        test_lg2("1.00..01):",  fxp_bin(1, 1));
        test_lg2("1):",         fxp(1));
        test_lg2("0.99...9):",  fxp(1) - 1);

        int exponent = (1u << (whole_bits - 1));
        int k = d2fxp(powl(2.0, -exponent ));
        if ((k > 0) && (exponent < FXP_INT_BITS)) {
                printf("\nWith only %d whole bit(s) (%d frac bits,)",
                            whole_bits, frac_bits);
                printf(" range for whole part is +/-%d,\n", FXP_whole_max);
                printf("lg2(x) must necessarily overflow returning -INF for\n");
                printf("all positive x <= k = ");
                print_fxp(k); printf("\n");
                test_lg2("k + tiniest):", k + 1);
                test_lg2("k):", k);
                test_lg2("k - tiniest):", k - 1);
                printf("\n");
        }

        test_lg2("tiniest):",   1);
        test_lg2("0):",         0);
        test_lg2("-INF):",      FXP_NEG_INF);
        test_lg2("UNDEF):",     FXP_UNDEF);

        printf("\n");
        test_ln("+INF):",       FXP_POS_INF);
        test_ln("largest):",    FXP_MAX);
        test_ln("100):",        fxp(100));
        test_ln("e+1):",        fxp_add(fxp_get_e(), fxp(1)));
        test_ln("e):",          fxp_get_e());
        test_ln("2.2):",        d2fxp(2.2));
        test_ln("2):",          fxp(2));
        test_ln("1.99...9):",   fxp_bin(1, FXP_frac_max));
        test_ln("1.5):",        d2fxp(1.5));
        test_ln("1.00..01):",   fxp_bin(1, 1));
        test_ln("1):",          fxp(1));
        test_ln("0.99...9):",   fxp(1) - 1);
        test_ln("0.50..01):",   FXP_half + 1);
        test_ln("0.5):",        FXP_half);
        test_ln("2*tiniest):",  2);
        test_ln("tiniest):",    1);
        test_ln("0):",          0);
        test_ln("-INF):",       FXP_NEG_INF);
        test_ln("UNDEF):",      FXP_UNDEF);

        printf("\n");
        test_lg10("+INF):",     FXP_POS_INF);
        test_lg10("largest):",  FXP_MAX);
        test_lg10("100):",      fxp(100));
        test_lg10("10):",       fxp(10));
        test_lg10("2.2):",      fxp_bin(2, FXP_frac_mask/5));
        test_lg10("2):",        fxp(2));
        test_lg10("1.99...9):", fxp_bin(1, FXP_frac_max));
        test_lg10("1.00..01):", fxp_bin(1, 1));
        test_lg10("1):",        fxp(1));
        test_lg10("0.99...9):", fxp(1) - 1);
        test_lg10("0.50..01):", FXP_half + 1);
        test_lg10("0.5):",      FXP_half);
        test_lg10("2*tiniest):",2);
        test_lg10("tiniest):",  1);
        test_lg10("0):",        0);
        test_lg10("-INF):",     FXP_NEG_INF);
        test_lg10("UNDEF):",    FXP_UNDEF);
}

void test_pow2(char * msg, int x)
{
        long double tgt = get_pow2_target(x);
        printf("pow2_l("); test_fxp(msg, tgt, fxp_pow2_l(x));
        printf("pow2(");   test_fxp(msg, tgt, fxp_pow2(x));
}

void test_exp(char * msg, int x)
{
        long double tgt = get_exp_target(x);
        printf("exp_l("); test_fxp(msg, tgt, fxp_exp_l(x));
        //printf("exp(");   test_fxp(msg, tgt, fxp_exp(x));
}

void test_pow10(char * msg, int x)
{
        long double tgt = get_pow10_target(x);
        printf("pow10_l("); test_fxp(msg, tgt, fxp_pow10_l(x));
        //printf("exp(");   test_fxp(msg, tgt, fxp_exp(x));
}



void test_powers()
{
        printf("\n");
        test_pow2("largest):",      FXP_MAX);
        test_pow2("14.999...):",    fxp_bin(14, FXP_frac_max));
        test_pow2("3.5):",          d2fxp(3.5));
        test_pow2("3):",            fxp(3));
        test_pow2("2.9):",          d2fxp(2.9));
        test_pow2("2):",            fxp(2));
        test_pow2("1+tiniest):",    fxp(1) + 1);
        test_pow2("1):",            fxp(1));
        test_pow2("1-tiniest):",    fxp(1) - 1);
        test_pow2("0.5):",          d2fxp(0.5));
        test_pow2("2*tiniest):",    2);
        test_pow2("tiniest):",      1);
        test_pow2("0):",            0);
        test_pow2("-tiniest):",     -1);
        test_pow2("-2*tiniest):",   -2);
        test_pow2("-4*tiniest):",   -4);
        test_pow2("-0.5):",         d2fxp(-0.5));
        test_pow2("-1+tiniest):",   fxp(-1) + 1);
        test_pow2("-1):",           fxp(-1));
        test_pow2("-1-tiniest):",   fxp(-1) - 1);
        test_pow2("-1.5):",         d2fxp(-1.5));
        test_pow2("-2):",           fxp(-2));
        test_pow2("-8):",           fxp(-8));
        test_pow2("-16):",          fxp(-16));
        test_pow2("-30):",          fxp(-30));
        test_pow2("-31):",          fxp(-31));
        test_pow2("-32):",          fxp(-32));
        test_pow2("-42):",          fxp(-42));
        test_pow2("most neg):",     -fxp_largest);
        test_pow2("-INF):",         FXP_NEG_INF);

        printf("\n");
        test_exp("largest):",      FXP_MAX);
        test_exp("14.999...):",    fxp_bin(14, FXP_frac_max));
        test_exp("3.5):",          d2fxp(3.5));
        test_exp("2):",            fxp(2));
        test_exp("1.5):",          d2fxp(1.5));
        test_exp("1+tiniest):",    fxp(1) + 1);
        test_exp("1):",            fxp(1));
        test_exp("1-tiniest):",    fxp(1) - 1);
        test_exp("0.5):",          d2fxp(0.5));
        test_exp("2*tiniest):",    2);
        test_exp("tiniest):",      1);
        test_exp("0):",            0);
        test_exp("-tiniest):",     -1);
        test_exp("-2*tiniest):",   -2);
        test_exp("-4*tiniest):",   -4);
        test_exp("-0.5):",         d2fxp(-0.5));
        test_exp("-1+tiniest):",   fxp(-1) + 1);
        test_exp("-1):",           fxp(-1));
        test_exp("-1-tiniest):",   fxp(-1) - 1);
        test_exp("-2):",           fxp(-2));
        test_exp("-3.5):",         d2fxp(-3.5));
        test_exp("-8):",           fxp(-8));
        test_exp("-16):",          fxp(-16));
        test_exp("-30):",          fxp(-30));
        test_exp("-32):",          fxp(-32));
        test_exp("most neg):",     -fxp_largest);

        printf("\n");
        test_pow10("largest):",      FXP_MAX);
        test_pow10("14.999...):",    fxp_bin(14, FXP_frac_max));
        test_pow10("3.5):",          d2fxp(3.5));
        test_pow10("2):",            fxp(2));
        test_pow10("1.5):",          d2fxp(1.5));
        test_pow10("1+tiniest):",    fxp(1) + 1);
        test_pow10("1):",            fxp(1));
        test_pow10("1-tiniest):",    fxp(1) - 1);
        test_pow10("0.5):",          d2fxp(0.5));
        test_pow10("2*tiniest):",    2);
        test_pow10("tiniest):",      1);
        test_pow10("0):",            0);
        test_pow10("-tiniest):",     -1);
        test_pow10("-2*tiniest):",   -2);
        test_pow10("-4*tiniest):",   -4);
        test_pow10("-0.5):",         d2fxp(-0.5));
        test_pow10("-1+tiniest):",   fxp(-1) + 1);
        test_pow10("-1):",           fxp(-1));
        test_pow10("-1-tiniest):",   fxp(-1) - 1);
        test_pow10("-2):",           fxp(-2));
        test_pow10("-3.5):",         d2fxp(-3.5));
        test_pow10("-8):",           fxp(-8));
        test_pow10("-16):",          fxp(-16));
        test_pow10("-30):",          fxp(-30));
        test_pow10("-32):",          fxp(-32));
        test_pow10("most neg):",     -fxp_largest);
}

void test_ops_with_rand_nums()
{
        int sign1, sign2, sign3, n1, n2, n3, n4, fxp1, fxp2;
        int vn1, vn2, vn3;
        long double ldx, ldy, ldz, tgt1, tgt2, tgt3;

        for (int i=0; i < MAX_RAND_NUMS; i++) {
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

                printf("\nRandom Test #%d, ", i);
                printf("frac bits: %d\n", frac_bits);

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
                test_fxp("n1: ", get_f_target(n1), f2fxp(fxp2f(n1)));
                test_fxp("n2: ", get_f_target(n2), f2fxp(fxp2f(n2)));
                test_fxp("n3: ", get_f_target(n3), f2fxp(fxp2f(n3)));

                fxp1 = fxp_add(n1, n2);
                tgt1 = get_target(ldx + ldy);
                test_fxp("add   (n1+n2)", tgt1, fxp1);
                fxp1 = fxp_add(n1, n3);
                tgt2 = get_target(ldx + ldz);
                test_fxp("add   (n1+n3)", tgt2, fxp1);

                fxp1 = fxp_mul(n1, n2);
                fxp2 = fxp_mul_l(n1, n2);
                tgt1 = get_target(ldx * ldy);
                test_fxp("mul   (n1*n2)", tgt1, fxp1);
                test_fxp("mul_l (n1*n2)", tgt1, fxp2);
                fxp1 = fxp_mul(n1, n3);
                fxp2 = fxp_mul_l(n1, n3);
                tgt2 = get_target(ldx * ldz);
                test_fxp("mul   (n1*n3)", tgt2, fxp1);
                test_fxp("mul_l (n1*n3)", tgt2, fxp2);

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

                // Make arguments positive for the log tests
                n1 = n1 < 0? -n1: n1;
                n2 = n2 < 0? -n2: n2;
                n3 = n3 < 0? -n3: n3;
                printf("Testing lg2 and pow2 implementations, ");
                printf("frac bits: %d", frac_bits);
                printf("\nn1 = "); print_fxp(n1);
                printf("\nn2 = "); print_fxp(n2);
                printf("\nn3 = "); print_fxp(n3); printf("\n");

                tgt1 = get_lg2_target(n1);
                tgt2 = get_lg2_target(n2);
                tgt3 = get_lg2_target(n3);
                if ((whole_bits >= 3) && (TEST_LG2_MUL_L)) {
                        test_lg2_mul_l("n1)", n1);
                        test_lg2_mul_l("n1)", n2);
                        test_lg2_mul_l("n1)", n3);
                }
                test_lg2("n1)", n1);
                test_lg2("n2)", n2);
                test_lg2("n3)", n3);

                test_ln("n1)", n1);
                test_ln("n2)", n2);
                test_ln("n3)", n3);

                test_lg10("n1)", n1);
                test_lg10("n2)", n2);
                test_lg10("n3)", n3);

                // Test powers with rand nums
                printf("\n");
                n4 = -fxp_get_bin_frac(n1);
                tgt1 = get_pow2_target(n4);
                tgt2 = get_pow2_target(n3);
                tgt3 = get_pow2_target(-n3);
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
}

int main(void)
{
        printf("\n%sFXP Tester run\n%s", DASHES, DASHES);
        print_sys_info();

        if (SET_RAND_SEED) srand((unsigned int) time(0));

        // Make sure starting trascendental constants still have exactly
        // the same values after resetting the # of frac bits to default
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

        printf("+INF_LD : %LE\n", FXP_PINF_LD);
        printf("+INF_D  : %lE\n", FXP_PINF_D);
        printf("+INF_F  : %E\n", FXP_PINF_F);
        printf("-INF_F  : %E\n", FXP_NINF_F);
        printf("-INF_D  : %lE\n", FXP_NINF_D);
        printf("-INF_LD : %LE\n", FXP_NINF_LD);
        printf("UNDEF_F : %E\n", FXP_UNDEF_F);
        printf("UNDEF_D : %lE\n", FXP_UNDEF_D);
        printf("UNDEF_LD: %LE\n", FXP_UNDEF_LD);

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
                fxp_p5 = fxp_bin(0, frac_max/2 + 1);
                fxp_halfmax = FXP_MAX / 2;
                fxp_halfp2 = fxp_add(fxp_halfmax, fxp_two);

                //warn_delta = ((long double) MIN(frac_bits/2, WDELTA)) / frac_max;
                warn_delta = ((long double) frac_bits) / 2.0;
                if (WDELTA < warn_delta) warn_delta = WDELTA;
                warn_delta /= ((long double) frac_max);

                max_warn_delta = lim_frac(warn_delta * WDELTA_MAX,
                                            frac_bits);
                warn_delta = lim_frac(warn_delta, frac_bits);

                nwarnings = 0;
                larger_delta = 0;
                printf("\n%sFXP configuration parameters:\n", DASHES);
                printf("frac bits : %d (requested was %d)\n",
                                frac_bits, fracbit_configs[nfb]);
                printf("whole bits: %d\n", whole_bits);
                printf("+INF      : %d (x%x,  wh,fr values: %d, %d,  Lf ~ %1.3LE)\n", \
                                FXP_POS_INF, FXP_POS_INF, \
                                fxp_get_whole_part(FXP_POS_INF), \
                                fxp_get_bin_frac(FXP_POS_INF), \
                                FXP_PINF_LD);
                printf("largest   : %d (fxp lf ~ %1.3le)\n", \
                                (FXP_MAX),
                                FXP_max_d);
                printf("largest <<: %ld (x%lX, largest l-shifted and rounded, for mul_l)\n", \
                                FXP_max_lshifted, FXP_max_lshifted);
                printf("whole max : %d\n", whole_max);
                printf("frac mask : %d\n", frac_mask);
                printf("frac max  : %d (->decimals: .%d)\n",
                                frac_max,
                                frac_max_dec);
                printf("whole min : %d\n", whole_min);
                printf("smallest  : %d (fxp lf ~ %1.3le)\n", \
                                (FXP_MIN),
                                FXP_min_d);
                printf("-INF      : %d (x%x,  wh,fr values: %d, %d,  Lf ~ %1.3LE)\n", \
                                FXP_NEG_INF, FXP_NEG_INF, \
                                fxp_get_whole_part(FXP_NEG_INF), \
                                fxp_get_bin_frac(FXP_NEG_INF), \
                                FXP_NINF_LD);
                printf("UNDEF     : %d (x%x,  wh,fr values: %d, %d,  Lf ~ %1.3LE)\n", \
                                FXP_UNDEF, FXP_UNDEF, \
                                fxp_get_whole_part(FXP_UNDEF), \
                                fxp_get_bin_frac(FXP_UNDEF), \
                                FXP_UNDEF_LD);
                dfxp_tiniest = fxp2ld(FXP_TINIEST);
                printf("tiniest   : "); print_fxp(FXP_TINIEST); printf("\n");
                printf("warn_delta: %LE (== %1.1Lfx tiniest)\n", warn_delta, warn_delta/dfxp_tiniest);
                printf("max_warn_d: %LE (== %1.1Lfx tiniest)\n", max_warn_delta, max_warn_delta/dfxp_tiniest);
                printf("FXP_max_f : %.12E\n", FXP_max_f);
                printf("FXP_min_f : %.12E\n", FXP_min_f);
                printf("FXP_max_d : %.12lE\n", FXP_max_d);
                printf("FXP_min_d : %.12lE\n", FXP_min_d);
                printf("FXP_max_ld: %.12LE\n", FXP_max_ld);
                printf("FXP_min_ld: %.12LE\n", FXP_min_ld);

                // Tests to run =============================
                tests_01();
                tests_02();
                tests_03();
                test_decbin_mappings();
                test_fracs();
                test_ops_with_whole_bits();
                test_ops_with_values_of_interest();
                test_logarithms();
                test_powers();
                if (TEST_WITH_RANDS)
                        test_ops_with_rand_nums();
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
        printf("Grand total of %d warnings checking %d configurations.\n", \
                    twarnings, nconfigs);
        if (twarnings > 0) {
            printf("Largest warning delta: %1.4LE (max allowed:%1.4LE, using %d frac bits)\n", \
                        largest_delta,
                        largest_madelta,
                        largest_delta_fbits);
        }
        return 0;
}
