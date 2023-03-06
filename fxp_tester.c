/* SPDX-License-Identifier: MIT */
/*
 * fxp_tester.c
 * Tests the implementation of binary fixed point numbers
 * and arithmetic operations (fxp.c)
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
#define SET_RAND_SEED 0

#define TEST_WITH_RANDS 1
#define MAX_RAND_NUMS 5
//#define MAX_RAND_NUMS 5000
#define TEST_LG2_MUL_L 0

#define WDELTA 2.0

// Testing the logs with bases e and 10 needs a more
// forgiving WDELTA_MAX factor, e.g. 6, otherwise the
// errors (getting those logs from the calculated
// lg2) can trigger the assert.
// When only log2 is being tested, this factor
// can be as low as 2 yet no assert gets triggered
#define WDELTA_MAX 6.0

static int fracbit_configs[] = {8, 11, 13, 16, 24, 30, 31};
//static int fracbit_configs[] = {13};
/*
static int fracbit_configs[] = {4, 8, 9, 10,   \
            11, 12, 13, 14, 15, 16, 17, 18, 19, 20,     \
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

void test_fxp(char *s, int fxp1, long double d_assert_val)
{
        printf("%s\n", s);
        printf("\tgot: "); print_fxp(fxp1);
        printf("\n\texp. ");

        long double avplim = lim_frac(d_assert_val, frac_bits);
        if (d_assert_val == FXP_UNDEF_LD)
            printf("UNDEF");
        else if (d_assert_val == FXP_NINF_LD)
            printf("-INF");
        else if (d_assert_val == FXP_PINF_LD)
            printf("+INF");
        else {
            printf("%.20LE", avplim);
            //printf("\n\t     %.20LE", d_assert_val);
            //printf("\n\t     %.20LE", FXP_PINF_LD);
            //printf("\n\tdiff 1: %LE", (FXP_PINF_LD - d_assert_val));
            //printf("\n\tdiff 2: %LE\n", (d_assert_val - FXP_PINF_LD));
        }

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
            printf("\n***** Warning %d: d=%1.2LE for %1.2LE (%1.2LE -- MAX %1.2LE allowed for %d f.bits)\n",
                        nwarnings, delta, df,
                        warn_delta,
                        max_warn_delta,
                        frac_bits);
            if (delta > larger_delta) {
                larger_delta = delta;
                printf("***** It's the largest delta so far: %LE\n", \
                        larger_delta);
            }
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
    if (y < FXP_min_ld) return FXP_NINF_LD;
    if (y > FXP_max_ld) return FXP_PINF_LD;
    return y;
}

static long double get_lg_target(int x, char base)
{
    if (x < 0) return FXP_UNDEF_LD;
    if (x == 0) return FXP_NINF_LD;
    if (x == FXP_POS_INF) return FXP_PINF_LD;
    long double l2 = get_target(log2l(fxp2ld(x)));
    if ((base == '2') || (l2 == FXP_UNDEF_LD) \
            || (l2 == FXP_NINF_LD) \
            || (l2 == FXP_PINF_LD))
            return l2;
    if (base == 'a')   // logs base 10
            return get_target(log10l(fxp2ld(x)));
    // if not base 2 or 10, then natural logs
    return get_target(logl(fxp2ld(x)));
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
                    target = lim_frac(ldx, frac_bits) / lim_frac(ldy, frac_bits);
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
        test_fxp("+Inf", FXP_POS_INF, FXP_PINF_LD);
        test_fxp("Largest", fxp_largest, fxp2ld(FXP_MAX));
        test_fxp("HalfMax", fxp_halfmax, fxp2ld(FXP_MAX)/2);
        test_fxp("Largest frac", frac_max, fxp2ld(fxp_bin(0, frac_max)));
        test_fxp("tiniest", FXP_TINIEST, dfxp_tiniest);
        test_fxp("0.5", fxp_p5, ((long double) 0.5));
        test_fxp("zero", fxp_bin(0, 0), 0.0);
        test_fxp("-tiniest", -FXP_TINIEST, fxp2ld(-FXP_TINIEST));
        test_fxp("-Largest frac", -frac_max, fxp2ld(-frac_max));
        test_fxp("Most negative", -fxp_largest, fxp2ld(FXP_MIN));
        test_fxp("-Inf", FXP_NEG_INF, FXP_NINF_LD);
        test_fxp("Undefined", FXP_UNDEF, FXP_UNDEF_LD);
}

void tests_02()
{
        printf("\nChecking extreme int values, part II, ");
        printf("frac bits: %d\n", frac_bits);

        test_fxp("Almost most negative",
                        fxp_add(-fxp_largest, FXP_TINIEST),
                        fxp2ld(-fxp_largest) + dfxp_tiniest);
        test_fxp(" Largest -Largest",
                        fxp_add(fxp_largest, -fxp_largest),
                        0.0);
        test_fxp("-Largest +Largest",
                        fxp_add(-fxp_largest, fxp_largest),
                        0.0);
        test_fxp("Largest + 0",
                        fxp_add(fxp_largest, zero),
                        fxp2ld(fxp_largest));
        test_fxp("-Largest - 0",
                        fxp_add(-fxp_largest, -zero),
                        fxp2ld(-fxp_largest));
        test_fxp("Largest - tiniest",
                        fxp_sub(fxp_largest, FXP_TINIEST),
                        fxp2ld(fxp_largest) - dfxp_tiniest);
        test_fxp("Largest + tiniest safe",
                        fxp_add(fxp_largest, FXP_TINIEST),
                        FXP_PINF_LD);
        test_fxp("Largest + tiniest unsafe",
                        fxp_unsafe_add(fxp_largest, FXP_TINIEST),
                        FXP_PINF_LD);
        test_fxp("-(+inf)",
                        -FXP_POS_INF,
                        FXP_NINF_LD);
        test_fxp("-(-inf)",
                        -FXP_NEG_INF,
                        FXP_PINF_LD);
        test_fxp("+inf + +inf",
                        fxp_add(FXP_POS_INF, FXP_POS_INF),
                        FXP_PINF_LD);
        test_fxp("-inf - +inf",
                        fxp_sub(FXP_NEG_INF, FXP_POS_INF),
                        FXP_NINF_LD);
        test_fxp("+inf + -inf",
                        fxp_add(FXP_POS_INF, FXP_NEG_INF),
                        FXP_UNDEF_LD);
        test_fxp("-inf + -inf",
                        fxp_add(FXP_NEG_INF, FXP_NEG_INF),
                        FXP_NINF_LD);
        test_fxp("-inf - -inf",
                        fxp_sub(FXP_NEG_INF, FXP_NEG_INF),
                        FXP_UNDEF_LD);
        test_fxp("+inf * -inf",
                        fxp_mul(FXP_POS_INF, FXP_NEG_INF),
                        FXP_NINF_LD);
        test_fxp("+inf - 0.5",
                        fxp_sub(FXP_POS_INF, fxp_p5),
                        FXP_PINF_LD);
        test_fxp("-inf + 0.5",
                        fxp_add(FXP_NEG_INF, fxp_p5),
                        FXP_NINF_LD);
        test_fxp("+num / zero",
                        fxp_div(fxp_largest, zero),
                        FXP_PINF_LD);
        test_fxp("zero / zero",
                        fxp_div(zero, zero),
                        FXP_UNDEF_LD);
        test_fxp("zero * zero",
                        fxp_mul(zero, zero),
                        0.0);
        test_fxp("zero + zero",
                        fxp_add(zero, zero),
                        0.0);
        test_fxp("zero - zero",
                        fxp_sub(zero, zero),
                        0.0);
        test_fxp("zero - undef",
                        fxp_sub(zero, FXP_UNDEF),
                        FXP_UNDEF_LD);
        test_fxp("-num / zero",
                        fxp_div(-fxp_largest, zero),
                        FXP_NINF_LD);
        test_fxp("zero * +inf",
                        fxp_mul(zero, FXP_POS_INF),
                        FXP_UNDEF_LD);
        test_fxp("zero * -inf",
                        fxp_mul(zero, FXP_NEG_INF),
                        FXP_UNDEF_LD);
        test_fxp("zero * undef",
                        fxp_mul(zero, FXP_UNDEF),
                        FXP_UNDEF_LD);
        test_fxp("-inf * undef",
                        fxp_mul(FXP_NEG_INF, FXP_UNDEF),
                        FXP_UNDEF_LD);
        test_fxp("+inf * undef",
                        fxp_mul(FXP_POS_INF, FXP_UNDEF),
                        FXP_UNDEF_LD);
        test_fxp("undef * undef",
                        fxp_mul(FXP_UNDEF, FXP_UNDEF),
                        FXP_UNDEF_LD);
        test_fxp("tiniest * inf",
                        fxp_mul(FXP_TINIEST, FXP_POS_INF),
                        FXP_PINF_LD);
}

void tests_03()
{
        printf("\nChecking extreme int values, part III, ");
        printf("frac bits: %d\n", frac_bits);
        test_fxp("Way Too Large whole part!",
                        fxp_add(fxp(whole_max), fxp(5)),
                        FXP_PINF_LD);
        test_fxp("Largest * 1",
                        fxp_mul(fxp_largest,  fxp_one),
                        fxp2ld(fxp_mul(fxp_largest, fxp(1))));
        test_fxp("Largest * -1",
                        fxp_mul(fxp_largest,  -fxp_one),
                        fxp2ld(fxp_mul(fxp_largest, fxp(-1))));
        test_fxp("Largest + two safe",
                        fxp_add(fxp_largest, fxp_two),
                        FXP_PINF_LD);
        test_fxp("Largest + two unsafe",
                        fxp_unsafe_add(fxp_largest, fxp_two),
                        fxp2ld(fxp_largest + fxp_two));
        test_fxp("Safe Too neg substraction",
                        fxp_sub( -fxp_largest, fxp_p5),
                        FXP_NINF_LD);
        test_fxp("Unsafe Too neg substraction",
                        fxp_unsafe_sub( -fxp_largest, fxp_p5),
                        fxp2ld(-fxp_largest - fxp_p5));
        test_fxp("Largest + 0.5",
                        fxp_add(fxp_largest, fxp_p5),
                        FXP_PINF_LD);
        test_fxp("-Largest - 0.5",
                        fxp_add(-fxp_largest, -fxp_p5),
                        FXP_NINF_LD);
        test_fxp("+HalfMax + HMaxp2",
                        fxp_add(fxp_halfmax, fxp_halfp2),
                        FXP_PINF_LD);
        test_fxp("-HalfMax - HMaxp2",
                        fxp_add(-fxp_halfmax, -fxp_halfp2),
                        FXP_NINF_LD);
        test_fxp("HalfMax + HalfMax",
                        fxp_add(fxp_halfmax, fxp_halfmax),
                        fxp2ld(fxp_halfmax + fxp_halfmax));
        test_fxp("FXP_MAX - HalfMax",
                        fxp_sub(FXP_MAX, fxp_halfmax),
                        fxp2ld(FXP_MAX - fxp_halfmax));
        test_fxp("HalfMax + FXP_MAX",
                        fxp_add(fxp_halfmax, FXP_MAX),
                        FXP_PINF_LD);
        test_fxp("-FXP_MAX - HalfMax",
                        fxp_sub(-FXP_MAX, fxp_halfmax),
                        FXP_NINF_LD);
        test_fxp("HalfMax * 2",
                        fxp_mul(fxp_halfmax, fxp(2)),
                        fxp2ld(fxp_mul(fxp_halfmax, fxp(2))));
        test_fxp("HalfMax * 2 (long)",
                        fxp_mul_l(fxp_halfmax, fxp(2)),
                        fxp2ld(fxp_mul_l(fxp_halfmax, fxp(2))));
        test_fxp("HalfMax * 3",
                        fxp_mul(fxp_halfmax, fxp(3)),
                        FXP_PINF_LD);
        test_fxp("-HalfMax * 3",
                        fxp_mul(-fxp_halfmax, fxp(3)),
                        FXP_NINF_LD);
        test_fxp("(HalfMax+0.5)*2",
                        fxp_mul(fxp_add(fxp_halfmax, fxp_p5), fxp(2)),
                        FXP_PINF_LD);
        test_fxp("(HalfMax+0.5)*2 (long)",
                        fxp_mul_l(fxp_add(fxp_halfmax, fxp_p5), fxp(2)),
                        FXP_PINF_LD);
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

        test_fxp("-0.(+)500", fxp_dec(-0, 500), 0.5);
        test_fxp("-0.(-)500", fxp_dec(-0, -500), -0.5);

        printf("\nChecking truncation of longer frac decimal arguments, ");
        printf("frac bits: %d\n", frac_bits);
        test_fxp("0.22222", fxp_dec(0, 22222), 0.222);
        test_fxp("0.4444444", fxp_dec(0, 4444444), 0.444);
        test_fxp("0.771999", fxp_dec(0, 771999), 0.771);
        test_fxp("0.9999999", fxp_dec(0, 9999999), 0.999);
        // Restoring original frac_max_dec and delta vars
        fxp_set_frac_max_dec(fmdec);
        warn_delta = bkp_wd;
        larger_delta = bkp_ld;
        max_warn_delta = bkp_mwd;

        printf("\nChecking extreme frac values, ");
        printf("frac bits: %d\n", frac_bits);
        test_fxp("maxfrac", frac_max, fxp2ld(frac_max));
        test_fxp("0.5 + maxfrac",
                fxp_add(fxp_p5, frac_max),
                whole_bits > 1? 0.5 + fxp2ld(frac_max):
                fxp2ld(FXP_POS_INF));
        test_fxp("maxfrac +tiniest",
                fxp_add(frac_max, FXP_TINIEST),
                fxp2ld(fxp_one));
        test_fxp("-maxfrac -tiniest",
                fxp_add(-frac_max, -FXP_TINIEST),
                fxp2ld(-fxp_one));
        test_fxp("maxfrac - maxfrac",
                fxp_add(frac_max, -frac_max),
                0.0);
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
                    fxp_add(fxp_one, fxp_one), fxp2ld(fxp2));
        test_fxp("-1 - 1",
                    fxp_add(-fxp_one, -fxp_one), fxp2ld(-fxp2));
        test_fxp("Ok sum == 2",
                    fxp_add(-fxp_halfmax, fxp_halfp2),
                    fxp2ld(fxp_add(-fxp_halfmax, fxp_halfp2)));
        test_fxp("Ok sum == -2",
                    fxp_add(fxp_halfmax, -fxp_halfp2),
                    fxp2ld(fxp_add(fxp_halfmax, -fxp_halfp2)));
        test_fxp(" num", fxp1,  dnum);
        test_fxp(" num +  2",
                    fxp_add(fxp1, fxp2), dnum + 2.0);
        test_fxp(" num + -2",
                    fxp_add(fxp1, -fxp2), dnum - 2.0);
        test_fxp("-num +  2",
                    fxp_add(-fxp1, fxp2), -dnum + 2.0);
        test_fxp("-num + -2",
                    fxp_add(-fxp1, -fxp2), -dnum - 2.0);
        test_fxp(" num -  2",
                    fxp_sub(fxp1, fxp2), dnum - 2.0);
        test_fxp(" num - -2",
                    fxp_sub(fxp1, -fxp2), dnum + 2.0);
        test_fxp("-num -  2",
                    fxp_sub(-fxp1, fxp2), -dnum - 2.0);
        test_fxp("-num - -2",
                    fxp_sub(-fxp1, -fxp2), -dnum + 2.0);
        test_fxp(" num *  2",
                    fxp_mul(fxp1, fxp2), dnum * 2.0);
        test_fxp(" num * -2",
                    fxp_mul(fxp1, -fxp2), dnum * -2.0);
        test_fxp("-num *  2",
                    fxp_mul(-fxp1, fxp2), -dnum * 2.0);
        test_fxp("-num * -2",
                    fxp_mul(-fxp1, -fxp2), -dnum * -2.0);
        test_fxp(" num *  2 (long)",
                    fxp_mul_l( fxp1, fxp2), dnum * 2.0);
        test_fxp(" num * -2 (long)",
                    fxp_mul_l( fxp1, -fxp2), dnum * -2.0);
        test_fxp("-num *  2 (long)",
                    fxp_mul_l(-fxp1, fxp2), -dnum * 2.0);
        test_fxp("-num * -2 (long)",
                    fxp_mul_l(-fxp1, -fxp2), -dnum * -2.0);
        test_fxp(" num /  2",
                    fxp_div( fxp1, fxp2), dnum / 2.0);
        test_fxp(" num / -2",
                    fxp_div(fxp1, -fxp2), dnum / -2.0);
        test_fxp("-num /  2",
                    fxp_div(-fxp1, fxp2), -dnum / 2.0);
        test_fxp("-num / -2",
                    fxp_div(-fxp1, -fxp2), -dnum / -2.0);
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
                printf("x : "); print_fxp(x); printf("\n");
                printf("f : %fE\n", fx);
                printf("x2: "); print_fxp(x2); printf("\n");

                test_fxp("\nx: ", x2, get_f_target(x));
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
                        test_fxp("\nmul  (x*y)", n1, tgt1);

                        //For division
                        tgt1 = get_div_target(x, y);
                        n1 = fxp_div(x, y);
                        test_fxp("div  (x/y)", n1, tgt1);
                }
        }
}

void test_lg()
{
        printf("\nTesting Transcendental constants as fxp's: ");
        printf("frac bits: %d\n", frac_bits);
        printf("e      : "); print_fxp(fxp_get_e()); printf("\n");
        printf("pi     : "); print_fxp(fxp_get_pi()); printf("\n");
        printf("ln(2)  : "); print_fxp(fxp_get_ln_2()); printf("\n");
        printf("lg10(2): "); print_fxp(fxp_get_lg10_2()); printf("\n\n");

        int x;
        if ((whole_bits >= 3) && (TEST_LG2_MUL_L)) {

                // To test logarithm calculation with fxp_log2_mul_l
                // we need an fxp configuration with at least 3 whole bits
                test_fxp("lg2_mul_l(+INF):",
                        fxp_lg2_mul_l(FXP_POS_INF),
                        get_lg_target(FXP_POS_INF, '2'));
                test_fxp("lg2_mul_l(largest):",
                        fxp_lg2_mul_l(FXP_MAX),
                        get_lg_target(FXP_MAX, '2'));
                x = fxp_bin(2, FXP_frac_mask/5);
                test_fxp("lg2_mul_l(2.2):",
                        fxp_lg2_mul_l(x),
                        get_lg_target(x, '2'));
                test_fxp("lg2_mul_l(2):",
                        fxp_lg2_mul_l(fxp(2)),
                        get_lg_target(fxp(2), '2'));
                x = fxp_bin(1, FXP_frac_max);
                test_fxp("lg2_mul_l(1.99...9):",
                        fxp_lg2_mul_l(x), get_lg_target(x, '2'));
                x = fxp_bin(1, 1);
                test_fxp("lg2_mul_l(1.00..01):",
                        fxp_lg2_mul_l(x), get_lg_target(x, '2'));
                x = fxp(1);
                test_fxp("lg2_mul_l(1):",
                        fxp_lg2_mul_l(x), get_lg_target(x, '2'));
                x = fxp(1) - 1;
                test_fxp("lg2_mul_l(0.99...):",
                        fxp_lg2_mul_l(x), get_lg_target(x, '2'));
                test_fxp("lg2_mul_l(0):",
                        fxp_lg2_mul_l(0), get_lg_target(0, '2'));
                test_fxp("lg2_mul_l(-INF):",
                        fxp_lg2_mul_l(FXP_NEG_INF),
                        get_lg_target(FXP_NEG_INF, '2'));
                test_fxp("lg2_mul_l(UNDEF):",
                        fxp_lg2_mul_l(FXP_UNDEF),
                        get_lg_target(FXP_UNDEF, '2'));
        } else {
                printf("Skipping lg2_mul_l() tests ");
                if (TEST_LG2_MUL_L) {
                        printf("(Only %d whole bit(s); 3 or more needed)\n", \
                                    whole_bits);
                } else {
                        printf("\n");
                }
        }

        test_fxp("\nlg2_l(+INF):",
                fxp_lg2_l(FXP_POS_INF),
                get_lg_target(FXP_POS_INF, '2'));
        test_fxp("lg2_l(largest):",
                fxp_lg2_l(FXP_MAX),
                get_lg_target(FXP_MAX, '2'));
        x = fxp_bin(2, FXP_frac_mask/5);
        test_fxp("lg2_l(2.2):",
                fxp_lg2_l(x),
                get_lg_target(x, '2'));
        test_fxp("lg2_l(2):",
                fxp_lg2_l(fxp(2)),
                get_lg_target(fxp(2), '2'));
        x = fxp_bin(1, FXP_frac_max);
        test_fxp("lg2_l(1.99...9):",
                fxp_lg2_l(x), get_lg_target(x, '2'));
        x = fxp_bin(1, 1);
        test_fxp("lg2_l(1.00..01):",
                fxp_lg2_l(x), get_lg_target(x, '2'));
        x = fxp(1);
        test_fxp("lg2_l(1):",
                fxp_lg2_l(x), get_lg_target(x, '2'));
        x = fxp(1) - 1;
        test_fxp("lg2_l(0.99...):",
                fxp_lg2_l(x), get_lg_target(x, '2'));
        test_fxp("lg2_l(0):",
                fxp_lg2_l(0), get_lg_target(0, '2'));
        test_fxp("lg2_l(-INF):",
                fxp_lg2_l(FXP_NEG_INF),
                get_lg_target(FXP_NEG_INF, '2'));
        test_fxp("lg2_l(UNDEF):",
                fxp_lg2_l(FXP_UNDEF),
                get_lg_target(FXP_UNDEF, '2'));

        printf("\nlg2(+INF)       : ");
                print_fxp(fxp_lg2(FXP_POS_INF)); printf("\n");
        printf("lg2(2)          : ");
                print_fxp(fxp_lg2(fxp(2))); printf("\n");
        printf("lg2(1)          : ");
                print_fxp(fxp_lg2(fxp(1))); printf("\n");
        printf("lg2(0.99...)    : ");
                print_fxp(fxp_lg2(x)); printf("\n");
        printf("lg2(0)          : ");
                print_fxp(fxp_lg2(0)); printf("\n");
        printf("lg2(-INF)       : ");
                print_fxp(fxp_lg2(FXP_NEG_INF)); printf("\n");
        printf("lg2(UNDEF)      : ");
                print_fxp(fxp_lg2(FXP_UNDEF)); printf("\n\n");

/*
        printf("ln(0)           : ");
                print_fxp(fxp_ln(0)); printf("\n");
        printf("ln(1)           : ");
                print_fxp(fxp_ln(fxp(1))); printf("\n");
        printf("ln(e)           : ");
                print_fxp(fxp_ln(fxp_get_e())); printf("\n\n");
        printf("lg10(0)         : ");
                print_fxp(fxp_lg10(0)); printf("\n");
        printf("lg10(1)         : ");
                print_fxp(fxp_lg10(fxp(1))); printf("\n");
        printf("lg10(10)        : ");
                print_fxp(fxp_lg10(fxp(10))); printf("\n");
*/
}

void test_pow()
{
        int x = 0;
        test_fxp("\npow2_l(0):", fxp_pow2_l(x),
                    get_target(powl(2.0, fxp2ld(x))));
        x = fxp_bin(0, FXP_frac_mask / 2);
        test_fxp("pow2_l(0.5):", fxp_pow2_l(x),
                    get_target(powl(2.0, fxp2ld(x))));
        x = fxp(1);
        test_fxp("pow2_l(1):", fxp_pow2_l(x),
                    get_target(powl(2.0, fxp2ld(x))));
        x = fxp(2);
        test_fxp("pow2_l(2):", fxp_pow2_l(x),
                    get_target(powl(2.0, fxp2ld(x))));
        x = fxp_bin(3, FXP_frac_max / 2);
        test_fxp("pow2_l(3.5):", fxp_pow2_l(x),
                    get_target(powl(2.0, fxp2ld(x))));
        x = fxp_bin(14, FXP_frac_max);
        test_fxp("pow2_l(14.999...):", fxp_pow2_l(x),
                    get_target(powl(2.0, fxp2ld(x))));

        printf("pow2_l(-0.999..): ");
                print_fxp(fxp_pow2_l(fxp_bin(0, FXP_frac_mask)));
                printf("\n");
        printf("pow2_l(-1)      : ");
                print_fxp(fxp_pow2_l(fxp(-1))); printf("\n");
        printf("pow2_l(-1.0..01): ");
                print_fxp(fxp_pow2_l(fxp_bin(-1, 1))); printf("\n");
        printf("pow2_l(-2)      : ");
                print_fxp(fxp_pow2_l(fxp(-2))); printf("\n");
        printf("pow2_l(-3.5)    : ");
                print_fxp(fxp_pow2_l(d2fxp(-3.5))); printf("\n");
        printf("pow2_l(most neg): ");
                print_fxp(fxp_pow2_l(-fxp_largest)); printf("\n");
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
                //printf("n1 = "); print_fxp(n1);
                //printf("\nn2 = "); print_fxp(n2);
                //printf("\nn3 = "); print_fxp(n3);
                //printf(" (n3 always in (-1,1))");
                //if (n3 == 0) printf(" == 0");
                //printf("\n");

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
                test_fxp("n1: ", f2fxp(fxp2f(n1)), get_f_target(n1));
                test_fxp("n2: ", f2fxp(fxp2f(n2)), get_f_target(n2));
                test_fxp("n3: ", f2fxp(fxp2f(n3)), get_f_target(n3));

                fxp1 = fxp_add(n1, n2);
                tgt1 = get_target(ldx + ldy);
                test_fxp("add   (n1+n2)", fxp1, tgt1);
                fxp1 = fxp_add(n1, n3);
                tgt2 = get_target(ldx + ldz);
                test_fxp("add   (n1+n3)", fxp1, tgt2);

                fxp1 = fxp_mul(n1, n2);
                fxp2 = fxp_mul_l(n1, n2);
                tgt1 = get_target(ldx * ldy);
                test_fxp("mul   (n1*n2)", fxp1, tgt1);
                test_fxp("mul_l (n1*n2)", fxp2, tgt1);
                fxp1 = fxp_mul(n1, n3);
                fxp2 = fxp_mul_l(n1, n3);
                tgt2 = get_target(ldx * ldz);
                test_fxp("mul   (n1*n3)", fxp1, tgt2);
                test_fxp("mul_l (n1*n3)", fxp2, tgt2);

                fxp1 = fxp_div_l(n1, n2);
                fxp2 = fxp_div(n1, n2);
                tgt1 = get_div_target(n1, n2);
                test_fxp("div_l (n1/n2)", fxp1, tgt1);
                test_fxp("div   (n1/n2)", fxp2, tgt1);
                fxp1 = fxp_div_l(n1, n3);
                fxp2 = fxp_div(n1, n3);
                tgt2 = get_div_target(n1, n3);
                test_fxp("div_l (n1/n3)", fxp1, tgt2);
                test_fxp("div   (n1/n3)", fxp2, tgt2);

                // Make arguments positive for the log tests
                n1 = n1 < 0? -n1: n1;
                n2 = n2 < 0? -n2: n2;
                n3 = n3 < 0? -n3: n3;
                printf("\nTesting lg2 implementations, ");
                printf("frac bits: %d", frac_bits);
                printf("\nn1 = "); print_fxp(n1);
                printf("\nn2 = "); print_fxp(n2);
                printf("\nn3 = "); print_fxp(n3); printf("\n");

                tgt1 = get_lg_target(n1, '2');
                tgt2 = get_lg_target(n2, '2');
                tgt3 = get_lg_target(n3, '2');
                if ((whole_bits >= 3) && (TEST_LG2_MUL_L)) {
                        test_fxp("lg2_mul_l(n1)", fxp_lg2_mul_l(n1), tgt1);
                        test_fxp("lg2_mul_l(n2)", fxp_lg2_mul_l(n2), tgt2);
                        test_fxp("lg2_mul_l(n3)", fxp_lg2_mul_l(n3), tgt3);
                }
                test_fxp("lg2(n1)", fxp_lg2(n1), tgt1);
                test_fxp("lg2(n2)", fxp_lg2(n2), tgt2);
                test_fxp("lg2(n3)", fxp_lg2(n3), tgt3);

                test_fxp("lg2_l(n1)", fxp_lg2_l(n1), tgt1);
                test_fxp("lg2_l(n2)", fxp_lg2_l(n2), tgt2);
                test_fxp("lg2_l(n3)", fxp_lg2_l(n3), tgt3);

                // Tests for other logarithms and pow still pending,
                // will enable after all-bits frac multiplication
                // (regardless of fxp configuration) is tested
                continue;

                tgt1 = get_lg_target(n1, 'e');
                tgt2 = get_lg_target(n2, 'e');
                tgt3 = get_lg_target(n3, 'e');
                test_fxp("ln(n1)", fxp_ln(n1), tgt1);
                test_fxp("ln(n2)", fxp_ln(n2), tgt2);
                test_fxp("ln(n3)", fxp_ln(n3), tgt3);

                test_fxp("ln_l(n1)", fxp_ln_l(n1), tgt1);
                test_fxp("ln_l(n2)", fxp_ln_l(n2), tgt2);
                test_fxp("ln_l(n3)", fxp_ln_l(n3), tgt3);

                tgt1 = get_lg_target(n1, 'a');
                tgt2 = get_lg_target(n2, 'a');
                tgt3 = get_lg_target(n3, 'a');
                test_fxp("lg10(n1)", fxp_lg10(n1), tgt1);
                test_fxp("lg10(n2)", fxp_lg10(n2), tgt2);
                test_fxp("lg10(n3)", fxp_lg10(n3), tgt3);

                test_fxp("lg10_l(n1)", fxp_lg10_l(n1), tgt1);
                test_fxp("lg10_l(n2)", fxp_lg10_l(n2), tgt2);
                test_fxp("lg10_l(n3)", fxp_lg10_l(n3), tgt3);
                /*
                tgt2 = get_target(powl(2, fxp2ld(n3)));
                tgt3 = get_target(powl(2, fxp2ld(-n3)));
                test_fxp("pow2_l( n3)", fxp_pow2_l(n3), tgt2);
                test_fxp("pow2_l(-n3)", fxp_pow2_l(-n3), tgt3);
                */
        }
}

void test_mul_wrefs()
{
        printf("\nTesting mul_wrefs\n");
        int n = fxp_bin(0, FXP_frac_max);
        int m = fxp_bin(0, 1 << (FXP_frac_bits - 1));
        //unsigned int f1 = ((unsigned int) fxp_get_bin_frac(n)) << FXP_WORD_BITS;
        //unsigned int f2 = ((unsigned int) fxp_get_bin_frac(m)) << FXP_WORD_BITS;
        unsigned int f1 = FXP_LWORD_MASK | FXP_RWORD_MASK;
        unsigned int f2 = f1;
        unsigned int r;
        //r = mul_fracs(f1, f2);
        //r = mul_fracs(16, 16);
        //printf("0.%X * 0.%X == 0.%X\n", f1, f2, r);
        n = fxp_dec(1133, 6699);
        m = fxp_dec(  24, 5577);
        printf("n: "); print_fxp(n); printf("\n");
        printf("m: "); print_fxp(m); printf("\n");
        int rw = 0, rf = 0;
        int prod = fxp_mul_wrefs(n, m, &rw, &rf);

        int p = fxp_mul(n, m);
        test_fxp("1133.6699 * 24.5577", p, get_mul_target(n, m));
}

int main(void)
{
        printf("\n%sFXP Tester run\n%s", DASHES, DASHES);
        print_sys_info();

        /*
        printf("\nTesting auto frac max dec:");
        for (int nfracbits = -1 ; nfracbits <= 33; nfracbits++) {
                fxp_set_frac_bits(nfracbits);
                fxp_set_auto_frac_max_dec();
                printf("\nAttempted frac bits: %d, effective: %d, ", \
                    nfracbits,
                    fxp_get_frac_bits());
                printf(" auto frac max dec: %d", fxp_get_frac_max_dec());
        }
        printf("\n");
        */

        /*
        printf("\nTesting function that counts # bits used by a number:");
        unsigned int n = 1;
        int nb = 1;
        do {
                printf("\nUnsigned number %u uses %d bits", n-1, \
                        fxp_nbits(n-1));
                printf("\nUnsigned Number %u uses %d bits", n, \
                        fxp_nbits(n));
                assert(fxp_nbits(n) == fxp_nbits_v0(n, FXP_INT_BITS));
                n = (n << 1);
                nb++;
        } while (nb <= FXP_INT_BITS);
        printf("\nUnsigned Number %u uses %d bits\n", UINT_MAX, \
                fxp_nbits(UINT_MAX));
        */

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
                printf("largest <<: %lX  (largest l-shifted and rounded)\n", \
                                FXP_max_lshifted);
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
                printf("FXP_max_f : %E\n", FXP_max_f);
                printf("FXP_min_f : %E\n", FXP_min_f);
                printf("FXP_max_d : %lE\n", FXP_max_d);
                printf("FXP_min_d : %lE\n", FXP_min_d);
                printf("FXP_max_ld: %LE\n", FXP_max_ld);
                printf("FXP_min_ld: %LE\n", FXP_min_ld);

                // Tests to run =============================
                tests_01();
                tests_02();
                tests_03();
                test_decbin_mappings();
                test_fracs();
                test_ops_with_whole_bits();
                //test_mul_wrefs();
                test_ops_with_values_of_interest();
                test_lg();
                //test_pow();
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
