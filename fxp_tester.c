/*
 * fxp_tester.c
 * Tests the implementation of binary fixed point numbers
 * and arithmetic operations (fxp.c)
 * By Raul Saavedra, 2022-11-13, Bonn Germany
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <float.h>
#include <time.h>
#include <math.h>
#include "fxp.h"
#include "fxp_aux.h"

#define DASHES "========================\n"

static int fracbit_configs[] = {0, 11, 14, 16, 24, 31};
//static int fracbit_configs[] = {4};
//static int fracbit_configs[] = {31};
static long double max_warn_delta = 0.0;
static long double larger_delta = 0.0;
static long double largest_delta = 0.0;
static long double largest_madelta = 0.0;
static int largest_delta_fbits = 0;
static long double warn_delta = 0.0;
static int twarnings = 0;
static int nwarnings = 0;

void test_fxp(char *s, int fxp1, long double d_assert_val)
{
        printf("%s\n", s);
        printf("\tgot: "); print_fxp(fxp1);
        printf("\n\texp. ");
        long double avplim = lim_frac(d_assert_val, fxp_get_frac_bits());
        if (d_assert_val == dfxp_pos_inf())
            printf("+INF");
        else if (d_assert_val == dfxp_neg_inf())
            printf("-INF");
        else if (d_assert_val == dfxp_undef())
            printf("UNDEF");
        else
            printf("%Lf", avplim);
        if (((fxp1 == FXP_UNDEF) && (d_assert_val == dfxp_undef())) \
            || ((fxp1 == FXP_POS_INF) && (d_assert_val == dfxp_pos_inf())) \
            || ((fxp1 == FXP_NEG_INF) && (d_assert_val == dfxp_neg_inf()))) {
            printf(" (~same)\n");
            return;
        }

        long avlimw = (long) avplim;
        long double df = dfxp(fxp1);
        long double delta = (df >= avplim)?
                        df - avplim:
                        avplim - df;

        if (delta <= warn_delta) {
            // No warning up to the warn_delta value
            printf(" (~same)\n");
        } else {
            nwarnings++;
            printf("\n***** Warning %d: delta=%LE (above:%LE, max allowed %LE)\n\n",
                        nwarnings, delta, warn_delta, max_warn_delta);
            if (delta > larger_delta) {
                larger_delta = delta;
                printf("Updating largest delta: %Lf\n", larger_delta);
            }
        }
        fflush( stdout );

        // Here assert that we never exceed the max_warn_delta
        assert( (delta <= max_warn_delta) &&
            ((df >= 0 && avplim >= 0) ||
                (df < 0 && avplim < 0)));
}

static long double get_target(long double x)
{
    long double y = lim_frac(x, fxp_get_frac_bits());
    if (y >= dfxp_pos_inf())
        return dfxp_pos_inf();
    else if (y <= dfxp_neg_inf())
        return dfxp_neg_inf();
    return y;
}

static long double get_div_target(int x, int y)
{
    long double target;
    if ((x == FXP_UNDEF) || y == (FXP_UNDEF) \
        || ((x == 0) && (y == 0)) \
        || ((x == FXP_POS_INF) && (y == FXP_POS_INF)) \
        || ((x == FXP_POS_INF) && (y == FXP_NEG_INF)) \
        || ((x == FXP_NEG_INF) && (y == FXP_POS_INF)) \
        || ((x == FXP_NEG_INF) && (y == FXP_NEG_INF))) {
        target = dfxp_undef();
    } else {
        if (((x > 0) && (y == 0)) \
            || ((x == FXP_POS_INF) && (y > 0)) \
            || ((x == FXP_NEG_INF) && (y < 0))) {
            target = dfxp_pos_inf();
        } else {
            if (((x < 0) && (y == 0)) \
                || ((x == FXP_NEG_INF) && (y > 0)) \
                || ((x == FXP_POS_INF) && (y < 0))) {
                target = dfxp_neg_inf();
            } else {
                if ((y == FXP_POS_INF) || (y == FXP_NEG_INF)) {
                    target = 0;
                } else {
                    long double ldx = dfxp(x);
                    long double ldy = dfxp(y);
                    int frbits = fxp_get_frac_bits();
                    target = lim_frac(ldx, frbits) / lim_frac(ldy, frbits);
                    target = lim_frac(target, frbits);
                    if (target >= dfxp_pos_inf()) {
                        target = dfxp_pos_inf();
                    } else {
                        if (target <= dfxp_neg_inf()) {
                            target = dfxp_neg_inf();
                        }
                    }
                }
            }
        }
    }
    return target;
}



int main(void) {
        printf("%sFXP Tester run\n%s", DASHES, DASHES);

        printf("Num type sizes in this system:\n");
        printf("char        has a size of %zd bytes.\n", sizeof(char));
        printf("int         has a size of %zd bytes.\n", sizeof(int));
        printf("long        has a size of %zd bytes.\n", sizeof(long));
        printf("long long   has a size of %zd bytes.\n", sizeof(long long));
        printf("float       has a size of %zd bytes.\n", sizeof(float));
        printf("double      has a size of %zd bytes.\n", sizeof(double));
        printf("long double has a size of %zd bytes.\n", sizeof(long double));

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

        // Randomize seed
        // (comment this out to be able to reproduce exactly between runs)
        //srand((unsigned int) time(0));

        int nconfigs = sizeof(fracbit_configs) / sizeof(fracbit_configs[0]);
        printf("\nRunning tests for frac-bit sizes: ");
        for (int nfb = 0; nfb < nconfigs; nfb++) {
            printf("%d", fracbit_configs[nfb]);
            if ((nfb + 1) < nconfigs) printf(", ");
        }
        printf("\n");

        for (int nfb = 0; nfb < nconfigs; nfb++) {
                fxp_set_frac_bits(fracbit_configs[nfb]);
                fxp_set_auto_frac_max_dec();

                int frac_bits = fxp_get_frac_bits();
                int whole_bits = fxp_get_whole_bits();
                int frac_mask = fxp_get_frac_mask();
                int frac_max = fxp_get_frac_max();
                int frac_max_dec = fxp_get_frac_max_dec();
                int whole_max = fxp_get_whole_max();
                int whole_min = fxp_get_whole_min();

                // Because of the precision loss in fxp_mul, let's
                // allow errors up to half the frac bits in use, but
                // start warning about calculation differences from target
                // when the error exceeds half the frac_bits
                warn_delta = ((long double) (1 << MAX(1, (frac_bits / 2)))) \
                                             / frac_max;
                warn_delta = lim_frac(warn_delta, frac_bits);
                // And assert/halt the tester if the error ever exceeds
                // 3/4 of the frac bits in use
                max_warn_delta = lim_frac(warn_delta * 1.5, frac_bits);

                nwarnings = 0;
                larger_delta = 0;
                printf("\n%sFXP configuration parameters:\n", DASHES);
                printf("frac bits   : %d (requested was %d)\n",

                                frac_bits, fracbit_configs[nfb]);
                printf("whole bits  : %d\n", whole_bits);
                printf("pos infinity: %d (Lf: %Lf)\n",
                                FXP_POS_INF, dfxp_pos_inf());
                printf("whole max   : %d\n", whole_max);
                printf("frac mask   : %d\n", frac_mask);
                printf("frac max    : %d (->decimals: .%d)\n",
                                frac_max,
                                frac_max_dec);
                printf("whole min   : %d\n", whole_min);
                printf("neg infinity: %d (Lf: %Lf)\n",
                                FXP_NEG_INF, dfxp_neg_inf());
                printf("undefined   : %d (Lf: %Lf)\n",
                                FXP_UNDEF, dfxp_undef());
                printf("warn_delta   : %LE\n", warn_delta);
                printf("max_warn_delta: %LE\n", max_warn_delta);

                /*
                printf("\nChecking conversion between fxp and lim-prec long doubles:");
                test_fxp("fraction", 0, int_to_frac(0));
                int max = fxp_get_frac_max() >> (frac_bits < FXP_INT_BITS_M1? 0: 1);
                for (int frac = 1; frac <= max; frac = frac << 1) {
                    int fxpfr = fxp_bin(0, frac);
                    test_fxp("fraction", fxpfr, int_to_frac(frac));
                }
                */

                printf("\nChecking extreme int values, part I:\n");
                int fxp_largest = FXP_MAX;
                int fxp_tiniest = 1;
                int fxp_ten = fxp_bin(10, 0);
                int fxp_two = fxp_bin(2, 0);
                int fxp_one = fxp_bin(1, 0);
                int fxp_p5 = fxp_bin(0, frac_max/2 + 1);
                int fxp_halfmax = FXP_MAX / 2;
                int fxp_halfp2 = fxp_add(fxp_halfmax, fxp_two);
                test_fxp("Infinity", FXP_POS_INF, dfxp_pos_inf());
                test_fxp("Largest", fxp_largest, dfxp(FXP_MAX));
                test_fxp("HalfMax", fxp_halfmax, dfxp(FXP_MAX)/2);
                test_fxp("Largest frac", frac_max, dfxp(fxp_bin(0, frac_max)));
                test_fxp("tiniest", fxp_tiniest, dfxp(fxp_tiniest));
                test_fxp("0.5", fxp_p5, ((long double) 0.5));
                test_fxp("zero", fxp_bin(0, 0), 0.0);
                test_fxp("-tiniest", -fxp_tiniest, dfxp(-fxp_tiniest));
                test_fxp("-Largest frac", -frac_max, dfxp(-frac_max));
                test_fxp("Most negative", -fxp_largest, dfxp(FXP_MIN));
                test_fxp("-Infinity", FXP_NEG_INF, dfxp_neg_inf());
                test_fxp("Undefined", FXP_UNDEF, dfxp_undef());


                printf("\nChecking extreme int values, part II:\n");
                test_fxp("Almost most negative",
                                fxp_add(-fxp_largest, fxp_tiniest),
                                dfxp(-fxp_largest) + dfxp(fxp_tiniest));
                test_fxp(" Largest -Largest",
                                fxp_add(fxp_largest, -fxp_largest),
                                0.0);
                test_fxp("-Largest +Largest",
                                fxp_add(-fxp_largest, fxp_largest),
                                0.0);
                test_fxp("Largest + 0",
                                fxp_add(fxp_largest, fxp(0)),
                                dfxp(fxp_largest));
                test_fxp("-Largest - 0",
                                fxp_add(-fxp_largest, -fxp(0)),
                                dfxp(-fxp_largest));
                test_fxp("Largest - tiniest",
                                fxp_sub(fxp_largest, fxp_tiniest),
                                dfxp(fxp_largest) - dfxp(fxp_tiniest));
                test_fxp("Largest + tiniest safe",
                                fxp_add(fxp_largest, fxp_tiniest),
                                dfxp_pos_inf());
                test_fxp("Largest + tiniest unsafe",
                                fxp_unsafe_add(fxp_largest, fxp_tiniest),
                                dfxp_pos_inf());
                test_fxp("+inf + +inf",
                                fxp_add(FXP_POS_INF, FXP_POS_INF),
                                dfxp_pos_inf());
                test_fxp("-inf - +inf",
                                fxp_sub(FXP_NEG_INF, FXP_POS_INF),
                                dfxp_neg_inf());
                test_fxp("+inf + -inf",
                                fxp_add(FXP_POS_INF, FXP_NEG_INF),
                                dfxp_undef());
                test_fxp("-inf + -inf",
                                fxp_add(FXP_NEG_INF, FXP_NEG_INF),
                                dfxp_neg_inf());
                test_fxp("-inf - -inf",
                                fxp_sub(FXP_NEG_INF, FXP_NEG_INF),
                                dfxp_undef());
                test_fxp("+inf * -inf",
                                fxp_mul(FXP_POS_INF, FXP_NEG_INF),
                                dfxp_neg_inf());
                test_fxp("+inf - 0.5",
                                fxp_sub(FXP_POS_INF, fxp_p5),
                                dfxp_pos_inf());
                test_fxp("-inf + 0.5",
                                fxp_add(FXP_NEG_INF, fxp_p5),
                                dfxp_neg_inf());
                test_fxp("+num / zero",
                                fxp_div(fxp_largest, fxp(0)),
                                dfxp_pos_inf());
                test_fxp("zero / zero",
                                fxp_div(fxp(0), fxp(0)),
                                dfxp_undef());
                test_fxp("zero * zero",
                                fxp_mul(fxp(0), fxp(0)),
                                0.0);
                test_fxp("zero + zero",
                                fxp_add(fxp(0), fxp(0)),
                                0.0);
                test_fxp("zero - zero",
                                fxp_sub(fxp(0), fxp(0)),
                                0.0);
                test_fxp("-num / zero",
                                fxp_div(-fxp_largest, fxp(0)),
                                dfxp_neg_inf());
                test_fxp("zero * +inf",
                                fxp_mul(fxp(0), FXP_POS_INF),
                                dfxp_undef());
                test_fxp("zero * -inf",
                                fxp_mul(fxp(0), FXP_NEG_INF),
                                dfxp_undef());
                test_fxp("zero * undef",
                                fxp_mul(fxp(0), FXP_UNDEF),
                                dfxp_undef());
                test_fxp("-inf * undef",
                                fxp_mul(FXP_NEG_INF, FXP_UNDEF),
                                dfxp_undef());
                test_fxp("+inf * undef",
                                fxp_mul(FXP_POS_INF, FXP_UNDEF),
                                dfxp_undef());
                test_fxp("undef * undef",
                                fxp_mul(FXP_UNDEF, FXP_UNDEF),
                                dfxp_undef());
                test_fxp("tiniest * inf",
                                fxp_mul(fxp_tiniest, FXP_POS_INF),
                                dfxp_pos_inf());


                printf("\nChecking extreme int values, part III:\n");
                test_fxp("Way Too Large whole part!",
                                fxp_add(fxp(whole_max), fxp(5)),
                                dfxp_pos_inf());
                test_fxp("Largest * 1",
                                fxp_mul(fxp_largest,  fxp_one),
                                dfxp(fxp_largest));
                test_fxp("Largest * -1",
                                fxp_mul(fxp_largest,  -fxp_one),
                                dfxp(-fxp_largest));
                test_fxp("Largest + two safe",
                                fxp_add(fxp_largest, fxp_two),
                                dfxp_pos_inf());
                test_fxp("Largest + two unsafe",
                                fxp_unsafe_add(fxp_largest, fxp_two),
                                dfxp(fxp_largest + fxp_two));
                test_fxp("Safe Too neg substraction",
                                fxp_sub( -fxp_largest, fxp_one),
                                dfxp_neg_inf());
                test_fxp("Unsafe Too neg substraction",
                                fxp_unsafe_sub( -fxp_largest, fxp_p5),
                                dfxp(-fxp_largest - fxp_p5));
                test_fxp("Largest + 0.5",
                                fxp_add(fxp_largest, fxp_p5),
                                dfxp_pos_inf());
                test_fxp("-Largest - 0.5",
                                fxp_add(-fxp_largest, -fxp_p5),
                                dfxp_neg_inf());
                test_fxp("+HalfMax + HMaxp2",
                                fxp_add(fxp_halfmax, fxp_halfp2),
                                dfxp_pos_inf());
                test_fxp("-HalfMax - HMaxp2",
                                fxp_add(-fxp_halfmax, -fxp_halfp2),
                                dfxp_neg_inf());
                test_fxp("HalfMax + HalfMax",
                                fxp_add(fxp_halfmax, fxp_halfmax),
                                dfxp(fxp_halfmax + fxp_halfmax));
                test_fxp("FXP_MAX - HalfMax",
                                fxp_sub(FXP_MAX, fxp_halfmax),
                                dfxp(FXP_MAX - fxp_halfmax));
                test_fxp("HalfMax + FXP_MAX",
                                fxp_add(fxp_halfmax, FXP_MAX),
                                dfxp_pos_inf());
                test_fxp("-FXP_MAX - HalfMax",
                                fxp_sub(-FXP_MAX, fxp_halfmax),
                                dfxp_neg_inf());
                test_fxp("HalfMax * 2",
                                fxp_mul(fxp_halfmax, fxp(2)),
                                dfxp(FXP_MAX));
                test_fxp("HalfMax * 2 (long)",
                                fxp_mul_l(fxp_halfmax, fxp(2)),
                                dfxp(FXP_MAX));
                test_fxp("HalfMax * 3",
                                fxp_mul(fxp_halfmax, fxp(3)),
                                dfxp_pos_inf());
                test_fxp("-HalfMax * 3",
                                fxp_mul(-fxp_halfmax, fxp(3)),
                                dfxp_neg_inf());
                test_fxp("(HalfMax+0.5)*2",
                                fxp_mul(fxp_add(fxp_halfmax, fxp_p5), fxp(2)),
                                dfxp_pos_inf());
                test_fxp("(HalfMax+0.5)*2 (long)",
                                fxp_mul_l(fxp_add(fxp_halfmax, fxp_p5), fxp(2)),
                                dfxp_pos_inf());

                printf("\nChecking decimal <=> bin mappings of frac ranges:");
                int fmbin = fxp_get_frac_max();
                int fmdec = fxp_get_frac_max_dec();
                printf("\nMax frac dec: %d (bin %d)", fmdec, fmbin);
                for (int i = 0; i <= 5; i++) {
                        printf("\nShowing fxp for 0.%3d: ", i);
                        int vf = fxp_dec(0, i);
                        print_fxp(vf);
                }
                printf("\n:");
                int m = (fxp_get_frac_max_dec() + 1) / 2;
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

                printf("\nChecking sign taken from frac when whole == 0:\n");
                // Temporary switching to a fixed frac_max_dec for these tests,
                // and relaxing the max allowed delta
                long double bkp_wd = warn_delta;
                long double bkp_ld = larger_delta;
                long double bkp_mwd = max_warn_delta;
                fxp_set_frac_max_dec(999);
                warn_delta = lim_frac(dfxp(fxp_dec(0,100)), frac_bits);
                max_warn_delta = lim_frac(warn_delta * 2, frac_bits);

                test_fxp("-0.(+)500", fxp_dec(-0, 500), 0.5);
                test_fxp("-0.(-)500", fxp_dec(-0, -500), -0.5);

                printf("\nTruncation of longer frac decimal arguments:\n");
                test_fxp("0.22222", fxp_dec(0, 22222), 0.222);
                test_fxp("0.4444444", fxp_dec(0, 4444444), 0.444);
                test_fxp("0.771999", fxp_dec(0, 771999), 0.771);
                test_fxp("0.9999999", fxp_dec(0, 9999999), 0.999);
                // Restoring original frac_max_dec and delta vars
                fxp_set_frac_max_dec(fmdec);
                warn_delta = bkp_wd;
                larger_delta = bkp_ld;
                max_warn_delta = bkp_mwd;

                printf("\nChecking extreme frac values:\n");
                test_fxp("maxfrac", frac_max, dfxp(frac_max));
                test_fxp("0.5 + maxfrac",
                        fxp_add(fxp_p5, frac_max),
                        whole_bits > 1? 0.5 + dfxp(frac_max):
                        dfxp(FXP_POS_INF));
                test_fxp("maxfrac +tiniest",
                        fxp_add(frac_max, fxp_tiniest),
                        dfxp(fxp_one));
                test_fxp("-maxfrac -tiniest",
                        fxp_add(-frac_max, -fxp_tiniest),
                        dfxp(-fxp_one));
                test_fxp("maxfrac - maxfrac",
                        fxp_add(frac_max, -frac_max),
                        0.0);


                if (whole_bits >= 3) {
                        printf("\nSimple ops when using 3 or more bits for the whole part:\n");
                        int whole = 0;
                        int bin_frac = frac_max / 2;  // == 0.500
                        int num = fxp_bin(whole, bin_frac);
                        int dec_frac = fxp_get_dec_frac(num);
                        int fxp1 = num;
                        int fxp2 = fxp(2);
                        long double dnum = dfxp(num);
                        test_fxp(" 1 + 1",
                                    fxp_add(fxp_one, fxp_one), dfxp(fxp2));
                        test_fxp("-1 - 1",
                                    fxp_add(-fxp_one, -fxp_one), dfxp(-fxp2));
                        test_fxp("Ok sum == 2",
                                    fxp_add(-fxp_halfmax, fxp_halfp2), dfxp(fxp2));
                        test_fxp("Ok sum == -2",
                                    fxp_add(fxp_halfmax, -fxp_halfp2), dfxp(-fxp2));
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


                int sign1, sign2, sign3, n1, n2, n3, n4, fxp1, fxp2;
                long double ldx, ldy, ldz, tgt1, tgt2;

                printf("\nVerifying divisions with values of interest:\n");
                int dividend[] = {FXP_UNDEF, FXP_POS_INF, FXP_MAX,
                        fxp(300000), fxp_dec(3000, 5), fxp(300), fxp(30),
                        fxp(3), fxp(2), fxp(1),
                        0, fxp(-1), fxp(-2), fxp(-3),
                        fxp(-30), fxp(-300), fxp_dec(-3000, 5), fxp(-300000)};
                int divisor[] = {FXP_UNDEF, FXP_POS_INF, FXP_NEG_INF, 0, fxp(1)};
                int x, y, posx, posy;
                int ndd = (int) (sizeof(dividend) / sizeof(int));
                int ndr = (int) (sizeof(divisor) / sizeof(int));
                for (int i=0; i<ndd; i++) {
                    for (int j=0; j<ndr; j++) {
                        x = dividend[i];
                        y = divisor[j];
                        ldx = dfxp(x);
                        ldy = dfxp(y);
                        tgt1 = get_div_target(x, y);
                        printf("x: ");
                        print_fxp(x);
                        printf("\ny: ");
                        print_fxp(y);
                        n1 = fxp_div(x, y);
                        test_fxp("\ndiv  (x/y)", n1, tgt1);
                    }
                }

                printf("\nVerifying different op ");
                printf("implementations using random numbers:\n");
                for (int i=0; i < 10; i++) {
                        sign1 = rand() % 2 == 1? -1: 1;
                        sign2 = rand() % 2 == 1? -1: 1;
                        sign3 = rand() % 2 == 1? -1: 1;
                        n1 = sign1 * rand();
                        n2 = sign2 * rand();
                        n3 = sign3 * rand();
                        if (frac_bits < FXP_INT_BITS_M1) {
                            // n3 always in (-1, 1)
                            n3 %= fxp_get_frac_max();
                        }
                        ldx = dfxp(n1);
                        ldy = dfxp(n2);
                        ldz = dfxp(n3);

                        printf("\nRandom Test #%d for %d frac bits:\n", i, frac_bits);
                        printf("n1 = "); print_fxp(n1);
                        printf("\nn2 = "); print_fxp(n2);
                        printf("\nn3 = "); print_fxp(n3);
                        printf("( n3 always in (-1,1) )");
                        if (n3 == 0) printf(" == 0");
                        printf("\n");

                        fxp1 = fxp_add(n1, n2);
                        fxp2 = fxp_add_l(n1, n2);
                        tgt1 = get_target(lim_frac(ldx + ldy, frac_bits));
                        test_fxp("add   (n1+n2)", fxp1, tgt1);
                        test_fxp("add_l (n1+n2)", fxp2, tgt1);
                        fxp1 = fxp_add(n1, n3);
                        fxp2 = fxp_add_l(n1, n3);
                        tgt2 = get_target(lim_frac(ldx + ldz, frac_bits));
                        test_fxp("add   (n1+n3)", fxp1, tgt2);
                        test_fxp("add_l (n1+n3)", fxp2, tgt2);

                        fxp1 = fxp_mul(n1, n2);
                        fxp2 = fxp_mul_l(n1, n2);
                        tgt1 = get_target(lim_frac(ldx * ldy, frac_bits));
                        test_fxp("mul   (n1*n2)", fxp1, tgt1);
                        test_fxp("mul_l (n1*n2)", fxp2, tgt1);
                        fxp1 = fxp_mul(n1, n3);
                        fxp2 = fxp_mul_l(n1, n3);
                        tgt2 = get_target(lim_frac(ldx * ldz, frac_bits));
                        test_fxp("mul   (n1*n3)", fxp1, tgt2);
                        test_fxp("mul_l (n1*n3)", fxp2, tgt2);
                        // Skipping testing of fxp_mul_d since it uses longs
                        // (like fxp_mul_l), but it is less efficient. No
                        // reason to chose it over fxp_mul_l

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
                }

                printf("\n%d Warnings for %d frac bits.\n", \
                            nwarnings, frac_bits);
                if (nwarnings > 0) {
                    printf("Largest delta was: %LE (max allowed: %LE)\n", \
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
            printf("Largest warning delta: %LE (max allowed:%LE, using %d frac bits)\n", \
                        largest_delta, largest_madelta, largest_delta_fbits);
        }
        return 0;
}
