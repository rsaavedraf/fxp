/*
 * fxp_tester.c
 * Tests the implementation of binary fixed point (fxp.c)
 * By Raul Saavedra, 2022-11-13, Bonn Germany
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <float.h>
#include <time.h>
#include "fxp.h"

#define DASHES "========================\n"

static int fracbit_configs[] = {0, 11, 14, 16, 24, 31};
static long double max_ok_delta = 0.0;
static long double larger_delta = 0.0;
static long double largest_delta = 0.0;
static int largest_delta_fbits = 0;
static long double max_delta = 0.0;
static int twarnings = 0;
static int nwarnings = 0;

void print_fxp(int fxp)
{
        if (fxp == FXP_POS_INF || fxp == FXP_NEG_INF || fxp == FXP_UNDEF)
                printf(fxp==FXP_UNDEF? "UNDEF":
                      (fxp==FXP_POS_INF? "+INF": "-INF"));
        else {
                int whole = fxp_get_whole_part(fxp);
                char * sign  = ((fxp < 0) && (whole == 0))? "-": "";
                //char * sign  = ""; // ((fxp<0) && (whole==0))? "-": "";
                printf("%s%d.%3d (bin frac=%d)",
                        sign,
                        whole,
                        fxp_get_dec_frac(fxp),
                        fxp_get_bin_frac(fxp));
        }
}

/*
 * Return a double value corresponding to a given fxp
 */
long double dfxp(int fxp)
{
    int wi = fxp_get_whole_part(fxp);
    int fi = fxp_get_bin_frac(fxp);
    //printf("\nfxpi w - f: %d - %d\n", wi, fi);
    long double wd = (long double) wi;
    long double fd = ((long double) fi) / fxp_get_frac_max();
    //printf("fxpd w - f: %Lf - %Lf\n", wd, fd);
    long double x = wd + fd ;
    return x;
}

void test_fxp(char *s, int fxp, long double d_assert_val)
{
        printf("%s, got ", s);
        print_fxp(fxp);
        printf(",\texpected %Lf", d_assert_val);
        long double df = dfxp(fxp);
        long double delta = (df >= d_assert_val)?
                        df - d_assert_val:
                        d_assert_val - df;

        if (delta <= max_ok_delta) {
            // No warning up to the max_ok_delta value
            printf(" (~same)\n");
        } else {
            nwarnings++;
            printf("\n***** Warning %d: delta=%LE\n\n",
                        nwarnings, delta);
        }
        if (delta > larger_delta)
            larger_delta = delta;

        fflush( stdout );
        assert( (delta <= max_delta) &&
            ((df >= 0 && d_assert_val >= 0) ||
                (df < 0 && d_assert_val < 0)));

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
        srand((unsigned int) time(0));

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

                int aux = (int) (frac_bits / 2.25);
                max_delta = ((long double) (frac_mask >> aux)) / frac_max;
                max_ok_delta = max_delta / 1.75;

                nwarnings = 0;
                larger_delta = 0;
                printf("\n%sFXP configuration parameters:\n", DASHES);
                printf("frac bits   : %d (requested was %d)\n",
                                frac_bits, fracbit_configs[nfb]);
                printf("whole bits  : %d\n", whole_bits);
                printf("pos infinity: %d (Lf: %Lf)\n",
                                FXP_POS_INF, dfxp(FXP_POS_INF));
                printf("whole max   : %d\n", whole_max);
                printf("frac mask   : %d\n", frac_mask);
                printf("frac max    : %d (->decimals: .%d)\n",
                                frac_max,
                                frac_max_dec);
                printf("whole min   : %d\n", whole_min);
                printf("neg infinity: %d (Lf: %Lf)\n",
                                FXP_NEG_INF, dfxp(FXP_NEG_INF));
                printf("undefined   : %d (Lf: %Lf)\n",
                                FXP_UNDEF, dfxp(FXP_UNDEF));
                printf("max_ok_delta: %LE\n", max_ok_delta);
                printf("max_delta   : %LE\n", max_delta);

                printf("\nChecking extreme int values, part I:\n");
                int fxp_largest = FXP_MAX;
                int fxp_tiniest = 1;
                int fxp_ten = fxp_bin(10, 0);
                int fxp_two = fxp_bin(2, 0);
                int fxp_one = fxp_bin(1, 0);
                int fxp_p5 = fxp_bin(0, frac_max/2 + 1);
                int fxp_zero = fxp_bin(0, 0);
                int fxp_halfmax = FXP_MAX / 2;
                int fxp_halfp2 = fxp_add(fxp_halfmax, fxp_two);
                test_fxp("Infinity", FXP_POS_INF, dfxp(FXP_POS_INF));
                test_fxp("Largest", fxp_largest, dfxp(FXP_MAX));
                test_fxp("HalfMax", fxp_halfmax, dfxp(FXP_MAX)/2);
                test_fxp("Largest frac", frac_max, dfxp(fxp_bin(0, frac_max)));
                test_fxp("tiniest", fxp_tiniest,   dfxp(fxp_tiniest));
                test_fxp("0.5", fxp_p5, ((long double) 0.5));
                test_fxp("zero", fxp_zero, dfxp(0));
                test_fxp("-tiniest", -fxp_tiniest, dfxp(-fxp_tiniest));
                test_fxp("-Largest frac", -frac_max, dfxp(-frac_max));
                test_fxp("Most negative", -fxp_largest, dfxp(FXP_MIN));
                test_fxp("-Infinity", FXP_NEG_INF, dfxp(FXP_NEG_INF));
                test_fxp("Undefined", FXP_UNDEF,   dfxp(FXP_UNDEF));

                printf("\nChecking extreme int values, part II:\n");
                test_fxp("Almost most negative",
                                fxp_add(-fxp_largest, fxp_tiniest),
                                dfxp(-fxp_largest) + dfxp(fxp_tiniest));
                test_fxp(" Largest -Largest",
                                fxp_add(fxp_largest, -fxp_largest),
                                dfxp(fxp_zero));
                test_fxp("-Largest +Largest",
                                fxp_add(-fxp_largest, fxp_largest),
                                dfxp(fxp_zero));
                test_fxp("Largest + 0",
                                fxp_add(fxp_largest, fxp_zero),
                                dfxp(fxp_largest));
                test_fxp("-Largest - 0",
                                fxp_add(-fxp_largest, -fxp_zero),
                                dfxp(-fxp_largest));
                test_fxp("Largest - tiniest",
                                fxp_sub(fxp_largest, fxp_tiniest),
                                dfxp(fxp_largest) - dfxp(fxp_tiniest));
                test_fxp("Largest + tiniest safe",
                                fxp_add(fxp_largest, fxp_tiniest),
                                dfxp(FXP_POS_INF));
                test_fxp("Largest + tiniest unsafe",
                                fxp_unsafe_add(fxp_largest, fxp_tiniest),
                                dfxp(FXP_POS_INF));
                test_fxp("+inf + +inf",
                                fxp_add(FXP_POS_INF, FXP_POS_INF),
                                dfxp(FXP_POS_INF));
                test_fxp("-inf - +inf",
                                fxp_sub(FXP_NEG_INF, FXP_POS_INF),
                                dfxp(FXP_NEG_INF));
                test_fxp("+inf + -inf",
                                fxp_add(FXP_POS_INF, FXP_NEG_INF),
                                dfxp(FXP_UNDEF));
                test_fxp("-inf + -inf",
                                fxp_add(FXP_NEG_INF, FXP_NEG_INF),
                                dfxp(FXP_NEG_INF));
                test_fxp("-inf - -inf",
                                fxp_sub(FXP_NEG_INF, FXP_NEG_INF),
                                dfxp(FXP_UNDEF));
                test_fxp("+inf * -inf",
                                fxp_mul(FXP_POS_INF, FXP_NEG_INF),
                                dfxp(FXP_NEG_INF));
                test_fxp("+inf - 0.5",
                                fxp_sub(FXP_POS_INF, fxp_p5),
                                dfxp(FXP_POS_INF));
                test_fxp("-inf + 0.5",
                                fxp_add(FXP_NEG_INF, fxp_p5),
                                dfxp(FXP_NEG_INF));
                test_fxp("+num / zero",
                                fxp_div(fxp_largest, fxp_zero),
                                dfxp(FXP_POS_INF));
                test_fxp("zero / zero",
                                fxp_div(fxp_zero, fxp_zero),
                                dfxp(FXP_UNDEF));
                test_fxp("zero * zero",
                                fxp_mul(fxp_zero, fxp_zero),
                                dfxp(fxp_zero));
                test_fxp("zero + zero",
                                fxp_add(fxp_zero, fxp_zero),
                                dfxp(fxp_zero));
                test_fxp("zero - zero",
                                fxp_sub(fxp_zero, fxp_zero),
                                dfxp(fxp_zero));
                test_fxp("-num / zero",
                                fxp_div(-fxp_largest, fxp_zero),
                                dfxp(FXP_NEG_INF));
                test_fxp("zero * +inf",
                                fxp_mul(fxp_zero, FXP_POS_INF),
                                dfxp(FXP_UNDEF));
                test_fxp("zero * -inf",
                                fxp_mul(fxp_zero, FXP_NEG_INF),
                                dfxp(FXP_UNDEF));
                test_fxp("zero * undef",
                                fxp_mul(fxp_zero, FXP_UNDEF),
                                dfxp(FXP_UNDEF));
                test_fxp("-inf * undef",
                                fxp_mul(FXP_NEG_INF, FXP_UNDEF),
                                dfxp(FXP_UNDEF));
                test_fxp("+inf * undef",
                                fxp_mul(FXP_POS_INF, FXP_UNDEF),
                                dfxp(FXP_UNDEF));
                test_fxp("undef * undef",
                                fxp_mul(FXP_UNDEF, FXP_UNDEF),
                                dfxp(FXP_UNDEF));
                test_fxp("tiniest * inf",
                                fxp_mul(fxp_tiniest, FXP_POS_INF),
                                dfxp(FXP_POS_INF));

                printf("\nChecking extreme int values, part III:\n");
                test_fxp("Way Too Large whole part!",
                                fxp_add(fxp(whole_max), fxp(5)),
                                dfxp(FXP_POS_INF));
                test_fxp("Largest * 1",
                                fxp_mul(fxp_largest,  fxp_one),
                                dfxp(fxp_largest));
                test_fxp("Largest * -1",
                                fxp_mul(fxp_largest,  -fxp_one),
                                dfxp(-fxp_largest));
                test_fxp("Largest + two safe",
                                fxp_add(fxp_largest, fxp_two),
                                dfxp(FXP_POS_INF));
                test_fxp("Largest + two unsafe",
                                fxp_unsafe_add(fxp_largest, fxp_two),
                                dfxp(fxp_largest + fxp_two));
                test_fxp("Safe Too neg substraction",
                                fxp_sub( -fxp_largest, fxp_one),
                                dfxp(FXP_NEG_INF));
                test_fxp("Unsafe Too neg substraction",
                                fxp_unsafe_sub( -fxp_largest, fxp_p5),
                                dfxp(-fxp_largest - fxp_p5));
                test_fxp("Largest + 0.5",
                                fxp_add(fxp_largest, fxp_p5),
                                dfxp(FXP_POS_INF));
                test_fxp("-Largest - 0.5",
                                fxp_add(-fxp_largest, -fxp_p5),
                                dfxp(FXP_NEG_INF));
                test_fxp("+HalfMax + HMaxp2",
                                fxp_add(fxp_halfmax, fxp_halfp2),
                                dfxp(FXP_POS_INF));
                test_fxp("-HalfMax - HMaxp2",
                                fxp_add(-fxp_halfmax, -fxp_halfp2),
                                dfxp(FXP_NEG_INF));
                test_fxp("HalfMax + HalfMax",
                                fxp_add(fxp_halfmax, fxp_halfmax),
                                dfxp(fxp_halfmax + fxp_halfmax));
                test_fxp("FXP_MAX - HalfMax",
                                fxp_sub(FXP_MAX, fxp_halfmax),
                                dfxp(FXP_MAX - fxp_halfmax));
                test_fxp("HalfMax + FXP_MAX",
                                fxp_add(fxp_halfmax, FXP_MAX),
                                dfxp(FXP_POS_INF));
                test_fxp("-FXP_MAX - HalfMax",
                                fxp_sub(-FXP_MAX, fxp_halfmax),
                                dfxp(FXP_NEG_INF));
                test_fxp("HalfMax * 2",
                                fxp_mul(fxp_halfmax, fxp(2)),
                                dfxp(FXP_MAX));
                test_fxp("HalfMax * 2 (long)",
                                fxp_mul_l(fxp_halfmax, fxp(2)),
                                dfxp(FXP_MAX));
                test_fxp("HalfMax * 3",
                                fxp_mul(fxp_halfmax, fxp(3)),
                                dfxp(FXP_POS_INF));
                test_fxp("-HalfMax * 3",
                                fxp_mul(-fxp_halfmax, fxp(3)),
                                dfxp(FXP_NEG_INF));
                test_fxp("(HalfMax+0.5)*2",
                                fxp_mul(fxp_add(fxp_halfmax, fxp_p5), fxp(2)),
                                dfxp(FXP_POS_INF));
                test_fxp("(HalfMax+0.5)*2 (long)",
                                fxp_mul_l(fxp_add(fxp_halfmax, fxp_p5), fxp(2)),
                                dfxp(FXP_POS_INF));

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
                // also relaxing the delta vars for them
                long double bkp_ld = larger_delta;
                long double bkp_mokd = max_ok_delta;
                long double bkp_md = max_delta;
                fxp_set_frac_max_dec(999);
                max_ok_delta = (max_ok_delta > 0.0015? max_ok_delta: 0.0015);
                max_delta = (max_delta > 0.002? max_delta: 0.002);

                test_fxp("-0.(+)500", fxp_dec(-0, 500), 0.5);
                test_fxp("-0.(-)500", fxp_dec(-0, -500), -0.5);

                printf("\nTruncation of longer frac decimal arguments:\n");
                test_fxp("0.222222", fxp_dec(0, 222222), 0.222);
                test_fxp("0.777777", fxp_dec(0, 777777), 0.777);
                test_fxp("0.991999", fxp_dec(0, 991999), 0.991);
                test_fxp("0.999999", fxp_dec(0, 999999), 0.999);
                // Restoring original frac_max_dec and delta vars
                fxp_set_frac_max_dec(fmdec);
                max_ok_delta = bkp_mokd;
                max_delta = bkp_md;
                larger_delta = bkp_ld;

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
                        fxp_add(frac_max, -frac_max), 0.0);

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

                printf("\nVerifying different op ");
                printf("implementations using random numbers:\n");
                int sign1, sign2, sign3, n1, n2, n3, fxp1, fxp2, fxp3, fxp4;
                for (int i=0; i < 10; i++) {
                        sign1 = rand() % 2 == 1? -1: 1;
                        sign2 = rand() % 2 == 1? -1: 1;
                        sign3 = rand() % 2 == 1? -1: 1;
                        n1 = sign1 * rand();
                        n2 = sign2 * rand();
                        n3 = sign3 * rand();
                        if (frac_bits < FXP_INT_BITS_M1) n3 %= fxp_get_frac_max();
                        printf("n1 = "); print_fxp(n1); printf("\n");
                        printf("n2 = "); print_fxp(n2); printf("\n");
                        printf("n3 = "); print_fxp(n3);
                        printf(" (n3 always random between -1 and 1)\n");

                        fxp1 = fxp_add(n1, n2);
                        fxp2 = fxp_add_l(n1, n2);
                        fxp3 = fxp_add(n1, n3);
                        fxp4 = fxp_add_l(n1, n3);
                        test_fxp("a. sum vs. sum_l", fxp1, dfxp(fxp2));
                        test_fxp("b. sum vs. sum_l", fxp3, dfxp(fxp4));

                        fxp1 = fxp_mul(n1, n2);
                        fxp2 = fxp_mul_l(n1, n2);
                        fxp3 = fxp_mul_d(n1, n2);
                        test_fxp("c. mul vs. mul_l", fxp1, dfxp(fxp2));
                        test_fxp("d. mul vs. mul_d", fxp1, dfxp(fxp3));
                        test_fxp("e. mul_l vs. mul_d", fxp2, dfxp(fxp3));
                        fxp1 = fxp_mul(n1, n3);
                        fxp2 = fxp_mul_l(n1, n3);
                        fxp3 = fxp_mul_d(n1, n3);
                        test_fxp("f. mul vs. mul_l", fxp1, dfxp(fxp2));
                        test_fxp("g. mul vs. mul_d", fxp1, dfxp(fxp3));
                        test_fxp("h. mul_l vs. mul_d", fxp2, dfxp(fxp3));

                        fxp1 = fxp_div(n1, n2);
                        fxp2 = fxp_div_l(n1, n2);
                        test_fxp("i. div vs. div_l", fxp1, dfxp(fxp2));
                }

                printf("\n%d Warnings for %d frac bits.\n", \
                            nwarnings, frac_bits);
                printf("Largest delta was: %LE\n", larger_delta);
                printf("All tests passed using %d-bit fracs, ", frac_bits);
                printf("and '%d' as max decimal frac.\n\n", frac_max_dec);
                twarnings += nwarnings;

                if (larger_delta > largest_delta) {
                        largest_delta = larger_delta;
                        largest_delta_fbits = frac_bits;
                }
        }
        printf("Grand total of %d warnings checking %d configurations.\n", \
                    twarnings, nconfigs);
        printf("Largest delta of all: %LE (using %d frac bits)\n", \
                    largest_delta, largest_delta_fbits);

        return 0;
}
