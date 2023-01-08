/*
 * fxp_tester.c
 * Tests the implementation of binary fixed point (fxp.c)
 * By Raul Saavedra, 2022-11-13
 */

#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include "fxp.h"

#define DASHES "========================\n"

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

void test_fxp(char *s, int fxp, int assert_val)
{
        printf("'%s', got ", s);
        print_fxp(fxp);
        printf(",\texpected ");
        if (fxp == assert_val)
            printf("same!");
        else
            print_fxp(assert_val);
        printf("\n");
        fflush( stdout );
        // Allowing 1 bit error (least significant bit in frac part)
        int delta = 0;
        if (fxp >= assert_val) {
                delta = fxp - assert_val;
                assert( (delta <= 1) &&
                    ((fxp >= 0 && assert_val >=0) ||
                        (fxp < 0 && assert_val < 0)));
        } else {
                delta = assert_val - fxp;
                assert( (delta <= 1) &&
                    ((fxp >= 0 && assert_val >=0) ||
                        (fxp < 0 && assert_val < 0)));
        }
        if (delta == 1) {
            nwarnings++;
            printf("***** Warning %d: delta == 1 found!\n", nwarnings);
        }
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

        printf("\nTesting function that counts in O(log(n)) the # bits used by a number:");
        unsigned int n = 1;
        int nb = 1;
        do {
                printf("\nUnsigned number %u uses %d bits", n-1, \
                        fxp_nbits(n-1, FXP_INT_BITS));
                printf("\nUnsigned Number %u uses %d bits", n, \
                        fxp_nbits(n, FXP_INT_BITS));
                n = (n << 1);
                nb++;
        } while (nb <= FXP_INT_BITS);
        printf("\nUnsigned Number %u uses %d bits", UINT_MAX, \
                fxp_nbits(UINT_MAX, FXP_INT_BITS));
        printf("\n");

        int afb[] = {0, 1, 10, 12, 15, 24, 31};
        for (int fb = 0; fb < 7; fb++) {

            fxp_set_frac_bits(afb[fb]);
            fxp_set_auto_frac_max_dec();

            int frac_bits    = fxp_get_frac_bits();
            int whole_bits   = fxp_get_whole_bits();
            int frac_mask    = fxp_get_frac_mask();
            int frac_max     = fxp_get_frac_max();
            int frac_max_dec = fxp_get_frac_max_dec();
            int whole_max    = fxp_get_whole_max();
            int whole_min    = fxp_get_whole_min();

            printf("\n%sFXP configuration parameters:\n", DASHES);
            printf("pos infinity: %d\n", FXP_POS_INF);
            printf("whole bits  : %d\n", whole_bits);
            printf("frac bits   : %d\n", frac_bits);
            printf("whole max   : %d\n", whole_max);
            printf("whole min   : %d\n", whole_min);
            printf("neg infinity: %d\n", FXP_NEG_INF);
            printf("undefined   : %d\n", FXP_UNDEF);
            printf("frac mask   : %d\n", frac_mask);
            printf("frac max    : %d (->decimals: .%d)\n",
                                        frac_max,
                                        frac_max_dec);

            printf("\nSome general tests:\n");
            int fxp_zero     = fxp_bin(0, 0);
            int fxp_tiniest  = frac_bits > 0? fxp_bin(0, 1): fxp_bin(1, 0);
            int fxp_mtiniest = -fxp_tiniest;
            int fxp_largest  = FXP_MAX;

            test_fxp("Infinity",  FXP_POS_INF, FXP_POS_INF);
            test_fxp("-Infinity", FXP_NEG_INF, FXP_NEG_INF);
            test_fxp("Undefined", FXP_UNDEF,   FXP_UNDEF);
            test_fxp("zero",     fxp_zero,      fxp_bin(0, 0));
            test_fxp("tiniest",  fxp_tiniest,   fxp_bin(0, 1));
            test_fxp("-tiniest", fxp_mtiniest,  -fxp_tiniest);
            test_fxp("Largest",  fxp_largest,
                                 fxp_bin(whole_max,
                                 frac_max - (whole_bits > 1? 1: 0)));
            test_fxp("Largest - tiniest", fxp_sub(fxp_largest, fxp_tiniest),
                                          fxp_largest - 1);
            test_fxp("Largest + tiniest safe",
                                fxp_sum(fxp_largest, fxp_tiniest),
                                FXP_POS_INF);
            test_fxp("Largest + tiniest unsafe",
                                fxp_unsafe_sum(fxp_largest, fxp_tiniest),
                                FXP_POS_INF);
            test_fxp("Most negative",
                                -fxp_largest,
                                FXP_MIN);
            test_fxp("Almost most negative",
                                fxp_sum(-fxp_largest, fxp_tiniest),
                                -fxp_largest + fxp_tiniest);
            test_fxp("+inf + +inf",
                                fxp_sum(FXP_POS_INF, FXP_POS_INF),
                                FXP_POS_INF);
            test_fxp("-inf - +inf",
                                fxp_sub(FXP_NEG_INF, FXP_POS_INF),
                                FXP_NEG_INF);
            test_fxp("+inf + -inf",
                                fxp_sum(FXP_POS_INF, FXP_NEG_INF),
                                FXP_UNDEF);
            test_fxp("-inf + -inf",
                                fxp_sum(FXP_NEG_INF, FXP_NEG_INF),
                                FXP_NEG_INF);
            test_fxp("-inf - -inf",
                                fxp_sub(FXP_NEG_INF, FXP_NEG_INF),
                                FXP_UNDEF);
            test_fxp("+inf * -inf",
                                fxp_mul(FXP_POS_INF, FXP_NEG_INF),
                                FXP_NEG_INF);
            test_fxp("+num / zero",
                                fxp_div(fxp_largest, fxp_zero),
                                FXP_POS_INF);
            test_fxp("zero / zero",
                                fxp_div(fxp_zero, fxp_zero),
                                FXP_UNDEF);
            test_fxp("zero * zero",
                                fxp_mul(fxp_zero, fxp_zero),
                                fxp_zero);
            test_fxp("zero + zero",
                                fxp_sum(fxp_zero, fxp_zero),
                                fxp_zero);
            test_fxp("zero - zero",
                                fxp_sub(fxp_zero, fxp_zero),
                                fxp_zero);
            test_fxp("-num / zero",
                                fxp_div(-fxp_largest, fxp_zero),
                                FXP_NEG_INF);
            test_fxp("zero * +inf",
                                fxp_mul(fxp_zero, FXP_POS_INF),
                                FXP_UNDEF);
            test_fxp("zero * -inf",
                                fxp_mul(fxp_zero, FXP_NEG_INF),
                                FXP_UNDEF);
            test_fxp("zero * undef",
                                fxp_mul(fxp_zero, FXP_UNDEF),
                                FXP_UNDEF);
            test_fxp("-inf * undef",
                                fxp_mul(FXP_NEG_INF, FXP_UNDEF),
                                FXP_UNDEF);
            test_fxp("+inf * undef",
                                fxp_mul(FXP_POS_INF, FXP_UNDEF),
                                FXP_UNDEF);
            test_fxp("undef * undef",
                                fxp_mul(FXP_UNDEF, FXP_UNDEF),
                                FXP_UNDEF);
            test_fxp("tiniest * inf",
                                fxp_mul(fxp_tiniest, FXP_POS_INF),
                                FXP_POS_INF);

            if (frac_bits > 0) {
                    printf("\nTests that apply when using at least one frac bit\n");
                    printf("\nChecking decimal <=> bin mappings of frac ranges:");
                    int fmbin = fxp_get_frac_max();
                    int fmdec = fxp_get_frac_max_dec();
                    printf("\nMax frac dec: %d (bin %d)", fmdec, fmbin);
                    for (int i = 0; i <= 5; i++) {
                            printf("\nShowing fxp for 0.%3d: ", i);
                            int vf = fxp_dec(0, i);
                            print_fxp(vf);
                    }
                    printf("\n");
                    int m = (fxp_get_frac_max_dec() + 1) / 2;
                    for (int i = (m >= 2 ? m - 2: 0); i <= m + 2; i++) {
                            printf("\nShowing fxp for 0.%3d: ", i);
                            int vf = fxp_dec(0, i);
                            print_fxp(vf);
                    }
                    printf("\n");
                    for (int i = (fmdec >= 5? fmdec - 5: 0);
                                i <= fmdec; i++) {
                            printf("\nShowing fxp for 0.%3d: ", i);
                            int vf = fxp_dec(0, i);
                            print_fxp(vf);
                    }
                    printf("\n");

                    printf("\nChecking sign taken from frac when whole==0\n");
                    fxp_set_frac_max_dec(999);
                    test_fxp("-0.(+)500",
                                        fxp_dec(-0, 500),
                                        fxp_dec(0, 500));
                    test_fxp("-0.(-)500",
                                        fxp_dec(-0, -500),
                                        fxp_dec(0, -500));

                    printf("\nTruncation of longer frac decimal arguments\n");
                    test_fxp("0.222222", fxp_dec(0, 222222), fxp_dec(0, 222));
                    test_fxp("0.777777", fxp_dec(0, 777777), fxp_dec(0, 777));
                    test_fxp("0.991999", fxp_dec(0, 991999), fxp_dec(0, 991));
                    test_fxp("0.999999", fxp_dec(0, 999999), fxp_dec(0, 999));
                    fxp_set_frac_max_dec(fmdec);
            }

            if (whole_bits > 1) {
                    printf("\nTests that apply when using at least two bits for whole part (1-bit sign, 1-bit value):\n");
                    test_fxp("Way Too Large whole part!",
                                        fxp_sum(fxp(whole_max), fxp(5)),
                                        FXP_POS_INF);
                    int halfmax = fxp_div(FXP_MAX, fxp(2));
                    test_fxp("Largest * 1", fxp_mul(fxp_largest,  fxp(1)), fxp_largest);
                    test_fxp("Largest * -1", fxp_mul(fxp_largest,  fxp(-1)), -fxp_largest);
                    int fxp_one = fxp(1);
                    test_fxp("one", fxp_one, fxp_dec(1, 0));
                    test_fxp("Largest + one safe",
                                        fxp_sum(fxp_largest, fxp_one),
                                        FXP_POS_INF);
                    test_fxp("Largest + one unsafe",
                                        fxp_unsafe_sum(fxp_largest, fxp_one),
                                        fxp_largest + fxp_one);
                    test_fxp("Safe Too neg substraction",
                                        fxp_sub( -fxp_largest, fxp_one),
                                        FXP_NEG_INF);
                    test_fxp("Unsafe Too neg substraction",
                                        fxp_unsafe_sub( -fxp_largest, fxp_one),
                                        -fxp_largest - fxp_one);
                    test_fxp("+inf - 1",
                                        fxp_sub(FXP_POS_INF, fxp_one),
                                        FXP_POS_INF);
                    test_fxp("-inf + 1",
                                        fxp_sum(FXP_NEG_INF, fxp_one),
                                        FXP_NEG_INF);
                    test_fxp("Half Max + Half Max",
                                        fxp_sum(halfmax, halfmax),
                                        halfmax+halfmax);
                    test_fxp("FXP_MAX - Half Max",
                                        fxp_sub(FXP_MAX, halfmax),
                                        FXP_MAX/2);
                    test_fxp("Half Max + FXP_MAX",
                                        fxp_sum(halfmax, FXP_MAX),
                                        FXP_POS_INF);
                    test_fxp("-FXP_MAX - Half Max",
                                        fxp_sub(-FXP_MAX, halfmax),
                                        FXP_NEG_INF);
                    test_fxp("Half Max * 2",
                                        fxp_mul(halfmax, fxp(2)),
                                        FXP_MAX);
                    test_fxp("Half Max * 2 (long)",
                                        fxp_mul_using_long(halfmax, fxp(2)),
                                        FXP_MAX);
                    test_fxp("Half Max * 3",
                                        fxp_mul(halfmax, fxp(3)),
                                        FXP_POS_INF);
                    test_fxp("-Half Max * 3",
                                        fxp_mul(-halfmax, fxp(3)),
                                        FXP_NEG_INF);
                    test_fxp("Half Max / 0.25",
                                        fxp_div(halfmax, fxp_dec(0, 250)),
                                        FXP_POS_INF);
                    test_fxp("-Half Max / 0.25",
                                        fxp_div(-halfmax, fxp_dec(0, 250)),
                                        FXP_NEG_INF);

                    test_fxp("(HalfMax+1)*2",
                                        fxp_mul(halfmax + fxp(1), fxp(2)),
                                        FXP_POS_INF);
                    test_fxp("(HalfMax+1)*2 (long)",
                                        fxp_mul_using_long(halfmax + fxp(1), fxp(2)),
                                        FXP_POS_INF);
            }

            if ((frac_bits > 0) && (whole_bits > 1)) {
                    printf("\nTests when using both whole and fraction parts:\n");
                    printf("\nSimple operations checking signs\n");
                    int whole    = whole_bits >= 4? 10: 0;
                    int bin_frac = frac_max / 4;  // == 0.250
                    int num      = fxp_bin(whole, bin_frac);
                    int dec_frac = fxp_get_dec_frac(num);
                    //printf("Whole %d, frac %d, num is ", whole, bin_frac);
                    print_fxp(num); printf("\n");
                    int fxp1    = num;
                    int fxp2    = fxp(2);
                    int nump2   = fxp_bin( whole + 2, bin_frac);    // num + 2
                    int numm2   = fxp_bin( whole - 2, bin_frac);    // num - 2
                    int mnumm2  = fxp_bin(-whole - 2, bin_frac);    // -num - 2
                    int mnump2  = fxp_bin(-whole + 2, bin_frac);    // -num + 2
                    int dnum    = fxp_bin( whole * 2, 2*bin_frac);  // num*2
                    int mdnum   = fxp_bin(-whole * 2, 2*bin_frac);  // -num*2
                    int hnum    = fxp_bin( whole / 2, bin_frac/2);  // half num
                    int mhnum   = fxp_bin(-whole / 2, bin_frac/2);  // - half num
                    test_fxp(" num +  2", fxp_sum( fxp1,  fxp2), nump2);
                    test_fxp(" num + -2", fxp_sum( fxp1, -fxp2), numm2);
                    test_fxp("-num +  2", fxp_sum(-fxp1,  fxp2), mnump2);
                    test_fxp("-num + -2", fxp_sum(-fxp1, -fxp2), mnumm2);
                    test_fxp(" num -  2", fxp_sub( fxp1,  fxp2), numm2);
                    test_fxp(" num - -2", fxp_sub( fxp1, -fxp2), nump2);
                    test_fxp("-num -  2", fxp_sub(-fxp1,  fxp2), mnumm2);
                    test_fxp("-num - -2", fxp_sub(-fxp1, -fxp2), mnump2);
                    test_fxp(" num *  2", fxp_mul( fxp1,  fxp2), dnum);
                    test_fxp(" num * -2", fxp_mul( fxp1, -fxp2), mdnum);
                    test_fxp("-num *  2", fxp_mul(-fxp1,  fxp2), mdnum);
                    test_fxp("-num * -2", fxp_mul(-fxp1, -fxp2), dnum);
                    test_fxp(" num *  2 (long)",
                            fxp_mul_using_long( fxp1,  fxp2), dnum);
                    test_fxp(" num * -2 (long)",
                            fxp_mul_using_long( fxp1, -fxp2), mdnum);
                    test_fxp("-num *  2 (long)",
                            fxp_mul_using_long(-fxp1,  fxp2), mdnum);
                    test_fxp("-num * -2 (long)",
                            fxp_mul_using_long(-fxp1, -fxp2), dnum);
                    test_fxp(" num /  2", fxp_div( fxp1,  fxp2), hnum);
                    test_fxp(" num / -2", fxp_div( fxp1, -fxp2), mhnum);
                    test_fxp("-num /  2", fxp_div(-fxp1,  fxp2), mhnum);
                    test_fxp("-num / -2", fxp_div(-fxp1, -fxp2), hnum);

                    if (frac_bits == 12) {
                            int m1 = fxp_bin(1, 2048); // == 1.5
                            printf("m1 is %u\n", m1);
                            test_fxp("1.5 * 1.5", fxp_mul(m1, m1),
                                    fxp_bin(2, 1024));
                            test_fxp("723.5 * 723.5",
                                    fxp_mul(fxp_bin(723, 2048), fxp_bin(723, 2048)),
                                    fxp_bin(523452, 1024));
                            test_fxp("-720.5 * -730.5",
                                    fxp_mul(fxp_dec(-720, 500), fxp_dec(-730, 500)),
                                    FXP_POS_INF);
                            test_fxp("0.999 * 0.005",
                                    fxp_mul(fxp_bin(0,4095), fxp_bin(0, 24)),
                                    fxp_bin(0, 23));
                    }
            }

            printf("\nTotal # of warnings: %d\n", nwarnings);
            printf("All tests passed using %d-bit fracs, ", frac_bits);
            printf("and '%d' as max decimal frac.\n\n", frac_max_dec);
        }

        return 0;
}
