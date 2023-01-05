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
                char * sign  = ((fxp<0) && (whole==0))? "-": "";
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
        print_fxp(assert_val);
        printf("\n");
        fflush( stdout );
        // Allowing 1 bit error (least significant bit in frac part)
        int delta = 0;
        if (fxp >= assert_val) {
                delta = fxp - assert_val;
                assert( delta <= 1);
        } else {
                delta = assert_val - fxp;
                assert( delta <= 1);
        }
        if (delta == 1) {
            nwarnings++;
            printf("***** Warning %d: delta == 1 found!\n", nwarnings);
        }
}

int main(void)
{
        printf("%sFXP Tester run\n%s", DASHES, DASHES);

        printf("Num type sizes in this system:\n");
        printf("char        has a size of %zd bytes.\n", sizeof(char));
        printf("int         has a size of %zd bytes.\n", sizeof(int));
        printf("long        has a size of %zd bytes.\n", sizeof(long));
        printf("long long   has a size of %zd bytes.\n", sizeof(long long));
        printf("float       has a size of %zd bytes.\n", sizeof(float));
        printf("double      has a size of %zd bytes.\n", sizeof(double));
        printf("long double has a size of %zd bytes.\n", sizeof(long double));

        printf("\nFXP configuration constants:\n");
        printf("frac bits   : %d\n", FXP_FRAC_BITS);
        printf("max fraction: %d (->decimals: .%d)\n", FXP_FRAC_MAX, FXP_FRAC_MAX_DEC );
        printf("max whole   : %d\n", FXP_WHOLE_MAX);
        printf("min whole   : %d\n", FXP_WHOLE_MIN);
        printf("pos infinity: %d\n", FXP_POS_INF);
        printf("neg infinity: %d\n", FXP_NEG_INF);
        printf("undefined   : %d\n", FXP_UNDEF);

        printf("\nTesting function that counts number of bits used by numbers:");
        int n = 1;
        int nb = 1;
        do {
            printf("\n Number %u uses %u bits", n-1, nbits(n-1, FXP_INT_BITS));
            printf("\n Number %u uses %u bits", n, nbits(n, FXP_INT_BITS));
            n = (n << 1);
            nb++;
        } while (nb <= FXP_INT_BITS);
        printf("\n");

        printf("\nChecking decimal<=>bin mappings of frac ranges:");
        for (int i=0; i<=10; i++) {
                printf("\nShowing fxp for 10.%3d: ", i);
                print_fxp(fxp_dec(10, i));
        }
        printf("\n");
        int m = (FXP_FRAC_MAX_DEC+1)/2;
        for (int i=m-5; i<=m+5; i++) {
                printf("\nShowing fxp for 10.%3d: ", i);
                print_fxp(fxp_dec(10, i));
        }
        printf("\n");
        for (int i=FXP_FRAC_MAX_DEC-10; i<=FXP_FRAC_MAX_DEC; i++) {
                printf("\nShowing fxp for 10.%3d: ", i);
                print_fxp(fxp_dec(10, i));
        }
        printf("\n");

        printf("\nSimple operations with no overflows, checking signs\n");
        int fxp1 = fxp_dec(10, 250);
        int fxp2 = fxp(2);
        int a12  = fxp_dec(12, 250);
        int am12 = fxp_dec(-12, -250);
        int a8   = fxp_dec(8, 250);
        int am8  = fxp_dec(-8, -250);
        int a20  = fxp_dec(20, 500);
        int am20 = fxp_dec(-20, -500);
        int a5   = fxp_dec(5, 125);
        int am5  = fxp_dec(-5, -125);
        test_fxp(" 10.25 +  2", fxp_sum( fxp1,  fxp2), a12);
        test_fxp(" 10.25 + -2", fxp_sum( fxp1, -fxp2), a8);
        test_fxp("-10.25 +  2", fxp_sum(-fxp1,  fxp2), am8);
        test_fxp("-10.25 + -2", fxp_sum(-fxp1, -fxp2), am12);
        test_fxp(" 10.25 -  2", fxp_sub( fxp1,  fxp2), a8);
        test_fxp(" 10.25 - -2", fxp_sub( fxp1, -fxp2), a12);
        test_fxp("-10.25 -  2", fxp_sub(-fxp1,  fxp2), am12);
        test_fxp("-10.25 - -2", fxp_sub(-fxp1, -fxp2), am8);
        test_fxp(" 10.25 *  2", fxp_mul( fxp1,  fxp2), a20);
        test_fxp(" 10.25 * -2", fxp_mul( fxp1, -fxp2), am20);
        test_fxp("-10.25 *  2", fxp_mul(-fxp1,  fxp2), am20);
        test_fxp("-10.25 * -2", fxp_mul(-fxp1, -fxp2), a20);
        test_fxp(" 10.25 /  2", fxp_div( fxp1,  fxp2), a5);
        test_fxp(" 10.25 / -2", fxp_div( fxp1, -fxp2), am5);
        test_fxp("-10.25 /  2", fxp_div(-fxp1,  fxp2), am5);
        test_fxp("-10.25 / -2", fxp_div(-fxp1, -fxp2), a5);

        printf("\nSome general tests\n");
        int fxp_zero     = fxp(0);
        int fxp_tiniest  = fxp_bin(0, 1);
        int fxp_mtiniest = -fxp_tiniest;
        int fxp_one      = fxp(1);
        int fxp_largest  = FXP_MAX;
        test_fxp("Infinity",  FXP_POS_INF, FXP_POS_INF);
        test_fxp("-Infinity", FXP_NEG_INF, FXP_NEG_INF);
        test_fxp("Undefined", FXP_UNDEF,   FXP_UNDEF);
        test_fxp("zero",     fxp_zero,      fxp_bin(0, 0));
        test_fxp("tiniest",  fxp_tiniest,   fxp_bin(0, 1));
        test_fxp("-tiniest", fxp_mtiniest,  -fxp_tiniest);
        test_fxp("one",      fxp_one,       fxp_dec(1, 0));
        test_fxp("Largest",  fxp_largest,
                             fxp_bin(FXP_WHOLE_MAX, FXP_FRAC_MAX-1));
        test_fxp("Largest - tiniest", fxp_sub(fxp_largest, fxp_tiniest),
                                      fxp_largest-1);
        test_fxp("Largest + tiniest safe",
                            fxp_sum(fxp_largest, fxp_tiniest),
                            FXP_POS_INF);
        test_fxp("Largest + tiniest unsafe",
                            fxp_unsafe_sum(fxp_largest, fxp_tiniest),
                            FXP_POS_INF);
        test_fxp("Largest + one safe",
                            fxp_sum(fxp_largest, fxp_one),
                            FXP_POS_INF);
        test_fxp("Largest + one unsafe",
                            fxp_unsafe_sum(fxp_largest, fxp_one),
                            fxp_largest + fxp_one);
        test_fxp("Way Too Large whole part!",
                            fxp(FXP_WHOLE_MAX*5),
                            FXP_POS_INF);
        test_fxp("Most negative",
                            -fxp_largest,
                            FXP_MIN);
        test_fxp("Almost most negative",
                            fxp_sum(-fxp_largest, fxp_tiniest),
                            fxp_bin(FXP_WHOLE_MIN, -FXP_FRAC_MAX+2));
        test_fxp("Safe Too neg substraction",
                            fxp_sub(-fxp_largest, fxp_one),
                            FXP_NEG_INF);
        test_fxp("Unsafe Too neg substraction",
                            fxp_unsafe_sub(-fxp_largest, fxp_one),
                            -fxp_largest - fxp_one);
        test_fxp("+inf + +inf",
                            fxp_sum(FXP_POS_INF, FXP_POS_INF),
                            FXP_POS_INF);
        test_fxp("-inf - +inf",
                            fxp_sub(FXP_NEG_INF, FXP_POS_INF),
                            FXP_NEG_INF);
        test_fxp("+inf - 1",
                            fxp_sub(FXP_POS_INF, fxp_one),
                            FXP_POS_INF);
        test_fxp("-inf + 1",
                            fxp_sum(FXP_NEG_INF, fxp_one),
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
                            fxp_div(fxp(10), fxp_zero),
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
                            fxp_div(fxp(-10), fxp_zero),
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

        printf("\nChecking sign taken from frac when whole==0\n");
        test_fxp("tiniest - 1",
                            fxp_sum(fxp_bin(0,1), fxp(-1)),
                            fxp_bin(0, -FXP_FRAC_MAX));
        test_fxp("-0.(+)500",
                            fxp_dec(-0, 500),
                            fxp_dec(0, 500));
        test_fxp("-0.(-)500",
                            fxp_dec(-0, -500),
                            fxp_dec(0, -500));

        printf("\nTruncation of longer frac decimal arguments\n");
        test_fxp("10.222222", fxp_dec(10, 222222), fxp_dec(10, 222));
        test_fxp("10.777777", fxp_dec(10, 777777), fxp_dec(10, 777));
        test_fxp("10.991999", fxp_dec(10, 991999), fxp_dec(10, 991));
        test_fxp("10.999999", fxp_dec(10, 999999), fxp_dec(10, 999));

        printf("\nFurther tests overflowing into infinities\n");
        test_fxp(" Max Fixed Point number", FXP_MAX, FXP_MAX);
        int halfmax = fxp_sub( fxp_div(FXP_MAX, fxp(2)), fxp_tiniest);
        test_fxp("~Half Max",
                            halfmax,
                            fxp_bin(FXP_WHOLE_MAX/2, FXP_FRAC_MAX-1));
        test_fxp("Half Max + 1000",
                            fxp_sum(halfmax, fxp(1000)),
                            fxp_bin(FXP_WHOLE_MAX/2+1000, FXP_FRAC_MAX-1));
        test_fxp("Half Max + Half Max",
                            fxp_sum(halfmax, halfmax),
                            halfmax+halfmax);
        test_fxp("FXP_MAX - Half Max",
                            fxp_sub(FXP_MAX, halfmax),
                            FXP_MAX/2 + fxp_tiniest);
        test_fxp("Half Max + FXP_MAX",
                            fxp_sum(halfmax, FXP_MAX),
                            FXP_POS_INF);
        test_fxp("-FXP_MAX - Half Max",
                            fxp_sub(-FXP_MAX, halfmax),
                            FXP_NEG_INF);
        test_fxp("Half Max * 2",
                            fxp_mul(halfmax, fxp(2)),
                            FXP_MAX - 2);
        test_fxp("Half Max * 2.001",
                            fxp_mul(halfmax, fxp_bin(2, 1)),
                            FXP_POS_INF);
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

        printf("\nTotal # of warnings: %d\n", nwarnings);
        printf("All tests passed using %d-bit and '%d' decimal fracs.\n",
                FXP_FRAC_BITS,
                FXP_FRAC_MAX_DEC);

        return 0;
}
