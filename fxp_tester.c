/*
 * fxp_tester.c
 * Tests the implementation of binary fixed point (fxp.c)
 * By Raul Saavedra, 2022-11-13
 */

#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include "fxp.h"

void print_fxp(int fxp)
{
        if (fxp == FXP_POS_INF || fxp == FXP_NEG_INF || fxp == FXP_UNDEF)
                printf(fxp==FXP_UNDEF? "UNDEF":
                      (fxp==FXP_POS_INF? "+INF": "-INF"));
        else
                printf("%d.%d (bin frac=%d)",
                        fxp_get_whole_part(fxp),
                        fxp_get_dec_frac(fxp),
                        fxp_get_bin_frac(fxp));
}

// Check if equal to the assert value
void test_fxp(char *s, int fxp, int assert_val)
{
        printf("For %s got ", s);
        print_fxp(fxp);
        printf(",\texpected ");
        print_fxp(assert_val);
        printf("\n");
        fflush( stdout );
        // Allowing 1 bit error (least significant bit in frac part)
        if (fxp >= assert_val)
            assert( fxp - assert_val <= 1);
        else
            assert( assert_val -fxp <= 1);
}


int main(void)
{
        printf("==============\nFXP Tester run\n==============\n");

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
        printf("max fraction: %d\n", FXP_FRAC_MAX);
        printf("max whole   : %d\n", FXP_WHOLE_MAX);
        printf("min whole   : %d\n", FXP_WHOLE_MIN);
        printf("pos infinity: %d\n", FXP_POS_INF);
        printf("neg infinity: %d\n", FXP_NEG_INF);
        printf("undefined   : %d\n", FXP_UNDEF);

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
        test_fxp("Almost most negative",
                            fxp_sum(-fxp_largest, fxp_tiniest),
                            fxp_bin(FXP_WHOLE_MIN, -FXP_FRAC_MAX+1));
        test_fxp("Most negative",
                            -fxp_largest,
                            FXP_MIN);
        test_fxp("Safe Too negative substraction",
                            fxp_sub(-fxp_largest, fxp_one),
                            FXP_NEG_INF);
        test_fxp("Unsafe Too negative substraction",
                            fxp_unsafe_sub(-fxp_largest, fxp_one),
                            -fxp_largest - fxp_one);
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

        printf("\nTrimming of large frac decimal arguments\n");
        print_fxp(fxp_dec(1000, 777)); printf("\n");
        print_fxp(fxp_dec(1000, 778)); printf("\n");
        test_fxp("1000.777777",
                            fxp_dec(1000, 777777),
                            fxp_dec(1000, 778));
        test_fxp("2000.222222",
                            fxp_dec(2000, 222222),
                            fxp_dec(2000, 222));
        test_fxp("9.991999 is:",
                            fxp_dec(9, 991999),
                            fxp_dec(9, 992));

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
                            fxp_bin(FXP_WHOLE_MAX, FXP_FRAC_MAX-2));
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

        printf("\nAll fxp tests passed for %d frac bits!\n", FXP_FRAC_BITS);

        return 0;
}
