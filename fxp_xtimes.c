/* SPDX-License-Identifier: MIT */
/*
 * fxp_xtimes.c
 * By Raul Saavedra, 2023-01-09
 *
 * Shows relavite execution times (smaller is better)
 * for the different fxp operations and functions.
 * Example output: file xtimes.txt in the repo
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "fxp.h"
#include "fxp_l.h"
#include "fxp_aux.h"

#define DASHES "=================================================\n"
#define MAX_NUMS 250
#define MAX_OPS  100

const long double PI_AS_LD = acosl(-1.0L);

// Returns a random angle in range (-pi, pi)
long double get_random_angle()
{
        return (PI_AS_LD * ((rand() % 2)? -1: 1)) \
                * (((long double) rand()) / INT_MAX);
}


int main(void) {

        int s1, s2, s3, n1, n2, n3, x, y, n, lastp;

        long double tadd, tadd_l, tmul, tmul_l, tmul_d, tdiv, tdiv_l;
        long double tlg2, tlg2_l, tlg2_mul, tlg2_mul_l, tln, tln_l;
        long double tpow2, tpow2_l, texp, texp_l, tpow10, tpow10_l;
        long double tlg10, tlg10_l, tsqrt, tsqrt_l, tpowxy, tpowxy_l;
        long double tcos_l, tsin_l, tcossin_l;
        long double avgadd, avgadd_l, avgmul, avgmul_l, avgmul_d;
        long double avgdiv, avgdiv_l, avglg2, avglg2_l, avglg2_mul;
        long double avglg2_mul_l, avgln, avgln_l, avglg10, avglg10_l;
        long double avgpow2, avgpow2_l, avgexp, avgexp_l;
        long double avgpow10, avgpow10_l, avgsqrt, avgsqrt_l;
        long double avgpowxy, avgpowxy_l;
        long double avgcos_l, avgsin_l, avgcossin_l;
        // Times for the system's native operations
        long double tadd_sys, tmul_sys, tdiv_sys;
        long double avgadd_sys, avgmul_sys, avgdiv_sys;

        clock_t t0, t1, dt;
        int val[] = {FXP_UNDEF, FXP_MIN, \
                        FXP_MIN/3, -3333333, -300000, -30000, -3000, -300,
                        -200, -100, -50, -30, \
                        -20, -15, -10, -5, -3, -1,
                        0,
                        1, 3, 5, 10, 15, 20,
                        30, 50, 100, 200,
                        300, 3000, 30000, 300000, 3333333, FXP_MAX/3, FXP_MAX};
        int nvals = (int) (sizeof(val) / sizeof(int));
        int fracbit_configs[] = {8, 12, 16, 20, 24, 28};
        int nconfigs = sizeof(fracbit_configs) / sizeof(fracbit_configs[0]);

        srand((unsigned int) time(0));  // randomize seed

        printf("\n%sRelative Execution Times of FXP operations\n%s", DASHES, DASHES);
        print_sys_info();

        avgadd = 0.0;
        avgmul = 0.0;
        avgmul_l = 0.0;
        avgmul_d = 0.0;
        avgdiv = 0.0;
        avgdiv_l = 0.0;
        avglg2 = 0.0;
        avglg2_l = 0.0;
        avglg2_mul = 0.0;
        avglg2_mul_l = 0.0;
        avgln = 0.0;
        avgln_l = 0.0;
        avglg10 = 0.0;
        avglg10_l = 0.0;
        avgpow2 = 0.0;
        avgpow2_l = 0.0;
        avgexp = 0.0;
        avgexp_l = 0.0;
        avgpow10 = 0.0;
        avgpow10_l = 0.0;
        avgpowxy = 0.0;
        avgpowxy_l = 0.0;
        avgsqrt = 0.0;
        avgsqrt_l = 0.0;
        avgadd_sys = 0.0;
        avgmul_sys = 0.0;
        avgdiv_sys = 0.0;
        avgcos_l = 0.0;
        avgsin_l = 0.0;
        avgcossin_l = 0.0;
        for (int nc = 0; nc < nconfigs; nc++) {
                int nfb = fracbit_configs[nc];
                fxp_set_frac_bits(nfb);
                printf("\nNumber of frac bits: %d\n", FXP_frac_bits);
                tadd = 0.0;
                tmul = 0.0;
                tmul_l = 0.0;
                tmul_d = 0.0;
                tdiv = 0.0;
                tdiv_l = 0.0;
                tlg2 = 0.0;
                tlg2_l = 0.0;
                tlg2_mul = 0.0;
                tlg2_mul_l = 0.0;
                tln = 0.0;
                tln_l = 0.0;
                tlg10 = 0.0;
                tlg10_l = 0.0;
                tpow2 = 0.0;
                tpow2_l = 0.0;
                texp = 0.0;
                texp_l = 0.0;
                tpow10 = 0.0;
                tpow10_l = 0.0;
                tpowxy = 0.0;
                tpowxy_l = 0.0;
                tsqrt = 0.0;
                tsqrt_l = 0.0;
                tcos_l = 0.0;
                tsin_l = 0.0;
                tcossin_l = 0.0;
                lastp = -1;
                tadd_sys = 0.0;
                tmul_sys = 0.0;
                tdiv_sys = 0.0;
                int nwb = FXP_INT_BITS - nfb;
                unsigned int mask_for_pow = (1u << ((nwb > 6? 6: nwb) + nfb - 1)) - 1;
                for (int n = 0; n < MAX_NUMS; n++) {
                        int p = (int) (((float) (n+1) * 100) / MAX_NUMS);
                        if (((p % 10) == 0) && (p != lastp) && (p != 0)) {
                                printf("%d%%  ", p);
                                fflush( stdout );
                                lastp = p;
                        }
                        // random signs s1, s2, s3
                        s1 = (rand() % 2)? -1: 1;
                        s2 = (rand() % 2)? -1: 1;
                        s3 = (rand() % 2)? -1: 1;
                        // random ints n1, n2
                        n1 = s1 * rand();
                        n2 = s2 * rand();
                        // n3 is a random number in (-1, 1)
                        n3 = s3 * (rand() % FXP_frac_max);

                        // System's native sum of ints
                        t0 = clock();
                        for (int i = 0; i < MAX_OPS; i++) {
                                x = n1 + n2;
                                x = n3 + n2;
                                for (int j = 0; j < nvals; j++) {
                                        y = val[j];
                                        x = n1 + y;
                                        x = n2 + y;
                                        x = n3 + y;
                                }
                        }
                        t1 = clock();
                        dt = ((double) t1 - t0);
                        tadd_sys += dt;
                        avgadd_sys += dt;

                        // Sums of fxp's
                        t0 = clock();
                        for (int i = 0; i < MAX_OPS; i++) {
                                x = fxp_add(n1, n2);
                                x = fxp_add(n3, n2);
                                for (int j = 0; j < nvals; j++) {
                                        y = val[j];
                                        x = fxp_add(n1, y);
                                        x = fxp_add(n2, y);
                                        x = fxp_add(n3, y);
                                }
                        }
                        t1 = clock();
                        dt = ((double) t1 - t0);
                        tadd += dt;
                        avgadd += dt;

                        // System's native multiplication of ints
                        t0 = clock();
                        for (int i = 0; i < MAX_OPS; i++) {
                                x = n1 * n2;
                                x = n3 * n2;
                                for (int j = 0; j < nvals; j++) {
                                        y = val[j];
                                        x = n1 * y;
                                        x = n2 * y;
                                        x = n3 * y;
                                }
                        }
                        t1 = clock();
                        dt = ((double) t1 - t0);
                        tmul_sys += dt;
                        avgmul_sys += dt;

                        // Multiplication of fxp's
                        t0 = clock();
                        for (int i = 0; i < MAX_OPS; i++) {
                                x = fxp_mul(n1, n2);
                                x = fxp_mul(n3, n2);
                                for (int j = 0; j < nvals; j++) {
                                        y = val[j];
                                        x = fxp_mul(n1, y);
                                        x = fxp_mul(n2, y);
                                        x = fxp_mul(n3, y);
                                }
                        }
                        t1 = clock();
                        dt = ((double) t1 - t0);
                        tmul += dt;
                        avgmul += dt;

                        // Multiplication of fxp's using longs
                        t0 = clock();
                        for (int i = 0; i < MAX_OPS; i++) {
                                x = fxp_mul_l(n1, n2);
                                x = fxp_mul_l(n3, n2);
                                for (int j = 0; j < nvals; j++) {
                                        y = val[j];
                                        x = fxp_mul_l(n1, y);
                                        x = fxp_mul_l(n2, y);
                                        x = fxp_mul_l(n3, y);
                                }
                        }
                        t1 = clock();
                        dt = ((double) t1 - t0);
                        tmul_l += dt;
                        avgmul_l += dt;

                        // Temporarily replace 0 with 1 in val[] for divisions
                        val[18] = 1;

                        // System's native division of ints
                        t0 = clock();
                        for (int i = 0; i < MAX_OPS; i++) {
                                x = n1 / n2;
                                x = n3 / n2;
                                for (int j = 0; j < nvals; j++) {
                                        y = val[j];
                                        x = n1 / y;
                                        x = n2 / y;
                                        x = n3 / y;
                                }
                        }
                        t1 = clock();
                        dt = ((double) t1 - t0);
                        tdiv_sys += dt;
                        avgdiv_sys += dt;

                        // Division of fxp's
                        t0 = clock();
                        for (int i = 0; i < MAX_OPS; i++) {
                                x = fxp_div(n1, n2);
                                x = fxp_div(n3, n2);
                                for (int j = 0; j < nvals; j++) {
                                        y = val[j];
                                        x = fxp_div(n1, y);
                                        x = fxp_div(n2, y);
                                        x = fxp_div(n3, y);
                                }
                        }
                        t1 = clock();
                        dt = ((double) t1 - t0);
                        tdiv += dt;
                        avgdiv += dt;

                        // Division of fxp's using longs
                        t0 = clock();
                        for (int i = 0; i < MAX_OPS; i++) {
                                x = fxp_div_l(n1, n2);
                                x = fxp_div_l(n3, n2);
                                for (int j = 0; j < nvals; j++) {
                                        y = val[j];
                                        x = fxp_div_l(n1, y);
                                        x = fxp_div_l(n2, y);
                                        x = fxp_div_l(n3, y);
                                }
                        }
                        t1 = clock();
                        dt = ((double) t1 - t0);
                        tdiv_l += dt;
                        avgdiv_l += dt;

                        // Divisions done, restore val[18] from 1 to 0
                        val[18] = 0;

                        // Calculation of pow2 of fxp's using
                        // BKM and only ints
                        n1 = n1 & mask_for_pow;
                        n2 = -n3;
                        t0 = clock();
                        for (int i = 0; i < MAX_OPS; i++) {
                                x = fxp_pow2(n1);
                                x = fxp_pow2(n2);
                                for (int j = 0; j < nvals; j++) {
                                        y = val[j];
                                        x = fxp_pow2(n1);
                                        x = fxp_pow2(n2);
                                        x = fxp_pow2(n3);
                                }
                        }
                        t1 = clock();
                        dt = ((double) t1 - t0);
                        tpow2 += dt;
                        avgpow2 += dt;

                        // Calculation of pow2_l of fxp's using
                        // BKM and longs
                        t0 = clock();
                        for (int i = 0; i < MAX_OPS; i++) {
                                x = fxp_pow2_l(n1);
                                x = fxp_pow2_l(n2);
                                for (int j = 0; j < nvals; j++) {
                                        y = val[j];
                                        x = fxp_pow2_l(n1);
                                        x = fxp_pow2_l(n2);
                                        x = fxp_pow2_l(n3);
                                }
                        }
                        t1 = clock();
                        dt = ((double) t1 - t0);
                        tpow2_l += dt;
                        avgpow2_l += dt;

                        // Calculation of exp of fxp's using
                        // pow2 (-> only ints)
                        n1 = n1 & mask_for_pow;
                        n2 = -n3;
                        t0 = clock();
                        for (int i = 0; i < MAX_OPS; i++) {
                                x = fxp_exp(n1);
                                x = fxp_exp(n2);
                                for (int j = 0; j < nvals; j++) {
                                        y = val[j];
                                        x = fxp_exp(n1);
                                        x = fxp_exp(n2);
                                        x = fxp_exp(n3);
                                }
                        }
                        t1 = clock();
                        dt = ((double) t1 - t0);
                        texp += dt;
                        avgexp += dt;

                        // Calculation of exp_l of fxp's using
                        // pow2_l (-> longs)
                        t0 = clock();
                        for (int i = 0; i < MAX_OPS; i++) {
                                x = fxp_exp_l(n1);
                                x = fxp_exp_l(n2);
                                for (int j = 0; j < nvals; j++) {
                                        y = val[j];
                                        x = fxp_exp_l(n1);
                                        x = fxp_exp_l(n2);
                                        x = fxp_exp_l(n3);
                                }
                        }
                        t1 = clock();
                        dt = ((double) t1 - t0);
                        texp_l += dt;
                        avgexp_l += dt;

                        // Calculation of pow10 of fxp's using
                        // pow2() (only ints)
                        t0 = clock();
                        for (int i = 0; i < MAX_OPS; i++) {
                                x = fxp_pow10(n1);
                                x = fxp_pow10(n2);
                                for (int j = 0; j < nvals; j++) {
                                        y = val[j];
                                        x = fxp_pow10(n1);
                                        x = fxp_pow10(n2);
                                        x = fxp_pow10(n3);
                                }
                        }
                        t1 = clock();
                        dt = ((double) t1 - t0);
                        tpow10 += dt;
                        avgpow10 += dt;

                        // Calculation of pow10_l of fxp's using
                        // pow2_l() (-> longs)
                        t0 = clock();
                        for (int i = 0; i < MAX_OPS; i++) {
                                x = fxp_pow10_l(n1);
                                x = fxp_pow10_l(n2);
                                for (int j = 0; j < nvals; j++) {
                                        y = val[j];
                                        x = fxp_pow10_l(n1);
                                        x = fxp_pow10_l(n2);
                                        x = fxp_pow10_l(n3);
                                }
                        }
                        t1 = clock();
                        dt = ((double) t1 - t0);
                        tpow10_l += dt;
                        avgpow10_l += dt;

                        // To measure the lg execution use only + arguments
                        if (n1 < 0) n1 = -n1;
                        if (n2 < 0) n2 = -n2;
                        if (n3 < 0) n3 = -n3;

                        // Calculation of lg2 of fxp's using only ints
                        t0 = clock();
                        for (int i = 0; i < MAX_OPS; i++) {
                                x = fxp_lg2(n1);
                                x = fxp_lg2(n2);
                                for (int j = 0; j < nvals; j++) {
                                        y = val[j];
                                        x = fxp_lg2(n1);
                                        x = fxp_lg2(n2);
                                        x = fxp_lg2(n3);
                                }
                        }
                        t1 = clock();
                        dt = ((double) t1 - t0);
                        tlg2 += dt;
                        avglg2 += dt;

                        // Calculation of lg2 of fxp's using longs
                        t0 = clock();
                        for (int i = 0; i < MAX_OPS; i++) {
                                x = fxp_lg2_l(n1);
                                x = fxp_lg2_l(n2);
                                for (int j = 0; j < nvals; j++) {
                                        y = val[j];
                                        x = fxp_lg2_l(n1);
                                        x = fxp_lg2_l(n2);
                                        x = fxp_lg2_l(n3);
                                }
                        }
                        t1 = clock();
                        dt = ((double) t1 - t0);
                        tlg2_l += dt;
                        avglg2_l += dt;

                        // Calculation of lg2 of fxp's using the
                        // multiplication-based algorithm, and longs
                        t0 = clock();
                        for (int i = 0; i < MAX_OPS; i++) {
                                x = fxp_lg2_mul_l(n1);
                                x = fxp_lg2_mul_l(n2);
                                for (int j = 0; j < nvals; j++) {
                                        y = val[j];
                                        x = fxp_lg2_mul_l(n1);
                                        x = fxp_lg2_mul_l(n2);
                                        x = fxp_lg2_mul_l(n3);
                                }
                        }
                        t1 = clock();
                        dt = ((double) t1 - t0);
                        tlg2_mul_l += dt;
                        avglg2_mul_l += dt;

                        // Calculation of ln of fxp's using ints
                        t0 = clock();
                        for (int i = 0; i < MAX_OPS; i++) {
                                x = fxp_ln(n1);
                                x = fxp_ln(n2);
                                for (int j = 0; j < nvals; j++) {
                                        y = val[j];
                                        x = fxp_ln(n1);
                                        x = fxp_ln(n2);
                                        x = fxp_ln(n3);
                                }
                        }
                        t1 = clock();
                        dt = ((double) t1 - t0);
                        tln += dt;
                        avgln += dt;

                        // Calculation of ln of fxp's using longs
                        t0 = clock();
                        for (int i = 0; i < MAX_OPS; i++) {
                                x = fxp_ln_l(n1);
                                x = fxp_ln_l(n2);
                                for (int j = 0; j < nvals; j++) {
                                        y = val[j];
                                        x = fxp_ln_l(n1);
                                        x = fxp_ln_l(n2);
                                        x = fxp_ln_l(n3);
                                }
                        }
                        t1 = clock();
                        dt = ((double) t1 - t0);
                        tln_l += dt;
                        avgln_l += dt;

                        // Calculation of lg10 of fxp's
                        t0 = clock();
                        for (int i = 0; i < MAX_OPS; i++) {
                                x = fxp_lg10(n1);
                                x = fxp_lg10(n2);
                                for (int j = 0; j < nvals; j++) {
                                        y = val[j];
                                        x = fxp_lg10(n1);
                                        x = fxp_lg10(n2);
                                        x = fxp_lg10(n3);
                                }
                        }
                        t1 = clock();
                        dt = ((double) t1 - t0);
                        tlg10 += dt;
                        avglg10 += dt;

                        // Calculation of lg10 of fxp's using longs
                        t0 = clock();
                        for (int i = 0; i < MAX_OPS; i++) {
                                x = fxp_lg10_l(n1);
                                x = fxp_lg10_l(n2);
                                for (int j = 0; j < nvals; j++) {
                                        y = val[j];
                                        x = fxp_lg10_l(n1);
                                        x = fxp_lg10_l(n2);
                                        x = fxp_lg10_l(n3);
                                }
                        }
                        t1 = clock();
                        dt = ((double) t1 - t0);
                        tlg10_l += dt;
                        avglg10_l += dt;

                        // Calculation of sqrt_l of fxp's using longs
                        t0 = clock();
                        for (int i = 0; i < MAX_OPS; i++) {
                                x = fxp_sqrt_l(n1);
                                x = fxp_sqrt_l(n2);
                                for (int j = 0; j < nvals; j++) {
                                        y = val[j];
                                        x = fxp_sqrt_l(n1);
                                        x = fxp_sqrt_l(n2);
                                        x = fxp_sqrt_l(n3);
                                }
                        }
                        t1 = clock();
                        dt = ((double) t1 - t0);
                        tsqrt_l += dt;
                        avgsqrt_l += dt;

                        // Calculation of sqrt of fxp's
                        t0 = clock();
                        for (int i = 0; i < MAX_OPS; i++) {
                                x = fxp_sqrt(n1);
                                x = fxp_sqrt(n2);
                                for (int j = 0; j < nvals; j++) {
                                        y = val[j];
                                        x = fxp_sqrt(n1);
                                        x = fxp_sqrt(n2);
                                        x = fxp_sqrt(n3);
                                }
                        }
                        t1 = clock();
                        dt = ((double) t1 - t0);
                        tsqrt += dt;
                        avgsqrt += dt;

                        // Calculation of powxy of fxp's using
                        // lg2 and pow2
                        // make n3 negative to enforce the calculations
                        n3 = -n3;
                        t0 = clock();
                        for (int i = 0; i < MAX_OPS; i++) {
                                x = fxp_powxy(n1, n3);
                                x = fxp_powxy(n2, n3);
                                for (int j = 0; j < nvals; j++) {
                                        y = val[j];
                                        x = fxp_powxy(n1, n3);
                                        x = fxp_powxy(n2, n3);
                                        x = fxp_powxy(n1, n3);
                                }
                        }
                        t1 = clock();
                        dt = ((double) t1 - t0);
                        tpowxy += dt;
                        avgpowxy += dt;

                        // Calculation of powxy_l of fxp's using
                        // lg2_l and pow2_l (-> longs)
                        t0 = clock();
                        for (int i = 0; i < MAX_OPS; i++) {
                                x = fxp_powxy_l(n1, n3);
                                x = fxp_powxy_l(n2, n3);
                                for (int j = 0; j < nvals; j++) {
                                        y = val[j];
                                        x = fxp_powxy_l(n1, n3);
                                        x = fxp_powxy_l(n2, n3);
                                        x = fxp_powxy_l(n1, n3);
                                }
                        }
                        t1 = clock();
                        dt = ((double) t1 - t0);
                        tpowxy_l += dt;
                        avgpowxy_l += dt;

                        // Trigonometrics
                        long double angle1 = get_random_angle();
                        long double angle2 = get_random_angle();
                        long double angle3 = get_random_angle();
                        int n1 = ld2fxp(angle1);
                        int n2 = ld2fxp(angle2);
                        int n3 = ld2fxp(angle3);

                        // Calculation of cos and sin of fxp's using
                        // cossin_l (-> longs)
                        fxptuple tup;
                        t0 = clock();
                        for (int i = 0; i < MAX_OPS; i++) {
                                tup = fxp_cos_sin_l(n1);
                                tup = fxp_cos_sin_l(n2);
                                for (int j = 0; j < nvals; j++) {
                                        y = val[j];
                                        tup = fxp_cos_sin_l(n1);
                                        tup = fxp_cos_sin_l(n2);
                                        tup = fxp_cos_sin_l(n3);
                                }
                        }
                        t1 = clock();
                        dt = ((double) t1 - t0);
                        tcossin_l += dt;
                        avgcossin_l += dt;
                }

                // Results for this configuration of frag bits
                printf("\nUsing only ints:\n");
                printf("\tadd      : %6.2lf\n", 1.0);
                printf("\tmul      : %6.2Lf\n", tmul / tadd);
                printf("\tdiv      : %6.2Lf\n", tdiv / tadd);
                printf("\tlg2      : %6.2Lf  (BKM-L, only ints)\n", \
                            tlg2 / tadd);
                printf("\tln       : %6.2Lf  (using lg2)\n", \
                            tln / tadd);
                printf("\tlg10     : %6.2Lf  (using lg2)\n", \
                            tlg10 / tadd);
                printf("\tpow2     : %6.2Lf  (BKM-E, only ints)\n", \
                            tpow2 / tadd);
                printf("\texp      : %6.2Lf  (about %5.2Lfx pow2, using pow2)\n", \
                            texp / tadd, texp / tpow2);
                printf("\tpow10    : %6.2Lf  (about %5.2Lfx pow2, using pow2)\n", \
                            tpow10 / tadd, tpow10 / tpow2);
                printf("\tsqrt     : %6.2Lf  (about %5.2Lfx pow2, using lg2 & pow2)\n", \
                            tsqrt / tadd, tsqrt / tpow2);
                printf("\tpowxy    : %6.2Lf  (about %5.2Lfx pow2, using lg2 & pow2)\n", \
                            tpowxy / tadd, tpowxy / tpow2);
                printf("Using longs:\n");
                printf("\tmul_l    : %6.2Lf\n", tmul_l / tadd);
                printf("\tdiv_l    : %6.2Lf\n", tdiv_l / tadd);
                printf("\tlg2_l    : %6.2Lf  (about %5.2Lfx lg2, using BKM-L and longs)\n", \
                            tlg2_l / tadd, tlg2_l / tlg2);
                printf("\tlg2_mul_l: %6.2Lf  (about %5.2Lfx lg2, using mult. and longs)\n", \
                            tlg2_mul_l / tadd, tlg2_mul_l / tlg2);
                printf("\tln_l     : %6.2Lf  (about %5.2Lfx lg2, using lg2_l)\n", \
                            tln_l / tadd, tln_l / tlg2);
                printf("\tlg10_l   : %6.2Lf  (about %5.2Lfx lg2, using lg2_l)\n", \
                            tlg10_l / tadd, tlg10_l / tlg2);
                printf("\tpow2_l   : %6.2Lf  (about %5.2Lfx pow2, using BKM-E and longs)\n", \
                            tpow2_l / tadd, tpow2_l / tpow2);
                printf("\texp_l    : %6.2Lf  (about %5.2Lfx pow2, using pow2_l)\n", \
                            texp_l / tadd, texp_l / tpow2);
                printf("\tpow10_l  : %6.2Lf  (about %5.2Lfx pow2, using pow2_l)\n", \
                            tpow10_l / tadd, tpow10_l / tpow2);
                printf("\tsqrt_l   : %6.2Lf  (about %5.2Lfx pow2, using lg2_l & pow2_l)\n", \
                            tsqrt_l / tadd, tsqrt_l / tpow2);
                printf("\tpowxy_l  : %6.2Lf  (about %5.2Lfx pow2, using lg2_l & pow2_l)\n", \
                            tpowxy_l / tadd, tpowxy_l / tpow2);
                printf("\tcossin_l : %6.2Lf\n", tcossin_l / tadd);
        }

        // Overall results
        printf("\n\n%sRelative Xtime averages for frac bit configurations {", DASHES);
        for (int nc = 0; nc < nconfigs; nc++) {
                printf("%d", fracbit_configs[nc]);
                printf((nc + 1 == nconfigs? "}": ", "));
        }
        printf("\n%s", DASHES);
        avgadd_sys /= nconfigs;
        avgmul_sys /= nconfigs;
        avgdiv_sys /= nconfigs;
        avgadd /= nconfigs;
        avgmul /= nconfigs;
        avgmul_l /= nconfigs;
        avgdiv /= nconfigs;
        avgdiv_l /= nconfigs;
        avglg2 /= nconfigs;
        avglg2_l /= nconfigs;
        avglg10 /= nconfigs;
        avglg10_l /= nconfigs;
        avgln /= nconfigs;
        avgln_l /= nconfigs;
        avglg2_mul_l /= nconfigs;
        avgpow2 /= nconfigs;
        avgpow2_l /= nconfigs;
        avgpow10 /= nconfigs;
        avgpow10_l /= nconfigs;
        avgexp /= nconfigs;
        avgexp_l /= nconfigs;
        avgsqrt /= nconfigs;
        avgsqrt_l /= nconfigs;
        avgpowxy /= nconfigs;
        avgpowxy_l /= nconfigs;
        avgcossin_l /= nconfigs;
        printf("\nUsing only ints:\n");
        printf("\tadd      : %6.2Lf  (%6.2Lfx system's native addition of ints)\n", \
                    1.0L, avgadd / avgadd_sys);
        printf("\tmul      : %6.2Lf  (%6.2Lfx system's native multiplication of ints)\n", \
                    avgmul / avgadd, avgmul / avgmul_sys);
        printf("\tdiv      : %6.2Lf  (%6.2Lfx system's native division of ints)\n", \
                    avgdiv / avgadd, avgdiv / avgdiv_sys);
        printf("\tlg2      : %6.2Lf  (BKM, only ints)\n", \
                    avglg2 / avgadd);
        printf("\tln       : %6.2Lf  (about %5.2Lfx lg2, using lg2)\n", \
                    avgln / avgadd, avgln / avglg2);
        printf("\tlg10     : %6.2Lf  (about %5.2Lfx lg2, using lg2)\n", \
                    avglg10 / avgadd, avglg10 / avglg2);
        printf("\tpow2     : %6.2Lf  (BKM, only ints)\n", \
                    avgpow2 / avgadd);
        printf("\texp      : %6.2Lf  (about %5.2Lfx pow2, using pow2)\n", \
                    avgexp / avgadd, \
                    avgexp / avgpow2);
        printf("\tpow10    : %6.2Lf  (about %5.2Lfx pow2, using pow2)\n", \
                    avgpow10 / avgadd, \
                    avgpow10 / avgpow2);
        printf("\tsqrt     : %6.2Lf  (about %5.2Lfx pow2, using lg2 & pow2)\n", \
                    avgsqrt / avgadd, \
                    avgsqrt / avgpow2);
        printf("\tpowxy    : %6.2Lf  (about %5.2Lfx pow2, using lg2 & pow2)\n", \
                    avgpowxy / avgadd, \
                    avgpowxy / avgpow2);

        printf("Using longs:\n");
        printf("\tmul_l    : %6.2Lf  (%6.2Lfx system's native multiplication of ints)\n", \
                    avgmul_l / avgadd, avgmul_l / avgmul_sys);
        printf("\tdiv_l    : %6.2Lf  (%6.2Lfx system's native division of ints)\n", \
                    avgdiv_l / avgadd, avgdiv_l / avgdiv_sys);
        printf("\tlg2_l    : %6.2Lf  (about %5.2Lfx lg2, using BKM and longs)\n", \
                    avglg2_l / avgadd, avglg2_l / avglg2);
        printf("\tlg2_mul_l: %6.2Lf  (about %5.2Lfx lg2, using mult. and longs)\n", \
                    avglg2_mul_l / avgadd, \
                    avglg2_mul_l / avglg2);
        printf("\tln_l     : %6.2Lf  (about %5.2Lfx lg2, using lg2_l)\n", \
                    avgln_l / avgadd, avgln_l / avglg2);
        printf("\tlg10_l   : %6.2Lf  (about %5.2Lfx lg2, using lg2_l)\n", \
                    avglg10_l / avgadd, avglg10_l / avglg2);
        printf("\tpow2_l   : %6.2Lf  (about %5.2Lfx pow2, using BKM and longs)\n", \
                    avgpow2_l / avgadd, \
                    avgpow2_l / avgpow2);
        printf("\texp_l    : %6.2Lf  (about %5.2Lfx pow2, using pow2_l)\n", \
                    avgexp_l / avgadd, \
                    avgexp_l / avgpow2);
        printf("\tpow10_l  : %6.2Lf  (about %5.2Lfx pow2, using pow2_l)\n", \
                    avgpow10_l / avgadd, \
                    avgpow10_l / avgpow2);
        printf("\tsqrt_l   : %6.2Lf  (about %5.2Lfx pow2, using lg2_l & pow2_l)\n", \
                    avgsqrt_l / avgadd, \
                    avgsqrt_l / avgpow2);
        printf("\tpowxy_l  : %6.2Lf  (about %5.2Lfx pow2, using lg2_l & pow2_l)\n", \
                    avgpowxy_l / avgadd, \
                    avgpowxy_l / avgpow2);
        printf("\tcossin_l : %6.2Lf\n", \
                    avgcossin_l / avgadd);

        printf("%s", DASHES);
        printf("(Keep in mind: compiler optimization options used/not used can ");
        printf("affect these measurements significantly.)\n\n");

        return 0;
}
