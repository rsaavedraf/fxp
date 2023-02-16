/*
 * fxp_xtimes.c
 * By Raul Saavedra, 2023-01-09
 *
 * Shows relavite execution times (smaller is better)
 * for the different fxp operations and functions.
 * Example part of the ouput, for a system with an
 * Intel i7-6700K cpu, run on 2023-02-15:
 *

Number of frac bits: 16
10%  20%  30%  40%  50%  60%  70%  80%  90%  100%
add      :   1.00
mul      :   8.90
mul_l    :   1.40
div      :  15.06
div_l    :   3.20
lg2      :  22.46  (BKM, only ints)
lg2_l    :  11.75  (about  0.52x lg2, using BKM and longs)
lg2_mul_l:  34.87  (about  1.55x lg2, using mult. and longs)

 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "fxp.h"
#include "fxp_aux.h"

#define DASHES "=================================================\n"
#define MAX_NUMS 1000
#define MAX_OPS  100

int main(void) {

        int s1, s2, s3, n1, n2, n3, x, y, n, lastp;
        double tadd, tadd_l, tmul, tmul_l, tmul_d;
        double tdiv, tdiv_l, tlg2, tlg2_l, tlg2_mul, tlg2_mul_l;
        clock_t t0, t1;
        int val[] = {FXP_UNDEF, FXP_MIN, \
                        FXP_MIN/3, -3333333, -300000, -30000, -3000, -300,
                        -200, -100, -50, -30, \
                        -20, -15, -10, -5, -3, -1,
                        0,
                        1, 3, 5, 10, 15, 20,
                        30, 50, 100, 200,
                        300, 3000, 30000, 300000, 3333333, FXP_MAX/3, FXP_MAX};
        int nvals = (int) (sizeof(val) / sizeof(int));
        int fracbit_configs[] = {8, 16, 24, 28};
        int nconfigs = sizeof(fracbit_configs) / sizeof(fracbit_configs[0]);

        srand((unsigned int) time(0));  // randomize seed

        printf("\n%sRelative Execution Times of FXP operations\n%s", DASHES, DASHES);
        print_sys_info();

        for (int nc = 0; nc < nconfigs; nc++) {
                int nfb = fracbit_configs[nc];
                fxp_set_frac_bits(nfb);
                printf("\nNumber of frac bits: %d\n", fxp_get_frac_bits());
                tadd = 0;
                tmul = 0;
                tmul_l = 0;
                tmul_d = 0;
                tdiv = 0;
                tdiv_l = 0;
                tlg2 = 0;
                tlg2_l = 0;
                tlg2_mul = 0;
                tlg2_mul_l = 0;
                lastp = -1;
                for (int n = 0; n < MAX_NUMS; n++) {
                        int p = (int) (((float) (n+1) * 100) / MAX_NUMS);
                        if (((p % 10) == 0) && (p != lastp) && (p != 0)) {
                                printf("%d%%  ", p);
                                fflush( stdout );
                                lastp = p;
                        }
                        s1 = rand() % 2 == 1? -1: 1;
                        s2 = rand() % 2 == 1? -1: 1;
                        s3 = rand() % 2 == 1? -1: 1;
                        n1 = s1 * rand();
                        n2 = s2 * rand();
                        // n3 is a number in (-1, 1)
                        n3 = s3 * (rand() % fxp_get_frac_max());
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
                        tadd += ((double) t1 - t0);

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
                        tmul += ((double) t1 - t0);

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
                        tmul_l += ((double) t1 - t0);

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
                        t1= clock();
                        tdiv += ((double) t1 - t0);

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
                        t1= clock();
                        tdiv_l += ((double) t1 - t0);

                        // To measure the lg execution use only + arguments
                        if (n1 < 0) n1 = -n1;
                        if (n2 < 0) n2 = -n2;
                        if (n3 < 0) n3 = -n3;

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
                        t1= clock();
                        tlg2 += ((double) t1 - t0);

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
                        t1= clock();
                        tlg2_l += ((double) t1 - t0);

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
                        t1= clock();
                        tlg2_mul_l += ((double) t1 - t0);

                }

                printf("\nadd      : %6.2lf\n", 1.0);
                printf("mul      : %6.2lf\n", tmul / tadd);
                printf("mul_l    : %6.2lf\n", tmul_l / tadd);
                printf("div      : %6.2lf\n", tdiv / tadd);
                printf("div_l    : %6.2lf\n", tdiv_l / tadd);
                printf("lg2      : %6.2lf  (BKM, only ints)\n", \
                            tlg2 / tadd);
                printf("lg2_l    : %6.2lf  (about %5.2lfx lg2, using BKM and longs)\n", \
                            tlg2_l / tadd, tlg2_l / tlg2);
                printf("lg2_mul_l: %6.2lf  (about %5.2lfx lg2, using mult. and longs)\n", \
                             tlg2_mul_l / tadd, tlg2_mul_l / tlg2);

        }
        return 0;
}
