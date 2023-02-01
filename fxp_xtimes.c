/*
 * fxp_xtimes.c
 * By Raul Saavedra, 2023-01-09
 *
 * Shows relavite execution times (smaller is better)
 * for the different fxp arithmetic operations.
 * Example end of the ouput, for a system with an
 * Intel i7-6700K cpu, run on 2023-01-31:
 *
 * add:    1.000000
 * add_l:  1.009193
 * mul:    8.523299
 * mul_l:  1.487134
 * div:    20.434878
 * div_l:  3.249009
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "fxp.h"

#define DASHES "================================================\n"
#define MAX_NUMS 10000
#define MAX_OPS  500

int main(void) {

        int s1, s2, s3, n1, n2, n3, x, y, n, lastp;
        double tadd, tadd_l, tmul, tmul_l, tmul_d, tdiv, tdiv_l;
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

        printf("%sRelative Execution Times of FXP operations\n%s", DASHES, DASHES);
        srand((unsigned int) time(0));  // randomize seed
        printf("MAX INT is %d\n", FXP_MAX);

        tadd = 0;
        tadd_l = 0;
        tmul = 0;
        tmul_l = 0;
        tmul_d = 0;
        tdiv = 0;
        tdiv_l = 0;
        lastp = -1;
        for (int n = 0; n < MAX_NUMS; n++) {
                int p = (int) (((float) (n+1) * 100) / MAX_NUMS);
                if (((p % 10) == 0) && (p != lastp)) {
                    printf("%3d %%\n", p);
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
                        x = fxp_add_l(n1, n2);
                        x = fxp_add_l(n3, n2);
                        for (int j = 0; j < nvals; j++) {
                            y = val[j];
                            x = fxp_add_l(n1, y);
                            x = fxp_add_l(n2, y);
                            x = fxp_add_l(n3, y);
                        }
                }
                t1 = clock();
                tadd_l += ((double) t1 - t0);

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

                /*
                t0 = clock();
                for (int i = 0; i < MAX_OPS; i++) {
                        x = fxp_mul_d(n1, n2);
                        x = fxp_mul_d(n3, n2);
                        for (int j = 0; j < nvals; j++) {
                            y = val[j];
                            x = fxp_mul_d(n1, y);
                            x = fxp_mul_d(n2, y);
                            x = fxp_mul_d(n3, y);
                        }
                }
                t1 = clock();
                tmul_d += ((double) t1 - t0);
                */

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
        }

        printf("add:\t%lf\n", 1.0);
        printf("add_l:\t%lf\n", tadd_l / tadd);
        printf("mul:\t%lf\n", tmul / tadd);
        printf("mul_l:\t%lf\n", tmul_l / tadd);
        //printf("mul_d:\t%lf\n", tmul_d / tadd);
        printf("div:\t%lf\n", tdiv / tadd);
        printf("div_l:\t%lf\n", tdiv_l / tadd);

        return 0;
}
