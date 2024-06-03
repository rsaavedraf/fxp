/* SPDX-License-Identifier: MIT */
/*
 * sqrt.c
 *
 * Initial implementations of square root using CORDIC and long doubles,
 * as explained in the Mathworks pages:
 * https://www.mathworks.com/help/fixedpoint/ug/compute-square-root-using-cordic.html#ComputeSquareRootUsingCORDICExample-4
 *
 * Written to test the implemented algorithm before tailoring it for fxp's
 *
 * By Raul Saavedra, Bonn Germany
 */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include "fxp_extern.h"
#include "print_as_bits.h"

//#define VERBOSE 1
#define DASHES "=========================================\n"

#define SET_RAND_SEED 0

const long double FXP_ZERO_LD = 1.0E-124L;
//const int LONG_BITS = sizeof(long) * 8;

// This scaler is the inverse of product of sqrt(1 - 2^(-2i)), but with element i in the sequence
// repeated for every 3k + 1 steps, k starting with 1. So operands with i = 4, 13, 40, 121...
// get repeated in the product to achieve result convergence.
// Here using the scaler for i=20
const long double CORDIC_SQRT_SCALER = 1.20749706776288909351346436569307772928L;
/* Scalers for different i values using binary128 precision:
i:10, prod: 0.82815949876289316388069518519831971464, 1/prod: 1.20749686684002614056886173133121622865
prod ulong value: 6A0121658778F32B  1/prod ulong value: 4D47A0F0466C9C9E
i:20, prod: 0.82815936096034116150136818623292654789, 1/prod: 1.20749706776288909351346436569307772928
prod ulong value: 6A01203D99B7E461  1/prod ulong value: 4D47A1C803AE2785
i:30, prod: 0.82815936096021562719591762738682837384, 1/prod: 1.20749706776307212870316438601668749092
prod ulong value: 6A01203D99A63987  1/prod ulong value: 4D47A1C803BB08CA
*/

/*
 * Calculates the square root using Cordic
 */
long double my_cordic_sqrt_kernel(long double xin, unsigned int loops)
{
        printf("\tkernel calculation for sqrt(%.15Lf)\n", xin);
        int k = 4; // Used for the repeated (3*k + 1) iteration steps
        long double x = xin + 0.25L;
        long double y = xin - 0.25L;
        long double npw2 = 0.5L, xtmp, ytmp;
        for (int idx = 1; idx <= loops; idx++) {
                xtmp = x * npw2;
                ytmp = y * npw2;
                if (y < 0.0L) {
                    x += ytmp;
                    y += xtmp;
                } else {
                    x -= ytmp;
                    y -= xtmp;
                }
                if (idx == k) {
                        xtmp = x * npw2;
                        ytmp = y * npw2;
                        if (y < 0.0L) {
                            x += ytmp;
                            y += xtmp;
                        } else {
                            x -= ytmp;
                            y -= xtmp;
                        }
                        k = 3*k + 1;
                }
                npw2 = npw2 / 2.0L;
        }
        //printf("%.5Lf\n", x);
        return x;
}

long double my_sqrt(long double v, unsigned int loops)
{
        if (v < 0.0L) return FXP_UNDEF_LD;
        if (v == 0.0L) return 0.0L;
        // Find an even integer c, so that:
        // v = u * 2^c (with 0.5 <= u < 2)
        // Then we calculate sqrt(v) = sqrt(u) * 2^(c/2)
        long double u = v, res;
        int c = 0;
        if (u >= 2.0L) {
                while (u >= 2.0L) {
                        c++;
                        u /= 2;
                }
                if (c % 2) {
                        c++;
                        u /= 2;
                }
                printf("\tu>=2, c:%d\n", c);
                res = my_cordic_sqrt_kernel(u, loops);
                if (c == 0) {
                        return res * CORDIC_SQRT_SCALER;
                } else {
                        printf("\tkernel result: %.15Lf\n", res);
                        long double scaled = res * CORDIC_SQRT_SCALER;
                        printf("\tscaled       : %.15Lf\n", scaled);
                        long double shifted = scaled * (1ul << (c >> 1));
                        printf("\tshifted      : %.15Lf\n", shifted);
                        return shifted;
                }
        } else {
                while (u < 0.5L) {
                        c++;
                        u *= 2;
                }
                if (c % 2) {
                        c++;
                        u *= 2;
                }
                //printf("\tu<2, c:-%d\n", c);
                res = my_cordic_sqrt_kernel(u, loops);
                return res * CORDIC_SQRT_SCALER / (1ul << (c >> 1));
        }
}


int main(void)
{
        if (SET_RAND_SEED) srand((unsigned int) time(0));

        printf("\n%sCalculating constants to use as gains for CORDIC square root:\n%s", DASHES, DASHES);
        long double gain = 1;
        int k = 4;
        long double product = 1.0L;
        for (int i=1; i<=30; i++) {
                // compute sqrt( 1 - 2^(-2i) )
                long double s = powl(2.0L, 2*i);
                long double p = sqrtl( 1.0L - 1.0L/s );
                product *= p;
                if (i == k) {
                        product *= p;
                        k = 3*k + 1;
                }
                if ((i == 1) || ((i % 10) == 0) || (i == 64)) {
                        long double invproduct = 1.0L / product;
                        printf("i:%2d, prod:%41.38Lf, 1/prod:%41.38Lf\n", i, product, invproduct);
                        //unsigned long ulvalue = get_ulong_bits_from_ldouble(product);
                        unsigned long ulvalue = inspect_long_double_aux(product, 0);
                        unsigned long ulinvvalue = get_ulong_bits_from_ldouble(invproduct);
                        printf("prod ulong value: %lX  ", ulvalue);
                        printf("1/prod ulong value: %lX\n", ulinvvalue);
                }
        }
        printf("\n\n");
        printf("\n%sTesting square root calculation using Cordic\n%s", DASHES, DASHES);

        // arguments of interest to check
        //long double aoi[] = {0.000000001L, 0.25L, 0.5L, 0.77L, 1.0L, 1.3L, 1.5L, 2.0L, 3.0L, 4.0L, 9.0L, 16.0L, 57.0L, FXP_MAX+0.5L};
        long double aoi[] = {9.0L};
        int naoi = (int) (sizeof(aoi) / sizeof(aoi[0]));

        //int n, loops = 10;
        for (int loops = 5; loops <= 20; loops += 5) {
                long double avg = 0.0L;
                printf("Calculations using %d loops:\n", loops);
                for (int i = 0; i < naoi; i++) {
                        long double x = aoi[i];
                        long double expected = sqrtl(x);
                        long double y = my_sqrt(x, loops);
                        long double delta = (y >= expected? y - expected: expected - y);
                        avg += delta;
                        printf("\tExpected   sqrt(%9.6Lf) = %41.38Lf\n", x, expected);
                        printf("\tCalculated sqrt(%9.6Lf) = %41.38Lf (delta: %5.3LE)\n\n", x, y, delta);
                }
                avg = avg/naoi;
                printf("Average delta for %d loops: %5.3LEf\n\n", loops, avg);
        }
}
