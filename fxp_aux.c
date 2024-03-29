/* SPDX-License-Identifier: MIT */
/*
 * fxp_aux.c
 *
 * By Raul Saavedra, Bonn Germany
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <time.h>
#include <math.h>
#include "fxp_extern.h"
#include "fxp_aux.h"
#include "print_as_bits.h"
//#include <assert.h>

void print_fxp_as_bin(int n)
{
        int an;
        int signspace = 0;
        if (n < 0) {
                an = -n;
                printf("-b");
        } else {
                an = n;
                printf(" b");
        }
        int nbn = fxp_nbits(an);
        int i = FXP_INT_BITS;
        while (i > 0) {
                //if (i == FXP_WORD_BITS) printf("_");
                if (i == FXP_frac_bits) printf(".");
                if (i <= nbn) {
                        int bit = (an >> (i - 1)) & 1;
                        printf("%d", bit);
                } else {
                        printf("0");
                }
                i--;
        }
}

void print_fxp(int x)
{
        if (x == FXP_POS_INF || x == FXP_NEG_INF \
                || x == FXP_UNDEF) {
                printf(x == FXP_UNDEF? "UNDEF":
                          (x == FXP_POS_INF? "+INF": "-INF"));
                return;
        }
        long double n = fxp2ld(x);
        if (x < 0) {
                int px = -x;
                printf("%17.10Le (fxp: %11d, -x%08X, ",
                        n, x, px);
        } else {
                printf("%17.10Le (fxp: %11d,  x%08X, ",
                        n, x, x);
        }
        print_fxp_as_bin(x);
        printf(")");
}

void print_fxp_div(int startmask, int nmaskbits, int n, int frac_bits)
{
        printf("   ");
        int i = startmask - 1;
        while (i > 0) {
                printf(" ");
                i--;
        }
        i = nmaskbits;
        while (i > 0) {
                printf(".");
                i--;
        }
        printf("v\n x:");
        int nbn = fxp_nbits(n);
        print_int_as_bin(n, 0);
        printf(" \t\t(%d, x%x, %d bits)\n   ", n, n, fxp_nbits(n));
        i = nbn;
        while (i > 0) {
                printf(" ");
                i--;
        }
        i = frac_bits;
        while (i > 0) {
                printf("0");
                i--;
        }
}

void trace_fxp_div( char * msg,
            int iteration, int frac_bits, int bindex, \
            int dividend, int divisor, \
            int quotient, int newqbit, \
            int mc, \
            int difference)
{
        int bmc = fxp_nbits(mc);
        int nbdividend = fxp_nbits(dividend);
        if (iteration == 0) {
                printf("\n%s frac_bits:%d  bindex:%2d\n", \
                        msg, frac_bits, bindex);
                print_fxp_div(bindex - bmc, bmc, dividend, FXP_frac_bits);
                printf("    \ty: ");
                print_int_as_bin(divisor, 0);
                printf(" (%d, x%x, %d bits)\n", \
                        divisor, (unsigned int) divisor, fxp_nbits(divisor));
        }
        printf(" m:");
        print_int_as_bin(mc, bindex);
        printf("  (x%x, %d bits) \tq:",
                (unsigned int) mc, bmc);
                print_fxp(quotient);
        printf("\n");

        int substraend = newqbit * divisor;
        printf("-s:");
        print_int_as_bin(substraend, bindex);
        printf("  (x%x, %d bits)\n",
                (unsigned int) substraend, fxp_nbits(substraend));
        int nxpdbit, shift;
        if (bindex < nbdividend) {
                shift = nbdividend - bindex - 1;
                nxpdbit = (dividend & (1 << shift)) >> shift;
        } else {
                nxpdbit = 0;
        }
        printf(" r:");
        print_int_as_bin(difference, bindex);
        printf("  (x%x, %d bits)   (bidx:%d, pd-bit:%d -> next m:x%x)\n",
                (unsigned int) difference, fxp_nbits(difference),
                bindex, nxpdbit, ((difference << 1) | nxpdbit));
}

void print_sys_info()
{
        fflush( stdout );
        printf("\nSystem details:\n");
        int r1 = system("hostnamectl | grep -e 'Operating System' -e 'Architecture'");
        int r2 = system("cat /proc/cpuinfo | grep -e 'CPU' -e 'model name' -e 'Model' | sort -r | head -1");

        printf("\nType sizes on this system (some might depend on compiler options):\n");
        printf("char        has a size of %zu bytes.\n", sizeof(char));
        printf("short       has a size of %zu bytes.\n", sizeof(short));
        printf("int         has a size of %zu bytes.\n", sizeof(int));
        printf("long        has a size of %zu bytes.\n", sizeof(long));
        printf("long long   has a size of %zu bytes.\n", sizeof(long long));
        printf("float       has a size of %zu bytes.\n", sizeof(float));
        printf("double      has a size of %zu bytes.\n", sizeof(double));
        printf("long double has a size of %zu bytes.\n", sizeof(long double));
}
