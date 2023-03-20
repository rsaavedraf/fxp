/* SPDX-License-Identifier: MIT */
/*
 * fxp_aux.c
 *
 * By Raul Saavedra, Bonn Germany
 */

#include <stdio.h>
#include <stdlib.h>
#include "print_as_bits.h"

const int INT_BITS = ((int) sizeof(int)) * 8;
const int LONG_BITS = ((int) sizeof(long)) * 8;

void print_int_as_bin(int n, int width)
{
        int an, neg_sign;
        if (n < 0) {
                an = -n;
                neg_sign = 1;
        } else {
                an = n;
                neg_sign = 0;
        }
        int nbn = INT_BITS - __builtin_clz(an);
        int margin = (nbn == 0? 1: nbn) - neg_sign;
        while (width > margin) {
                printf(" ");
                width--;
        }
        if (neg_sign) printf("-");
        nbn = (nbn == 0)? 1: nbn;
        int i = nbn;
        while (i > 0) {
                int bit = (an >> (i - 1)) & 1;
                printf("%d", bit);
                i--;
        }
}

void print_long_as_bin(long n)
{
        int an;
        if (n < 0) {
                an = -n;
                printf("-");
        } else {
                an = n;
        }
        int i = LONG_BITS;
        while (i > 0) {
                int bit = (an >> (i - 1)) & 1;
                printf("%d", bit);
                i--;
        }
}

void print_uint_as_bin(unsigned int n)
{
        int i = INT_BITS;
        while (i > 0) {
                int bit = (n >> (i - 1)) & 1;
                printf("%d", bit);
                i--;
        }

}

void print_ulong_as_bin(unsigned long n)
{
        int i = LONG_BITS;
        while (i > 0) {
                int bit = (n >> (i - 1)) & 1;
                printf("%d", bit);
                i--;
        }
}
