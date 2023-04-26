/* SPDX-License-Identifier: MIT */
/*
 * ulongy_tester.c
 *
 * By Raul Saavedra, Bonn Germany
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <float.h>
#include <time.h>
#include <math.h>

#include "fxp.h"
#include "fxp_extern.h"
//#include "fxp_l.h"
//#include "fxp_aux.h"
//#include "fxp_conv.h"
#include "print_as_bits.h"
#include "ulongy.h"

#define DASHES "========================\n"
#define RULER "----*----1----*----2----*----3-][-*----4----*----5----*----6----\n"


static void test_ulongy(unsigned long expected, ulongy actual)
{
        unsigned int ahi = ulongy_get_hi(actual);
        unsigned int alo = ulongy_get_lo(actual);
        unsigned long vactual = (((unsigned long) ahi) << FXP_INT_BITS) | alo;
        fflush( stdout );
        assert (expected == vactual);
}

static void print_ulong(char * msg, unsigned long x)
{
        printf("%s: %lX\n", msg, x);
        print_ulong_as_bin(x);
        printf("\n");
}

static void print_ulongy(char * msg, ulongy u)
{
        unsigned int uhi, ulo;
        uhi = ulongy_get_hi(u);
        ulo = ulongy_get_lo(u);
        printf("%s: {%X,%X}\n", msg, uhi, ulo);
        print_uint_as_bin(uhi);
        print_uint_as_bin(ulo);
        printf("\n");
}


int main(void)
{
        printf("\n%sTester of ulongy:\n%s", DASHES, DASHES);

        unsigned int x = FXP_POS_INF - 1, y = 10;
        unsigned long tgt = (((unsigned long) x) << FXP_INT_BITS) | y;
        printf(RULER);
        printf("x: %X,  y: %X\n", x, y);
        print_ulong("\nUnsigned long from x and y", tgt);
        ulongy u = ulongy_create(x, y);
        print_ulongy("Ulongy from x and y", u);
        test_ulongy(tgt, u);

        unsigned long sum = tgt + tgt;
        print_ulong("\nSum", sum);
        ulongy usum = ulongy_add(u, u);
        print_ulongy("Ulongy sum", usum);
        test_ulongy(sum, usum);

        unsigned long rs = sum >> 16;
        print_ulong("\nRshifted", rs);
        ulongy urs = ulongy_rshift(usum, 16);
        print_ulongy("Ulongy rshifted", urs);
        test_ulongy(rs, urs);

        unsigned long zero = 0lu;
        unsigned long one = 1lu;
        unsigned long two = 2lu;
        ulongy uzero = ulongy_create(0u, 0u);
        ulongy uone = ulongy_create(0u, 1u);
        ulongy utwo = ulongy_create(0u, 2u);
        unsigned long d0 = zero - two;
        ulongy ud0 = ulongy_sub(uzero, utwo);
        print_ulong("\nDifference 0", d0);
        print_ulongy("Ulongy difference 0", ud0);
        test_ulongy(d0, ud0);

        unsigned long d1 = sum - rs;
        print_ulong("\nDifference 1", d1);
        ulongy ud1 = ulongy_sub(usum, urs);
        print_ulongy("Ulongy difference 1", ud1);
        test_ulongy(d1, ud1);

        unsigned long d2 = rs - sum;
        print_ulong("\nDifference 2", d2);
        ulongy ud2 = ulongy_sub(urs, usum);
        print_ulongy("Ulongy difference 2", ud2);
        test_ulongy(d2, ud2);

        printf("\n\tAll tests passed!\n\n");
}
