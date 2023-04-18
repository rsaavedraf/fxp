/* SPDX-License-Identifier: MIT */
/*
 * dmul_tester.c
 * Tests the distributive multiplication functions.
 *
 * The distributive multiplication turned out to be quite
 * critical and needed in several places once trying to offer
 * ultimate precision for all fxp transcendental functions
 * (ln, exp, powxy etc) when using only ints. Also depending on
 * where the dmul functions are getting called from, different
 * combinations of argument types and result types are needed,
 * so separate implementations (so far 6) required to cover all
 * use cases:
 *
 * Using only ints:
 *     uint x uint --> uint
 *     uint x uint --> ulongy
 *     ulongy x uint --> ulongy
 *     ulongy x ulongy --> ulongy
 * Using longs:
 *     ulong x uint --> ulong
 *     ulong x ulong --> ulong
 *
 * This program tests all of these dmul implementations.
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
#include "fxp_l.h"
#include "fxp_aux.h"
#include "fxp_conv.h"
#include "print_as_bits.h"
#include "ulongy.h"

#define DASHES "========================\n"

#define SET_RAND_SEED 0
#define MAX_RAND_LOOPS 1

static void dump_uint(char * msg, unsigned int x)
{
        printf("%s: %10u (x%8X, \n\tb", msg, x, x);
        print_uint_as_bin(x);
        printf(")\n");
}

static void dump_ulong(char * msg, unsigned long x)
{
        printf("%s: %20lu  (x%16lX,\n\tb", msg, x, x);
        print_ulong_as_bin(x);
        printf(")\n");
}

static void dump_ulongy(char * msg, ulongy pong)
{
        unsigned int hi = ulongy_get_hi(pong);
        unsigned int lo = ulongy_get_lo(pong);
        printf("%s: %u|%u  (x%X|%X) \n\tb", msg, hi, lo, hi, lo);
        print_uint_as_bin(hi); printf("|");
        print_uint_as_bin(lo);
        printf(")\n");
}

static void test_uint(unsigned int expected, unsigned int actual)
{
        dump_uint("\texpected ", expected);
        dump_uint("\tactual   ", actual);
        fflush(stdout);
        assert(expected == actual);
        printf("\t(same)\n");
}

static void test_ulong(unsigned long expected, unsigned long actual)
{
        dump_ulong("\texpected ", expected);
        dump_ulong("\tactual   ", actual);
        fflush(stdout);
        assert(expected == actual);
        printf("\t(same)\n");
}

static void test_ulongy(ulongy expected, ulongy actual)
{
        dump_ulongy("\texpected ", expected);
        dump_ulongy("\tactual   ", actual);
        fflush(stdout);
        assert( ulongy_compare(expected, actual) == 0 );
        printf("\t(same)\n");
}

static void test_uints( unsigned int x, unsigned int y,
                        unsigned int expected)
{
        printf("\nTesting dmul_into_uint(uint, uint) --> uint:\n");
        dump_uint("\tx (uint) ", x);
        dump_uint("\ty (uint) ", y);
        //unsigned int actual = dmul_uints(x, y);
        unsigned int actual = dmul_into_uint(x, y);
        test_uint(expected, actual);
}

static void test_ulong_x_uint(  unsigned long x, unsigned int y,
                                unsigned long expected)
{
        printf("\nTesting dmul_ulong_x_uint(ulong, uint) --> ulong:\n");
        dump_ulong("\tx (ulong)", x);
        dump_uint( "\ty (uint) ", y);
        unsigned long actual = dmul_ulong_x_uint(x, y);
        test_ulong(expected, actual);
}

static void test_ulongs(unsigned long x, unsigned long y,
                        unsigned long expected)
{
        printf("\nTesting dmul_ulongs(ulong, ulong) --> ulong:\n");
        dump_ulong("\tx (ulong)", x);
        dump_ulong("\ty (ulong)", y);
        unsigned long actual = dmul_ulongs(x, y);
        test_ulong(expected, actual);
}

static void test_ulongy_from_dmul(  unsigned int x, unsigned int y,
                                    ulongy expected)
{
        printf("\nTesting ulongy_from_dmul(uint, uint) --> ulongy:\n");
        dump_uint("\tx (uint) ", x);
        dump_uint("\ty (uint) ", y);
        ulongy actual = ulongy_from_dmul(x, y);
        test_ulongy(expected, actual);
}

static void test_ulongy_x_uint( ulongy x, unsigned int y,
                                ulongy expected)
{
        printf("\nTesting dmul_ulongy_x_uint(ulongy, uint) --> ulongy:\n");
        dump_ulongy("\tx (ulongy)", x);
        dump_uint  ("\ty (uint)  ", y);
        ulongy actual = dmul_ulongy_x_uint(x, y);
        test_ulongy(expected, actual);
}

static void test_ulongys(ulongy x, ulongy y, ulongy expected)
{
        printf("\nTesting dmul_ulongys(ulongy, ulongy) --> ulongy:\n");
        dump_ulongy("\tx (ulongy)", x);
        dump_ulongy("\ty (ulongy)", y);
        ulongy actual = dmul_ulongys(x, y);
        test_ulongy(expected, actual);
}

static ulongy ulongy_from_ulong(unsigned long x)
{
        return ulongy_create(
                        (unsigned int) (x >> FXP_INT_BITS),
                        (unsigned int) (x & FXP_RINT_MASK));
}

int main(void)
{
        printf("\n%sTesting Distributive Multiplication functions (dmul) \n%s",
                DASHES, DASHES);
        // Multiplying pi by a number with a 1 in the
        // lowest bit of the highest word of the second operand.
        unsigned int x = FXP_PI_I32;
        unsigned int y = 1u << FXP_WORD_BITS;
        unsigned long xl = FXP_PI_I64;
        unsigned long yl = ((unsigned long) y) << FXP_INT_BITS;
        ulongy xlongy = ulongy_from_ulong(FXP_PI_I64);
        ulongy ylongy = ulongy_create(y, 0u);

        // The shifted and rounded expected results
        unsigned int rbit = (FXP_PI_I32 >> FXP_WORD_BITS_M1) & 1u;
        unsigned int expected_uint = (FXP_PI_I32 >> FXP_WORD_BITS) + rbit;
        rbit = (unsigned int) ((FXP_PI_I64 >> FXP_WORD_BITS_M1) & 1lu);
        unsigned long expected_ulong = (FXP_PI_I64 >> FXP_WORD_BITS) + rbit;
        ulongy expected_ulongy_1 = ulongy_from_ulong(expected_ulong);
        unsigned long xy = ((unsigned long) x) * y;
        ulongy expected_ulongy_2 = ulongy_from_ulong(xy);

        // Tests
        test_uints(x, y, expected_uint);
        test_ulongs(xl, yl, expected_ulong);
        test_ulong_x_uint(xl, y, expected_ulong);
        test_ulongy_x_uint(xlongy, y, expected_ulongy_1);
        test_ulongy_from_dmul(x, y, expected_ulongy_2);
        test_ulongys(xlongy, ylongy, expected_ulongy_1);

        // Test with random numbers
        if (SET_RAND_SEED) srand((unsigned int) time(0));
        for (int i = 0; i < MAX_RAND_LOOPS; i++) {
                printf("\nTest #%d with random numbers:\n", i + 1);
                unsigned int n1, n2, n3, n4, l3hi, l3lo, l4hi, l4lo;
                unsigned long l1, l2, l3, l4;
                ulongy p1, p2, p3, p4;
                n1 = rand();
                n2 = rand();
                n3 = rand();
                n4 = rand();
                printf("Testing mul of uints into ulong == dmul of uints into ulongy\n");
                test_ulongy(
                        ulongy_from_dmul(n1, n2),
                        ulongy_from_ulong(((unsigned long) n1) * n2));

                printf("Testing that dmul of ulongs == dmul of ulongys");
                l1 = (((unsigned long) n1) << FXP_INT_BITS) | n2;
                l2 = (((unsigned long) n3) << FXP_INT_BITS) | n4;
                l3 = dmul_ulongs(l1, l2);
                p1 = ulongy_create(n1, n2);
                p2 = ulongy_create(n3, n4);
                p3 = ulongy_from_ulong(l3);
                test_ulongys(p1, p2, p3);

                printf("Testing that dmul of ulong x uint == dmul of ulongy x uint");
                l4 = dmul_ulong_x_uint(l1, n3);
                p4 = ulongy_from_ulong(l4);
                test_ulongy_x_uint(p1, n3, p4);
        }

        printf("\nAll dmul tests passed successfully!!!\n\n");

        return 0;
}
