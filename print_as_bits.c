/* SPDX-License-Identifier: MIT */
/*
 * print_as_bits.c
 *
 * Utility functions to print ints, longs, and long doubles as binaries
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
        unsigned long an;
        if (n < 0) {
                an = -n;
                printf("-");
        } else {
                an = n;
        }
        int i = LONG_BITS;
        while (i > 0) {
                int bit = (int) ((an >> (i - 1)) & 1ul);
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
                int bit = (n >> (i - 1)) & 1ul;
                printf("%d", bit);
                i--;
        }
}

void shift_ulongs(unsigned long *hi, unsigned long *lo, int shift)
{
        unsigned long auxhi, auxlo;
        if (shift > 0) {
                if (shift >= 128) {
                        *hi = 0ul;
                        *lo = 0ul;
                        return;
                }
                if (shift >= 64) {
                        *hi = 0ul;
                        *lo = *hi >> (shift - 64);
                } else {
                        auxhi = *hi;
                        *hi >>= shift;
                        *lo = (auxhi << (64 - shift)) | (*lo >> shift);
                }
        } else if (shift < 0) {
                if (shift <= -128) {
                        *hi = 0ul;
                        *lo = 0ul;
                        return;
                }
                if (shift <= -64) {
                        *lo = 0ul;
                        *hi = *lo << (-shift - 64);
                } else {
                        auxlo = *lo;
                        *lo <<= -shift;
                        *hi = (*hi << -shift) | (auxlo >> (64 + shift));
                }
        }
}

/*
 * Print out a long double (quadruple precision) following IEEE-754 standard
 * 1 sign bit
 * 15 exponent bits
 * 112 bits explicitely stored, but 113 bits of significant precision
 * True exponent is the stored one minus the offset of 16383
 * Exponents 0000 and 7FFF are intrepreted specially:
 * 0000        -> subnormal numbers:  (-1)^signbit *        2^(-16382)        * 0.significandbits
 * 0001 - 7FFE -> normalized value:   (-1)^signbit * 2^(exponentbits - 16383) * 1.significandbits
 * 7FFF        -> NaN
 * Pi: 4000 921f b544 42d1 8469 898c c517 01b8
 * -2: c000 0000 0000 0000 0000 0000 0000 0000
 *  1: 3fff 0000 0000 0000 0000 0000 0000 0000
 *  0: 0000 0000 0000 0000 0000 0000 0000 0000
 * -0: 8000 0000 0000 0000 0000 0000 0000 0000
 * For more details:
 * https://en.wikipedia.org/wiki/Quadruple-precision_floating-point_format
 *
 * Now also supporting long double in x86 Extended Precision format, which is
 * not exactly the same as IEEE-754. For more details:
 * https://en.wikipedia.org/wiki/Extended_precision
 */
void inspect_long_double(long double x, int VERBOSE)
{
        if (VERBOSE) printf("%40.38Lf\n", x);
        void *px = &x;
        unsigned char *pc = px;
        // Assuming quadruple-precision/binary128 long doubles (16 bytes)
        unsigned char bytes[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        int ai = 0;
        int sum = 0;
        int exbias = -16383;
        int implicit_bit = 0;
        int shift_base = -14;
        #ifdef __ARM_ARCH
                for (int i = 15; i >= 0; i--) {
                        unsigned char c = *(pc + i);
                        bytes[ai] = c;
                        ai++;
                        sum += c;
                }
        #elif __x86_64__
                // x86 does not quite use IEEE-754 encoding, but its own 80-bit
                // floating point format for long doubles
                for (int i = 9; i >= 0; i--) {
                        unsigned char c = *(pc + i);
                        bytes[ai] = c;
                        ai++;
                        sum += c;
                }
                shift_base = -15;
        #endif
        if (VERBOSE) {
                printf("\tAs stored : ");
                for (int i = 0; i < 16; i++) {
                        printf("%02x", bytes[i]);
                        if ((i > 0) && (i % 2)) printf(" ");
                }
        }
        int sign = ((bytes[0] & 0x80) == 0x80);
        int ebits = ((((int) (bytes[0] & 0x7F)) << 8) | bytes[1]);
        int exponent = ebits + exbias;
        if (VERBOSE) {
                printf("\n\tsign      : %1d", sign);
                printf("\n\texponent  : %d", exponent);
                printf("\n\tmagnitude : ");
        }
        unsigned long hi = 0ul;
        unsigned long lo = 0ul;
        if (ebits == 0x7FFF) {
                // Special cases
                if (VERBOSE) {
                        if (sum == 0) {
                                printf("<inf>");
                        } else {
                                printf("<NaN>");
                        }
                }
        } else {
                #ifdef __ARM_ARCH
                        // Include the implicit bit as in IEEE-754
                        if (ebits == 0x0) {
                                // Special case for subnormals
                                if (VERBOSE) printf("(0). ");
                        } else {
                                // normalized values
                                if (VERBOSE) printf("(1). ");
                                hi = 1ul;
                        }
                #elif __x86_64__
                        // No implicit/hidden bit in x86's Extended Precision format
                #endif
                int hiroomleft = 48; // 2 bytes in the "hi" ulong used already
                // fraction bits
                for (int i = 2; i < 16; i++) {
                        unsigned char c = bytes[i];
                        if (VERBOSE) {
                                printf("%02x", c);
                                if ((i > 0) && (i % 2)) printf(" ");
                        }
                        if (hiroomleft >= 8) {
                                hi = (hi << 8) | c;
                                hiroomleft -= 8;
                        } else {
                                lo = (lo << 8) | c;
                        }
                }
        }
        if (VERBOSE) {
                printf("\n\tAs ulongs : 0x000%lx 0x%016lx", hi, lo);
                printf("\n\tAs bits   : ");
                print_ulong_as_bin(hi); printf(" "); print_ulong_as_bin(lo);
        }
        // L-shift the bits
        shift_ulongs(&hi, &lo, shift_base + ((exponent < 0)? -exponent - 1: 0));
        if (VERBOSE) {
                printf("\n\tL-shifted : ");
                print_ulong_as_bin(hi); printf(" "); print_ulong_as_bin(lo);
                printf("\n\t            ");
        }
        #ifdef __ARM_ARCH
                // Round last bit in hi ulong
                if ((lo >> 63) & 1ul) hi++;
        #elif __x86_64__
                // Only bits 0 to 62 hold the fractional part in the x86 Extended
                // Precision format, so ignore bit 63.
        #endif
        printf("0x%016lX", hi);
        if (VERBOSE) {
                printf("\n");
        } else {
                printf(", ");
        }
}
