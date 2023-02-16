/*
 * test_quad.c
 *
 * Trying to get Euler's number 'e' from my Intel i7-6700K-based PC
 * with the same precision that I can already get (with absolutely
 * no extra effort) from a tiny Raspberry Pi 4b just using long double.
 *
 * Based on: https://gcc.gnu.org/onlinedocs/libquadmath/quadmath_005fsnprintf.html#quadmath_005fsnprintf
 *
 * By Raul Saavedra, Bonn Germany
 * 2023-02-16
 *
 */

#include <quadmath.h>
#include <stdlib.h>
#include <stdio.h>

#define STR_E_DEC  "2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274"

int main ()
{
        __float128 r;
        int prec = 35;
        int width = 50;
        char buf[128];

        printf("Testing usage of quadmath.h on an Intel-based PC.\n");
        printf("33 decimal digits of e should now be correct here:\n\n");
        printf("            ----*----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8----*----9----*----0\n");
        printf("'True' e: %s\n", STR_E_DEC);

        r = M_Eq;
        int n = quadmath_snprintf (buf, sizeof buf, "%+-#*.50Qe", width, r);
        if ((size_t) n < sizeof buf)
        printf ("M_Eq    :%s\n\n", buf);

        return 0;
}
