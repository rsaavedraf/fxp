/* SPDX-License-Identifier: MIT */
/*
 * fxp_tconst.c
 *
 * Auxiliary program to generate fxp versions of some
 * important/transcendental constants. It also shows
 * the actual bit-depth precision of long doubles
 * in the current system, comparing expl() vs. the
 * fractional part of e in binary as reference
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
#include "fxp_aux.h"

#define DASHES "========================\n"

#define LONGLONG_BITS ((int) (sizeof(long long) * 8))

#define STR_TEST1 "0000000000000000011111100000000000"
#define STR_TEST2 "0000000000000000100100000000000000"
#define STR_TEST3 "0000000000000001111100000000000000"

// Euler's constant (e) in decimal and binary
//                    ----*----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8----*----9----*----0
#define STR_E_DEC  "2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274"
#define STR_E_BIN "10.101101111110000101010001011000101000101011101101001010100110101010111111011100010101100010000000100111001111010011110011110001"
// For more details, or even more precise versions of e:
//  https://apod.nasa.gov/htmltest/gifcity/e.2mil
//  https://www.math.utah.edu/~pa/math/e.html
//  https://rosettacode.org/wiki/Calculating_the_value_of_e#C
//  https://www.exploringbinary.com/pi-and-e-in-binary/

// Pi in decimal and binary
//                     ----*----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8----*----9----*----0
#define STR_PI_DEC  "3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679"
#define STR_PI_BIN "11.0010010000111111011010101000100010000101101000110000100011010011000100110001100110001010001011100000"
// For more details, or even more precise versions of Pi:
//  https://oeis.org/A004601/constant
//  https://oeis.org/A000796/constant
//  http://www.befria.nu/elias/pi/binpi.html
//  https://www.exploringbinary.com/pi-and-e-in-binary/

// ln(2)
#define STR_LN_2_DEC ".693147180559945309417232121458176568075500134360255254120680009493393621969694715605863326996418687"
// https://oeis.org/A002162

// log10(2)
#define STR_LG10_2_DEC ".301029995663981195213738894724493026768189881462108541310427461127108189274424509486927252118186172040684"
//  https://oeis.org/A007524

// log2(e)
#define STR_LG2_E_DEC   "1.4426950408889634073599246810018921374266459541529859341354494069311092191811850798855266228935063444"

// log2(10)
#define STR_LG2_10_DEC "3.32192809488736234787031942948939017586483139302458061205475639581593477660862521585013974335937015"

// Sqrt(2)
//#define STR_SQRT2_BIN "1.01101010000010011110011001100111111100111011110011001001000010001011001011111011000100110110011011"
//#define STR_SQRT2_DEC "1.41421356237309504880168872420969807856967187537694807317667973799073247846210703885038753432764157"
//  https://community.wolfram.com/groups/-/m/t/1063480
//  https://oeis.org/A002193
//  https://oeis.org/A004539


//const long double FXP_ZERO_LD = 1.0E-124;

/*
 * Return binary expansion (up to nbits) of a string representation
 * of a binary number, with last bit either verbatim or rounded
 */
unsigned long long bex_from_bin(char * pbinnum, int wbits, int fbits, int rounded)
{
        char * p = pbinnum;
        if ((wbits <= 0) || (wbits > 32)) wbits = 32;
        if ((fbits <= 0) || (fbits + wbits > 64)) fbits = 64 - wbits;
        int roomleft = wbits + fbits;
        unsigned long long bnum = 0;
        while ((roomleft > 0) && (*p != '\0')) {
                char cbit = *p;
                if (cbit != '.') {
                        bnum = (bnum << 1) | ((cbit == '1')? 1: 0);
                        roomleft--;
                }
                p++;
        }
        if (roomleft > 0) {
                bnum = bnum << roomleft;
                // Source string with too few bits to need any rounding,
                // so return already
                return bnum;
        }
        if (rounded) {
                //printf("  (Next bit:");
                int nextbit = 0;
                char cbit = *p;
                //printf("%c", cbit);
                int rbit = ((cbit  == '1')? 1: 0);
                //printf(") -> adding %d", rbit);
                bnum += rbit;
        }
        return bnum;
}

/*
 * Return binary expansion (up to wbits.fbits) of a string representation
 * of a decimal number, with last bit either verbatim or rounded
 */
unsigned long long bex_from_dec(char * pdecnum, int wbits, int fbits, int rounded)
{
        char * p = pdecnum;
        if ((wbits < 0) || (wbits > 31)) wbits = 31;
        if ((fbits < 0) || (fbits + wbits > 64)) fbits = 63 - wbits;
        // Process whole part first
        int nwd = 0;
        int ptens = 1;
        while (*p != '.') {
                p++;
                ptens *= 10;
                nwd++;
        }
        int wnum = 0;
        p = pdecnum;
        while (*p != '.') {
                char cdigit = *p;
                ptens /= 10;
                wnum += ((int) cdigit - '0') * ptens;
                p++;
        }
        //printf("Whole part is d:%d (%d dec digits)\n", wnum, nwd);
        long long binnum = 0;
        int nb = 0;
        while (wnum > 0) {
                int r = wnum % 2;
                //printf("Appending a %d to the whole part\n", r);
                binnum = (r << nb) | binnum;
                nb++; // here count the whole bits
                wnum /= 2;
        }
        //printf("binnum: (x:%llx. , %d bits)\n", binnum, nb);
        // Process fractional part
        p++; // skip the '.'
        int ndd = nwd;
        long double ftens = 1.0/10.0;
        long double frac = 0;
        int ndig = 0;
        while (*p != '\0') { //&& (ftens > ZERO)) {
                char cdigit = *p;
                frac += ((int) cdigit - '0') * ftens;
                ftens /= 10;
                ndig++;
                ndd++;
                p++;
        }
        //printf("From %d decimal frac digits, fraction is: %1.20Lf \n", ndig, frac);
        int roomleft = fbits;
        long double m;
        while (roomleft > 0) {
                m = frac * 2.0;
                if (m >= 1.0) {
                        binnum = ((binnum << 1) | 1);
                        frac = m - 1.0;
                        //printf("1");
                } else {
                        binnum = (binnum << 1);
                        frac = m;
                        //printf("0");
                }
                roomleft--;
        }
        if (rounded) {
                m = frac * 2.0;
                int rbit = ((m >= 1.0)? 1: 0);
                binnum += rbit;
                //printf(" (+ round bit %d)", rbit);
        }
        //printf("\n");
        return binnum;
}

/*
 * Return a long double from a str representation of a decimal number
 */
long double Lf_from_dec(char * pdecnum)
{
        char * p = pdecnum;
        // Process whole part first
        int nwd = 0;
        int ptens = 1;
        while (*p != '.') {
                p++;
                ptens *= 10;
                nwd++;
        }
        long double num = 0.0;
        p = pdecnum;
        while (*p != '.') {
                char cdigit = *p;
                ptens /= 10;
                num += ((int) cdigit - '0') * ptens;
                p++;
        }
        // Process fractional part
        p++; // skip the '.'
        long double tenth = 1.0 / 10.0;
        long double ftens = tenth;
        while (*p != '\0') {
                char cdigit = *p;
                num += ((int) cdigit - '0') * ftens;
                ftens *= tenth;
                p++;
        }
        return num;
}

long double Lf_from_bin(char * pbinnum)
{
        char * p = pbinnum;
        // Process whole part first
        int nwd = 0;
        int ptwo = 1;
        while (*p != '.') {
                p++;
                ptwo *= 2;
                nwd++;
        }
        long double num = 0.0;
        p = pbinnum;
        while (*p != '.') {
                char cdigit = *p;
                ptwo /= 2;
                num += ((int) cdigit - '0') * ptwo;
                p++;
        }
        // Process fractional part
        p++; // skip the '.'
        long double half = 0.5;
        long double frac = half;
        while (*p != '\0') {
                char cdigit = *p;
                num += ((int) cdigit - '0') * frac;
                frac *= half;
                p++;
        }
        return num;
}

int main(void)
{
        printf("\n%sfxp_tconst.c\n%s", DASHES, DASHES);
        print_sys_info();

        //printf("\ntest1 %llx\n", bex_from_bin(STR_TEST1, 16, 1));
        //printf("\ntest2 %llx\n", bex_from_bin(STR_TEST2, 16, 1));
        //printf("\ntest3 %llx\n", bex_from_bin(STR_TEST3, 16, 1));

        printf("\nPrecisions of e on this system:\n");
        printf("Frac digits   :   ----*----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8----*----9----*----0\n");
        const long double E_EXP = (long double) exp(1);
        const long double MY_E_DEC = Lf_from_dec(STR_E_DEC);
        const long double MY_E_BIN = Lf_from_bin(STR_E_BIN);
        const long double E_EXPL = (long double) expl(1);
        printf("e from M_E    : %.64LE (- expl(): %1.4LE)\n",
                (long double) M_E, (M_E - E_EXPL));
        printf("e from exp()  : %.64LE (- expl(): %1.4LE)\n",
                E_EXP, (E_EXP - E_EXPL));
        printf("My e (dec str): %.64LE (- expl(): %1.4LE)\n",
                MY_E_DEC, (MY_E_DEC - E_EXPL));
        printf("My e (bin str): %.64LE (- expl(): %1.4LE)\n",
                MY_E_BIN, (MY_E_BIN - E_EXPL));
        printf("e from expl() : %.64LE (- expl(): %1.4LE)\n",
                E_EXPL, (E_EXPL - E_EXPL));
        printf("'True' e (dec): %s\n", STR_E_DEC);
        printf("'True' e (bin):%s\n", STR_E_BIN);
        printf("bin expl()    :");
        int intbits = sizeof(unsigned int) * 8;
        unsigned int first_mask = 1u << (intbits - 1);
        long double eshifted = E_EXPL * powl(2.0, intbits -2);
        int e_ref_index = 0;
        int binpoint_seen = 0;
        int match_count = 0;
        int checking = 1;
        do {
                //printf("eshifted is: %Lf\n", eshifted);
                unsigned int e_chunk = truncl(eshifted);
                //printf("e, %d-bits chunk: %x \n", intbits, e_chunk);
                // Compare bits in e_chunk vs. the reference
                unsigned int mask = first_mask;
                while (mask > 0) {
                    unsigned int ebit = e_chunk & mask;
                    char c_refbit = STR_E_BIN[e_ref_index];
                    if ((!binpoint_seen) && (c_refbit == '.')) {
                        c_refbit = STR_E_BIN[++e_ref_index];
                        binpoint_seen = 1;
                        printf(".");
                    }
                    unsigned int refbit = (c_refbit == '1')? 1u: 0u;
                    if ((refbit && (ebit == 0)) || ((!refbit) && (ebit > 0))) {
                        printf("(%d)\n", (ebit == 0? 0: 1));
                        printf("Exact match up to %d frac bits (%d-bit total match)\n", \
                                match_count - 2, match_count);
                        long double ddigits = (match_count - 2) / (logl(10)/logl(2));
                        printf("Corresponds to exact match up to %d decimal frac digits.\n\n", (int) truncl(ddigits));
                        checking = 0;
                        break;
                    }
                    printf((ebit > 0)? "1": "0");
                    match_count++;
                    mask >>= 1;
                    e_ref_index++;
                }
                eshifted = eshifted - ((long double) e_chunk);
                eshifted = eshifted * powl(2, intbits);
        } while (checking);

        //printf("bin e (expl()): %lx\n\n", bin_e);

        int frbits = 30;
        printf("e as binary: %s\n", STR_E_BIN);
        printf("e as fxp (%d frac bits) %llX\n", \
                    frbits, bex_from_bin(STR_E_BIN, 2, frbits, 1));
        printf("e as fxp (%d frac bits) %llX\n\n", \
                    32 + frbits, bex_from_bin(STR_E_BIN, 2, 32 + frbits, 1));

        printf("pi as binary: %s\n", STR_PI_BIN);
        printf("pi as fxp (%d frac bits) %llX\n", \
                    frbits, bex_from_bin(STR_PI_BIN, 2, frbits, 1));
        printf("pi as fxp (%d frac bits) %llX\n\n", \
                    32 + frbits, bex_from_bin(STR_PI_BIN, 2, 32 + frbits, 1));

        printf("lg2(10) as decimal: %s\n", STR_LG2_10_DEC);
        printf("lg2(10)  as fxp (%d frac bits) %llX\n", \
                    frbits, bex_from_dec(STR_LG2_10_DEC, 0, frbits, 1));
        printf("lg2(10)  as fxp (%d frac bits) %llX\n\n", \
                    frbits, bex_from_dec(STR_LG2_10_DEC, 0, 32 + frbits, 1));

        frbits = 31;
        printf("lg2(e)  as decimal: %s\n", STR_LG2_E_DEC);
        printf("lg2(e)  as fxp (%d frac bits) %llX\n", \
                    frbits, bex_from_dec(STR_LG2_E_DEC, 0, frbits, 1));
        printf("lg2(e)  as fxp (%d frac bits) %llX\n\n", \
                    frbits, bex_from_dec(STR_LG2_E_DEC, 0, 32 + frbits, 1));

        frbits = 32;
        printf("ln(2) as decimal: %s\n", STR_LN_2_DEC);
        printf("ln(2) as fxp (%d frac bits) %llX\n", \
                    frbits, bex_from_dec(STR_LN_2_DEC, 0, frbits, 1));
        printf("ln(2) as fxp (%d frac bits) %llX\n\n", \
                    32 + frbits, bex_from_dec(STR_LN_2_DEC, 0, 32 + frbits, 1));

        printf("lg10(2) as decimal: %s\n", STR_LG10_2_DEC);
        printf("lg10(2) as fxp (%d frac bits) %llX\n", \
                    frbits, bex_from_dec(STR_LG10_2_DEC, 0, frbits, 1));
        printf("lg10(2) as fxp (%d frac bits) %llX\n\n", \
                    32 + frbits, bex_from_dec(STR_LG10_2_DEC, 0, 32 + frbits, 1));

        return 0;
}
