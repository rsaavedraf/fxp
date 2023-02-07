/*
 * fxp_constants.c
 *
 * By Raul Saavedra, Bonn Germany
 * v1. 2023-02-05
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
#include "fxp_constants.h"

#define DASHES "========================\n"

#define STR_TEST1 "0000000000000000011111100000000000"
#define STR_TEST2 "0000000000000000100100000000000000"
#define STR_TEST3 "0000000000000001111100000000000000"

//                  ----*----1----*----2----*----3----*----4----*----5----*----6----*
#define STR_E_BIN "1010110111111000010101000101100010100010101110110100101010011010101011111101110001010110001000000010011100111101"
#define STR_E_DEC "2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274274663923"

#define STR_PI_BIN "11001001000011111101101010100010001000010110100011000010001101001100010011000110011000101000101110000000110111"
#define STR_PI_DEC "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214"

#define STR_SQRT2_BIN "101101010000010011110011001100111111100111011110011001001000010001011001011111011000100110110011011"
#define STR_SQRT2_DEC "1.41421356237309504880168872420969807856967187537694807317667973799073247846210703885038753432764157"

// log2(e)
#define STR_LOG2E_DEC "1.4426950408889634073599246810019"

// ln(2)
#define STR_LN2_DEC ".69314718055994530941723212145818"


/*
 * Return binary expansion (up to nbits) of a string representation
 * of a binary number, with last bit either verbatim or rounded depending
 * on the following bit
 */
unsigned long long bex_from_bin(char * pbinnum, int nbits, int rounded)
{
        char * p = pbinnum;
        unsigned long long bnum = 0;
        if (nbits <= 0) nbits = 1;
        int i = 0;
        while (i < nbits) {
            char cbit = *p;
            printf("%c", cbit);
            bnum = (bnum << 1) | ((cbit == '1')? 1: 0);
            i++;
            p++;
        }
        if (rounded) {
            printf("  (Next bit:");
            int nextbit = 0;
            char cbit = *p;
            printf("%c", cbit);
            int rbit = ((cbit  == '1')? 1: 0);
            printf(") -> adding %d", rbit);
            bnum += rbit;
        }
        return bnum;
}

unsigned long long bex_from_dec(char * pdecnum, int nbits, int rounded)
{
        char * p = pdecnum;
        if (nbits <= 0) nbits = 1;
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
        printf("Whole part is d:%d (%d dec digits)\n", wnum, nwd);
        long long binnum = 0;
        int nb = 0;
        while (wnum > 0) {
            int r = wnum % 2;
            //printf("Appending a %d to the whole part\n", r);
            binnum = (r << nb) | binnum;
            nb++;
            wnum /= 2;
        }
        //printf(" (x:%llx. , %d bits)\n", binnum, nb);

        // Process fractional part
        p++; // skip the point
        int ndd = nwd;
        int maxd = nbits / 2 + 1;
        long double ftens = 1.0/10.0;
        long double frac = 0;
        while ((ndd < maxd) && (*p != '\0')) {
            char cdigit = *p;
            frac += ((int) cdigit - '0') * ftens;
            ftens /= 10;
            ndd++;
            p++;
        }
        printf("Fraction part is %Lf \n", frac);
        long double m;
        while (nb < nbits) {
            m = frac * 2.0;
            if (m >= 1.0) {
                binnum = ((binnum << 1) | 1);
                frac = m - 1.0;
                //printf("Appending a 1\n");
            } else {
                binnum = (binnum << 1);
                frac = m;
                //printf("Appending a 0\n");
            }
            nb++;
        }
        printf("\nFull bin number (not rounded): d:%lld, x:%llx\n", binnum, binnum);
        if (rounded) {
            m = frac * 2.0;
            int rbit = ((m >= 1.0)? 1: 0);
            printf("  (Next bit: %d)\n", rbit);
            binnum += rbit;
        }
        return binnum;
}

int main(void)
{
        printf("%sFXP Constants\n%s", DASHES, DASHES);

        printf("\ntest1 %llx\n", bex_from_bin(STR_TEST1, 16, 1));
        printf("\ntest2 %llx\n", bex_from_bin(STR_TEST2, 16, 1));
        printf("\ntest3 %llx\n", bex_from_bin(STR_TEST3, 16, 1));

        printf("\ne  (31 bits) %llx\n\n", bex_from_bin(STR_E_BIN, 31, 1));
        printf("\ne  (63 bits) %llx\n\n", bex_from_bin(STR_E_BIN, 63, 1));

        printf("\npi (31 bits) %llx\n\n", bex_from_bin(STR_PI_BIN, 31, 1));
        printf("\npi (63 bits) %llx\n\n", bex_from_bin(STR_PI_BIN, 63, 1));

        printf("\nlog2(e) (31 bits) %llx\n\n", bex_from_dec(STR_LOG2E_DEC, 31, 1));
        printf("\nlog2(e) (63 bits) %llx\n\n", bex_from_dec(STR_LOG2E_DEC, 63, 1));

        printf("\nln(2) (31 bits) %llx\n\n", bex_from_dec(STR_LN2_DEC, 31, 1));
        printf("\nln(2) (63 bits) %llx\n\n", bex_from_dec(STR_LN2_DEC, 63, 1));

        printf("\nsqrt2 (31 bits) %llx\n\n", bex_from_bin(STR_SQRT2_BIN, 31, 1));
        printf("\nsqrt2 (63 bits) %llx\n\n", bex_from_bin(STR_SQRT2_BIN, 63, 1));

        return 0;
}
