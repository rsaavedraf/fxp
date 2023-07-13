
/* SPDX-License-Identifier: MIT */
/*
 * trig.c
 *
 * Initial implementations of trigonometric functions using long doubles.
 * Written to test the CORDIC algorithm implementation,
 * before tailoring it for fxp's
 *
 * For more details about the CORDIC algorithm:
 * https://en.wikipedia.org/wiki/CORDIC
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

#define DASHES "=========================================\n"

#define SET_RAND_SEED 1

const long double PI_AS_LD = acosl(-1.0L);
const long double HALF_PI_AS_LD = PI_AS_LD / 2.0L;
const long double RAD_TO_GRAD = 180.0L / PI_AS_LD;
const long double GRAD_TO_RAD = PI_AS_LD / 180.0L;

const int MAX_NUMS = 500;

const int LONG_BITS = sizeof(long) * 8;

const long double LOG2_10 = log2l(10.0L);

// ANGLE[n] = atanl(2^-n), with n from 0 to 63
// These angles all in radians, 1st element is pi/4
static const long double ANGLE[] = {
0.78539816339744830961566084581987569937L,  // [0]
0.46364760900080611621425623146121439713L,
0.24497866312686415417208248121127580642L,
0.12435499454676143503135484916387102417L,
0.06241880999595734847397911298550511298L,
0.03123983343026827625371174489249097774L,
0.01562372862047683080280152125657031928L,
0.00781234106010111129646339184219928134L,
0.00390623013196697182762866531142438715L,
0.00195312251647881868512148262507671385L,
0.00097656218955931943040343019971729082L,
0.00048828121119489827546923962564484868L,
0.00024414062014936176401672294325965999L,
0.00012207031189367020423905864611795630L,
0.00006103515617420877502166256917382915L,
0.00003051757811552609686182595343853602L,  // [15]
0.00001525878906131576210723193581269788L,  // [16]
0.00000762939453110197026338848234010509L,
0.00000381469726560649628292307561637299L,
0.00000190734863281018703536536930591724L,
0.00000095367431640596087942067068992311L,
0.00000047683715820308885992758382144925L,
0.00000023841857910155798249094797721893L,
0.00000011920928955078068531136849713792L,
0.00000005960464477539055441392106214179L,
0.00000002980232238769530367674013276771L,
0.00000001490116119384765514709251659596L,
0.00000000745058059692382798713656457450L,
0.00000000372529029846191404526707057181L,
0.00000000186264514923095702909588382148L,
0.00000000093132257461547851535573547768L,
0.00000000046566128730773925777884193471L,  // [31]
0.00000000023283064365386962890204274184L,  // [32]
0.00000000011641532182693481445259909273L,
0.00000000005820766091346740722649676159L,
0.00000000002910383045673370361327303270L,
0.00000000001455191522836685180663959784L,
0.00000000000727595761418342590332018410L,
0.00000000000363797880709171295166014020L,
0.00000000000181898940354585647583007612L,
0.00000000000090949470177292823791503881L,
0.00000000000045474735088646411895751950L,
0.00000000000022737367544323205947875976L,
0.00000000000011368683772161602973937988L,
0.00000000000005684341886080801486968994L,
0.00000000000002842170943040400743484497L,
0.00000000000001421085471520200371742249L,
0.00000000000000710542735760100185871124L,
0.00000000000000355271367880050092935562L,
0.00000000000000177635683940025046467781L,
0.00000000000000088817841970012523233891L,
0.00000000000000044408920985006261616945L,
0.00000000000000022204460492503130808473L,
0.00000000000000011102230246251565404236L,
0.00000000000000005551115123125782702118L,
0.00000000000000002775557561562891351059L,
0.00000000000000001387778780781445675530L,
0.00000000000000000693889390390722837765L,
0.00000000000000000346944695195361418882L,
0.00000000000000000173472347597680709441L,
0.00000000000000000086736173798840354721L,  // [60]
0.00000000000000000043368086899420177360L,
0.00000000000000000021684043449710088680L,
0.00000000000000000010842021724855044340L   // [63]
};

// (*) Using IEEE-754's long double precision, the precision of
// the calculation of the CORDIC K scaling factor increases
// up to n = 56. For all n > 56, K[n] == K[56] ==
const long double CORDIC_KVALUE = 0.60725293500888125616944675250492850317L;

const unsigned int MAX_LOOPS = 64u;

/*
 * Calculates simultaneously the sine and cosine of the argument x
 * using the CORDIC algorithm:
 * https://en.wikipedia.org/wiki/CORDIC
 */
void my_sincos1(long double x, long double *mysin, long double *mycos, unsigned int loops)
{
        // x should be in [pi/2, -pi/2]
        long double a = 0.0L;
        long double c = 1.0L;
        long double s = 0.0L;
        if (loops > MAX_LOOPS) loops = MAX_LOOPS;
        //printf("Start: angle is %Lf, c is %Lf, s is %Lf\n", a/PI_AS_LD*180.0L, c, s);
        // CORDIC implementation
        for (int i = 0; i < loops; i++) {
                long double angle = ANGLE[i];
                long double tangent = 1.0L / ((long double) (1ul << i));
                if (a < x) {
                        //printf("+ rotation\n");
                        a += angle;
                        // These products would get implemented with simple shifts
                        long double newc = c - s*tangent;
                        s += c*tangent;
                        c = newc;
                } else {
                        //printf("- rotation\n");
                        a -= angle;
                        // These products would get implemented with simple shifts
                        long double newc = c + s*tangent;
                        s -= c*tangent;
                        c = newc;
                }
                //printf("Iteration %2d: angle change %Lf, new c is %Lf, new s is %Lf, angle' %14.12Lf\n", \
                //            i, angle/PI_AS_LD*180.0L, c, s, a/PI_AS_LD*180.0L);
        }
        *mycos = c * CORDIC_KVALUE;
        *mysin = s * CORDIC_KVALUE;
        return;
}

/*
 * Alternative CORDIC implementation avoids the final multiplication with the scaling factor
 * by directly starting from point (CORDIC_KVALUE, 0) instead of from (1, 0)
 */
void my_sincos2(long double x, long double *mysin, long double *mycos, unsigned int loops)
{
        // x should be in [pi/2, -pi/2]
        long double a = 0.0L;
        long double c = CORDIC_KVALUE;
        long double s = 0.0L;
        if (loops > MAX_LOOPS) loops = MAX_LOOPS;
        //printf("Start: angle is %Lf, c is %Lf, s is %Lf\n", a/PI_AS_LD*180.0L, c, s);
        // CORDIC implementation
        for (int i = 0; i < loops; i++) {
                long double angle = ANGLE[i];
                long double tangent = 1.0L / ((long double) (1ul << i));
                if (a < x) {
                        //printf("+ rotation\n");
                        a += angle;
                        // These products would get implemented with simple shifts
                        long double newc = c - s*tangent;
                        s += c*tangent;
                        c = newc;
                } else {
                        //printf("- rotation\n");
                        a -= angle;
                        // These products would get implemented with simple shifts
                        long double newc = c + s*tangent;
                        s -= c*tangent;
                        c = newc;
                }
                //printf("Iteration %2d: angle change %Lf, new c is %Lf, new s is %Lf, angle' %14.12Lf\n", \
                //            i, angle/PI_AS_LD*180.0L, c, s, a/PI_AS_LD*180.0L);
        }
        *mycos = c;
        *mysin = s;
        return;
}

int main(void)
{
        if (SET_RAND_SEED) srand((unsigned int) time(0));

        printf("\n%sTests the trigonometric functions using long doubles\n%s", DASHES, DASHES);
        printf("K scaling factor for Cordic: %41.38Lf\n", CORDIC_KVALUE);

        long double errors1[] = {0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, \
                                 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L};
        long double errors2[] = {0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, \
                                 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L};
        int n;
        for (int i = 0; i < MAX_NUMS; i++) {
                int sign = ((rand() % 2) == 1)? -1: 1;
                const long double angle = sign * HALF_PI_AS_LD * ((long double) rand() / INT_MAX);
                //const long double angle = 44.95L / 180.0L * PI_AS_LD;
                if (i % 50 == 0)
                        printf("Angle is: %9.6Lf rad (%9.6Lf °)\n", angle, angle * RAD_TO_GRAD);
                const long double s = sinl(angle);
                const long double c = cosl(angle);
                n = 0;
                for (unsigned int loops = 8; loops < 65; loops += 8) {
                        long double mysin1 = 0.0L, mycos1 = 0.0L;
                        long double mysin2 = 0.0L, mycos2 = 0.0L;
                        my_sincos1(angle, &mysin1, &mycos1, (unsigned int) loops);
                        my_sincos2(angle, &mysin2, &mycos2, (unsigned int) loops);
                        // Accumulate the largest errors for method 1 and 2
                        long double errsin1 = fabsl(mysin1 - s);
                        long double errcos1 = fabsl(mycos1 - c);
                        long double maxerr1 = (errsin1 > errcos1)? errsin1: errcos1;
                        if (maxerr1 > errors1[n]) errors1[n] = maxerr1;
                        long double errsin2 = fabsl(mysin2 - s);
                        long double errcos2 = fabsl(mycos2 - c);
                        long double maxerr2 = (errsin2 > errcos2)? errsin2: errcos2;
                        if (maxerr2 > errors2[n]) errors2[n] = maxerr2;
                        if ((loops == 32) && (i % 50 == 0)) {
                                printf("\tUsing %u Cordic loops: \n", loops);
                                printf("\t\tclib sin(%9.6Lf) = %41.38Lf\n", angle, s);
                                printf("\t\tmy  sin1(%9.6Lf) = %41.38Lf (delta: %5.3LE)\n", angle, mysin1, errsin1);
                                printf("\t\tmy  sin2(%9.6Lf) = %41.38Lf (delta: %5.3LE)\n", angle, mysin2, errsin2);
                                printf("\t\tclib cos(%9.6Lf) = %41.38Lf\n", angle, c);
                                printf("\t\tmy  cos1(%9.6Lf) = %41.38Lf (delta: %5.3LE)\n", angle, mycos1, errcos1);
                                printf("\t\tmy  cos2(%9.6Lf) = %41.38Lf (delta: %5.3LE)\n", angle, mycos2, errcos2);
                        }
                        n++;
                }
        }
        n = 0;
        for (unsigned int loops = 8; loops < 65; loops += 8) {
                long double err1 = errors1[n];
                long double err2 = errors2[n];
                long double maxerr = (err1 > err2)? err1: err2;
                long double decprec = fabsl(log10l(maxerr));
                printf("For %2u Cordic loops, max error: %5.3LE  ", loops, maxerr);
                printf("-> Dec. Prec: %5.2Lf, Bits of precision: %5.2Lf\n", \
                        decprec, decprec * LOG2_10);
                n++;
        }
}
