/*
 * fxp_constants.c
 *
 * Auxiliary program to generate fxp versions of some
 * important/transcendental constants
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

//                 100+ bits of E
//                 ----*----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8----*----9----*----0----*----1--
#define STR_E_BIN "1010110111111000010101000101100010100010101110110100101010011010101011111101110001010110001000000010011100111101"
#define STR_E_DEC "2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274274663923"

//                 100+ bits of PI
#define STR_PI_BIN "11001001000011111101101010100010001000010110100011000010001101001100010011000110011000101000101110000000110111"
#define STR_PI_DEC "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214"

#define STR_SQRT2_BIN "101101010000010011110011001100111111100111011110011001001000010001011001011111011000100110110011011"
#define STR_SQRT2_DEC "1.41421356237309504880168872420969807856967187537694807317667973799073247846210703885038753432764157"

// log2(e)
#define STR_LOG2E_DEC "1.4426950408889634073599246810019"

// ln(2)
#define STR_LN2_DEC ".69314718055994530941723212145818"

#define ZERO 1.0E-124


// Decimal values from table A_2 in:
// https://en.wikipedia.org/wiki/BKM_algorithm
// The values here are simply: a[k] = log2(1 + 0.5^k)
static char * alog2[] = {
// 76 decimal digits for each number here
//0----*----1----*----2----*----3----*----4----*----5----*----6----*----7----*-
"1.0000000000000000000000000000000000000000000000000000000000000000000000000000",
"0.5849625007211561814537389439478165087598144076924810604557526545410982276485",
"0.3219280948873623478703194294893901758648313930245806120547563958159347765589",
"0.1699250014423123629074778878956330175196288153849621209115053090821964552970",
"0.0874628412503394082540660108104043540112672823448206881266090643866965081686",
"0.0443941193584534376531019906736094674630459333742491317685543002674288465967",
"0.0223678130284545082671320837460849094932677948156179815932199216587899627785",
"0.0112272554232541203378805844158839407281095943600297940811823651462712311786",
"0.0056245491938781069198591026740666017211096815383520359072957784732489771013",
"0.0028150156070540381547362547502839489729507927389771959487826944878598909400",
"0.0014081943928083889066101665016890524233311715793462235597709051792834906001",
"0.0007042690112466432585379340422201964456668872087249334581924550139514213168",
"0.0003521774803010272377989609925281744988670304302127133979341729842842377649",
"0.0001760994864425060348637509459678580940163670081839283659942864068257522373",
"0.0000880524301221769086378699983597183301490534085738474534831071719854721939",
"0.0000440268868273167176441087067175806394819146645511899503059774914593663365",
"0.0000220136113603404964890728830697555571275493801909791504158295359319433723",
"0.0000110068476674814423006223021573490183469930819844945565597452748333526464",
"0.0000055034343306486037230640321058826431606183125807276574241540303833251704",
"0.0000027517197895612831123023958331509538486493412831626219340570294203116559",
"0.0000013758605508411382010566802834037147561973553922354232704569052932922954",
"0.0000006879304394358496786728937442939160483304056131990916985043387874690617",
"0.0000003439652607217645360118314743718005315334062644619363447395987584138324",
"0.0000001719826406118446361936972479533123619972434705828085978955697643547921",
"0.0000000859913228686632156462565208266682841603921494181830811515318381744650",
"0.0000000429956620750168703982940244684787907148132725669106053076409624949917",
"0.0000000214978311976797556164155504126645192380395989504741781512309853438587",
"0.0000000107489156388827085092095702361647949603617203979413516082280717515504",
"0.0000000053744578294520620044408178949217773318785601260677517784797554422804",
"0.0000000026872289172287079490026152352638891824761667284401180026908031182361",
"0.0000000013436144592400232123622589569799954658536700992739887706412976115422",
"0.0000000006718072297764289157920422846078078155859484240808550018085324187007",
};

static const long double BKM_LOGS[] = {
// Skipping first one with a one, we won't really use it in our bkm,
1.0000000000000000000000000000000000000000000000000000000000000000000000000000,
0.5849625007211561814537389439478165087598144076924810604557526545410982276485,
0.3219280948873623478703194294893901758648313930245806120547563958159347765589,
0.1699250014423123629074778878956330175196288153849621209115053090821964552970,
0.0874628412503394082540660108104043540112672823448206881266090643866965081686,
0.0443941193584534376531019906736094674630459333742491317685543002674288465967,
0.0223678130284545082671320837460849094932677948156179815932199216587899627785,
0.0112272554232541203378805844158839407281095943600297940811823651462712311786,
0.0056245491938781069198591026740666017211096815383520359072957784732489771013,
0.0028150156070540381547362547502839489729507927389771959487826944878598909400,
0.0014081943928083889066101665016890524233311715793462235597709051792834906001,
0.0007042690112466432585379340422201964456668872087249334581924550139514213168,
0.0003521774803010272377989609925281744988670304302127133979341729842842377649,
0.0001760994864425060348637509459678580940163670081839283659942864068257522373,
0.0000880524301221769086378699983597183301490534085738474534831071719854721939,
0.0000440268868273167176441087067175806394819146645511899503059774914593663365,
0.0000220136113603404964890728830697555571275493801909791504158295359319433723,
0.0000110068476674814423006223021573490183469930819844945565597452748333526464,
0.0000055034343306486037230640321058826431606183125807276574241540303833251704,
0.0000027517197895612831123023958331509538486493412831626219340570294203116559,
0.0000013758605508411382010566802834037147561973553922354232704569052932922954,
0.0000006879304394358496786728937442939160483304056131990916985043387874690617,
0.0000003439652607217645360118314743718005315334062644619363447395987584138324,
0.0000001719826406118446361936972479533123619972434705828085978955697643547921,
0.0000000859913228686632156462565208266682841603921494181830811515318381744650,
0.0000000429956620750168703982940244684787907148132725669106053076409624949917,
0.0000000214978311976797556164155504126645192380395989504741781512309853438587,
0.0000000107489156388827085092095702361647949603617203979413516082280717515504,
0.0000000053744578294520620044408178949217773318785601260677517784797554422804,
0.0000000026872289172287079490026152352638891824761667284401180026908031182361,
0.0000000013436144592400232123622589569799954658536700992739887706412976115422,
0.0000000006718072297764289157920422846078078155859484240808550018085324187007,
};

static const unsigned long BKM_LOGSX[] = {
0x000000,   0x7ef5a2,   0xcd1b8b,   0xfdeb43,   0xb244c5,
0xb1dfb0,   0x5c2d22,   0xb77842,   0xb5eab1,   0xfe14ea,
0x3743f4,   0x17bcab,   0xb0a86f,   0x623f4a,   0x93b17d,
0x828017,   0xef6a3e,   0x833fb7,   0x448283,   0xa2f9eb,
0x51ab20,   0xa8e11a,   0x547370,   0xaa3a70,   0x551d66,
0x2a8ebe,   0x154762,   0x8aa3b1,   0xc551d9,   0xe2a8ec,
0x715476,   0xb8aa3b,
};



/*
 * Return binary expansion (up to nbits) of a string representation
 * of a binary number, with last bit either verbatim or rounded
 */
unsigned long long bex_from_bin(char * pbinnum, int wbits, int fbits, int rounded)
{
        char * p = pbinnum;
        unsigned long long bnum = 0;
        if (wbits <= 0) wbits = 31;
        if (fbits <= 0) fbits = 32;
        int roomleft = wbits + fbits;
        while ((roomleft > 0) && (*p != '\0')) {
            char cbit = *p;
            //printf("%c", cbit);
            bnum = (bnum << 1) | ((cbit == '1')? 1: 0);
            roomleft--;
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
        if (wbits <= 0) wbits = 31;
        if (fbits <= 0) fbits = 32;
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
        long double tenth = 1.0/10.0;
        long double ftens = tenth;
        while (*p != '\0') {
            char cdigit = *p;
            num += ((int) cdigit - '0') * ftens;
            ftens *= tenth;
            p++;
        }
        return num;
}

#define LONGLONG_BITS ((int) (sizeof(long long) * 8))

// Calculate the log2 of a long double number
// (Implementing and testing the logarithm algorithm first
// with floating points, before implementing it for fxp's)
//
// Based on the general algorithm to calculate binary logarithms
// explained by Clay Turner in IEEE Signal Processing Magazine,
// Sep/2010. His short article has D. E. Knuth's "The Art of
// Computer Programming Vol 2: Seminumerical Algorithms",
// 2nd ed., 1981 (pages 441 - 446) as the only reference.
// Additional useful resources: CORDIC and/or BKM algorithms
long double my_log2(long double x)
{
        if (x < 0.0) return dfxp_undef();
        if (x <= ZERO) return dfxp_neg_inf();
        long double z = x;
        int c = 0; // characteristic
        while (z >= 2.0) {
            c++;
            z /= 2;
        }
        while (z < 1.0) {
            c--;
            z *= 2;
        }
        // Here we have already calculated the log characteristic c, and
        // we have z satisfying: 1 <= z < 2, so we can use it to calculate
        // the mantissa
        //printf("log c: %d\n", c);
        long double m = 0; // mantissa
        long double b = 0.5;
        int nb = 60; // desired number of mantissa bits to process
        while (nb > 0) {
            z = z * z;
            if (z >= 2.0) {
                z /= 2;
                m += b;
            }
            b = b/2;
            nb--;
        }
        // Here we have already calculated the mantissa down to
        // nb bits of precision
        //printf("log m: %Lf\n", m);
        long double full_log = ((long double) c) + m;
        //printf("log(%Lf) = c + m = %Lf\n", x, full_log);
        return full_log;
}

// Implementation of the bkm algorithm:
// https://en.wikipedia.org/wiki/BKM_algorithm
long double my_log2_bkm(long double x, int nbits)
{
        if (x < 0.0) return dfxp_undef();
        if (x <= ZERO) return dfxp_neg_inf();
        long double z = x;
        int c = 0; // characteristic
        while (z >= 2.0) {
            c++;
            z /= 2;
        }
        while (z < 1.0) {
            c--;
            z *= 2;
        }
        // Here we have already calculated the log characteristic c, and
        // we have z satisfying: 1 <= z < 2
        long double xx = 1.0;
        long double yy = 0.0;
        // Notice ss starting with 0.5 because our first table entry is not
        // the 1.00 but the next one. Skipping the 1.0 because we know
        // zz < 2, so that first test if (zz < z) when ss == 1.0 would be
        // false for sure
        long double ss = 0.5;
        long double zz;
        for (int k = 0; k < nbits; k++) {
                zz = xx + xx * ss;
                if (zz < z) {
                        xx = zz;
                        yy += BKM_LOGS[k];
                }
                ss *= 0.5;
        }
        long double full_log = ((long double) c) + yy;
        //printf("log(%Lf) = c + m = %Lf\n", x, full_log);
        return full_log;
}



int main(void)
{
        printf("%sFXP Constants\n%s", DASHES, DASHES);

        //printf("\ntest1 %llx\n", bex_from_bin(STR_TEST1, 16, 1));
        //printf("\ntest2 %llx\n", bex_from_bin(STR_TEST2, 16, 1));
        //printf("\ntest3 %llx\n", bex_from_bin(STR_TEST3, 16, 1));

        printf("e  (31 bits) %llx\n", bex_from_bin(STR_E_BIN, 2, 29, 1));
        printf("e  (63 bits) %llx\n\n", bex_from_bin(STR_E_BIN, 2, 61, 1));

        printf("pi (31 bits) %llx\n", bex_from_bin(STR_PI_BIN, 2, 29, 1));
        printf("pi (63 bits) %llx\n\n", bex_from_bin(STR_PI_BIN, 2, 61, 1));

        printf("log2(e) (31 bits) %llx\n", bex_from_dec(STR_LOG2E_DEC, 1, 30, 1));
        printf("log2(e) (63 bits) %llx\n\n", bex_from_dec(STR_LOG2E_DEC, 1, 62, 1));

        printf("ln(2) (31 bits) %llx\n", bex_from_dec(STR_LN2_DEC, 0, 31, 1));
        printf("ln(2) (63 bits) %llx\n\n", bex_from_dec(STR_LN2_DEC, 0, 63, 1));

        //printf("\nsqrt2 (31 bits) %llx\n\n", bex_from_bin(STR_SQRT2_BIN, 31, 1));
        //printf("\nsqrt2 (63 bits) %llx\n\n", bex_from_bin(STR_SQRT2_BIN, 63, 1));

        // Generate fxp's for the logs table of numbers to be used by BKM algorithm
        int nfbits = 62;
        printf("\nPrecalculated logs for BKM (frac bits: %d):\n",
                nfbits);
        fxp_set_frac_bits(nfbits);
        int n = sizeof(alog2) / sizeof(alog2[0]);
        int j = 1;
        long sumlogs2 = 0;
        for (int i=0; i < n; i++) {
            //printf("log2(1 + 0.5^%i) is ", i);
            //printf(alog2[i]);
            //printf("\n");
            //printf("In binary: \t");
            unsigned long long num = bex_from_dec(alog2[i], 0, nfbits, 0);
            //unsigned int x = (unsigned int) num;
            printf("0x%llX,\n", num);
            printf("0x%lX,\n\n", BKM_LOGSX[i]);
            sumlogs2 += BKM_LOGSX[i];
            if (j % 5 == 0) printf("\n");
            j++;
        }
        printf("Sum BKM_LOGSX is:%ld\n", sumlogs2);
        printf("Shifted 24 bits:%ld\n", (sumlogs2 >> 24));
        /*
        long double bkm_logs[n];
        for (int i=0; i < n; i++) {
            bkm_logs[i] = Lf_from_dec(alog2[i]);
            //int x = (int) num;
            //printf("0x%x,\n", x);
            printf("%5.24Lf\n", bkm_logs[i]);
        }
        */
        long double nums[] = {8.0, 1025.55, 0.5, 0.01};
        long double invln2 = 1 / log(2);
        //printf("Size of long long is %d\n", LONGLONG_BITS);
        n = sizeof(nums) / sizeof(nums[0]);
        for (int i=0; i < n; i++) {
                printf("\nlg2(%Lf) = \t%5.20Lf\n",
                        nums[i],
                        log(nums[i]) * invln2);
                printf("my_Lg2(%Lf) = \t%5.20Lf\n",
                        nums[i],
                        my_log2(nums[i]));
                printf("my_Lg2_bkm(%Lf) =\t%5.20Lf\n",
                        nums[i],
                        my_log2_bkm(nums[i], 32));
        }

        return 0;
}
