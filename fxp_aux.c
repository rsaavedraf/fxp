/*
 * fxp_aux.c
 *
 * By Raul Saavedra, 2023-01-28, Bonn Germany
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

long double dfxp_max_ld()
{
    return dfxp(FXP_MAX);
}

long double dfxp_min_ld()
{
    return dfxp(FXP_MIN);
}

/*
 * Returns the exact fraction value from an fxp as
 * a long double
 */
long double int_to_frac(long frac_value)
{
    long double bfsum = 0.0;
    long twopower = pow(2, fxp_get_frac_bits());
    long abs_frac = (frac_value < 0)? -frac_value: frac_value;
    while (abs_frac > 0) {
        long bit = abs_frac & 1;
        if (bit)
            bfsum += ((long double) 1.0) / twopower;
        abs_frac = abs_frac >> 1;
        twopower = twopower >> 1;
    }
    return (frac_value < 0)? -bfsum: bfsum;
}

/*
 * Limits the precision of a double to the specified
 * number of fractional bits. This eases the
 * comparisons between results from (precision limited)
 * fxp operations vs. long double operations used
 * as the reference
 */
long double lim_frac(long double x, int fbp)
{
    if (x == FXP_UNDEF_LD) return FXP_UNDEF_LD;
    if (x < dfxp_min_ld()) return FXP_NEG_INF_LD;
    if (x > dfxp_max_ld()) return FXP_POS_INF_LD;
    long shift = pow(2, fbp);
    long double px = (x < 0)? -x: x;
    long pxwhole = trunc(px);
    long double dfrac = px - ((long double) pxwhole);
    long lfrac = trunc(dfrac * shift);
    //printf("Whole part: %ld, frac part: %ld (fbp:%d, shift:%ld)\n", pxwhole, lfrac, fbp, shift);
    long double bfsum = int_to_frac( lfrac );
    long double xpreclim;
    if (x < 0) {
        xpreclim = -pxwhole - bfsum;
    } else {
        xpreclim = pxwhole + bfsum;
    }
    //printf("Orig. long double: %Lf\n", x);
    //printf("Limited frac prec: %Lf (%d frac bits)\n", xpreclim, fbp);
    return xpreclim;
}

/*
 * Creates an fxp given a long double value
 */
int fxp_from_ldouble(long double x)
{
    if (x <= FXP_UNDEF_LD) return FXP_UNDEF;
    if (x < dfxp_min_ld()) return FXP_NEG_INF;
    if (x > dfxp_max_ld()) return FXP_POS_INF;
    long shift = pow(2, fxp_get_frac_bits());
    long double px = (x < 0)? -x: x;
    long plwhole = trunc(px);
    long double pdfrac = px - ((long double) plwhole);
    long plfrac = trunc(pdfrac * shift);
    if (plwhole > 0)
        return fxp_bin((int) plwhole, (int) plfrac);
    else if (plwhole < 0)
        return fxp_bin((int) -plwhole, (int) plfrac);
    else
        return fxp_bin(0, (int) -plfrac);
}


/*
 * Returns a long double value corresponding to a given fxp
 */
long double dfxp(int fxp)
{
    if (fxp == FXP_UNDEF) return FXP_UNDEF_LD;
    if (fxp == FXP_NEG_INF) return FXP_NEG_INF_LD;
    if (fxp == FXP_POS_INF) return FXP_POS_INF_LD;
    long wi = fxp_get_whole_part(fxp);
    long fi = fxp_get_bin_frac(fxp);
    //printf("\nfxp long  whole & frac parts: %ld, %ld\n", wi, fi);
    long double wd = (long double) wi;
    long double fd = int_to_frac(fi);
    //printf("fxp Lf whole & frac parts: %Lf, %Lf\n", wd, fd);
    long double x = wd + fd;
    return x;
}

void print_int_as_bin(int n, int width)
{
    int an;
    int signspace = 0;
    if (n < 0) {
        an = -n;
        signspace = 1;
    } else {
        an = n;
    }
    int nbn = fxp_nbits(an);
    int margin = (nbn == 0? 1: nbn) - signspace;
    while (width > margin) {
        printf(" ");
        width--;
    }
    if (signspace) printf("-");
    nbn = (nbn == 0)? 1: nbn;
    int i = nbn;
    while (i > 0) {
        int bit = (an & (1 << (i - 1))) >> (i - 1);
        printf("%d", bit);
        i--;
    }
}

void print_fxp_as_bin(int n, int width)
{
    int an;
    int signspace = 0;
    if (n < 0) {
        an = -n;
        signspace = 1;
    } else {
        an = n;
    }
    int nbn = fxp_nbits(an);
    int margin = (nbn == 0? 1: nbn) - signspace;
    while (width > margin) {
        printf(" ");
        width--;
    }
    if (signspace) printf("(-)");
    int frbits = fxp_get_frac_bits();
    int i = MAX(frbits, (nbn == 0)? 1: nbn);
    while (i > 0) {
        if (i == frbits) printf(".");
        if (i <= nbn) {
            int bit = (an & (1 << (i - 1))) >> (i - 1);
            printf("%d", bit);
        } else {
            printf("0");
        }
        i--;
    }
}

/*
 * Prints out an fxp
 */
void print_fxp(int fxp)
{
        if (fxp == FXP_POS_INF || fxp == FXP_NEG_INF || fxp == FXP_UNDEF) {
                printf(fxp==FXP_UNDEF? "UNDEF":
                          (fxp==FXP_POS_INF? "+INF": "-INF"));
                return;
        }
        //printf("\ninput fxp is: %d", fxp);
        long double n = lim_frac(dfxp(fxp), fxp_get_frac_bits());
        //printf("long double n = %Lf\n", n);
        int whole = fxp_get_whole_part(fxp);
        int nbits;
        if (fxp_get_frac_max_dec() <= 9999) {
            if (fxp < 0) {
                int pfxp = -fxp;
                nbits = fxp_nbits(pfxp);
                printf("%d.%5d (=Lf:%Lf = %d = x(-)%x = b",
                        whole,
                        fxp_get_dec_frac(fxp),
                        n, fxp, pfxp);
            } else {
                nbits = fxp_nbits(fxp);
                printf("%d.%5d (=Lf:%Lf = %d = x%x = b",
                        whole,
                        fxp_get_dec_frac(fxp),
                        n, fxp, fxp);
            }
        } else {
            if (fxp < 0) {
                int pfxp = -fxp;
                nbits = fxp_nbits(pfxp);
                printf("%d.%8d (=LE:%1.2LE = %d = x(-)%x = b",
                        whole,
                        fxp_get_dec_frac(fxp),
                        n, fxp, pfxp);
            } else {
                nbits = fxp_nbits(fxp);
                printf("%d.%8d (=LE:%1.2LE = %d = x%x = b",
                        whole,
                        fxp_get_dec_frac(fxp),
                        n, fxp, fxp);
            }
        }
        print_fxp_as_bin(fxp, 0);
        printf(", %d bits)", nbits);
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
        print_fxp_div(bindex - bmc, bmc, dividend, fxp_get_frac_bits());
        printf("    \ty: ");
        print_int_as_bin(divisor, 0);
        printf(" (%d, x%x, %d bits)\n", \
                divisor, (unsigned int) divisor, fxp_nbits(divisor));
    }
    printf(" m:");
    print_int_as_bin(mc, bindex);
    //printf("  (x%x, %d bits) \tq:%d (x%x, %d bits,  last q bit:%d)\n",
    //        (unsigned int) mc, bmc, \
    //        quotient, (unsigned int) quotient, \
    //        fxp_nbits(quotient), lastqbit);
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
        system("hostnamectl | grep -e 'Operating System' -e 'Architecture'");
        system("cat /proc/cpuinfo | grep -e 'CPU' -e 'model name' -e 'Model' | sort -r | head -1");

        printf("\nNum type sizes:\n");
        printf("char        has a size of %zu bytes.\n", sizeof(char));
        printf("int         has a size of %zu bytes.\n", sizeof(int));
        printf("long        has a size of %zu bytes.\n", sizeof(long));
        printf("long long   has a size of %zu bytes.\n", sizeof(long long));
        printf("float       has a size of %zu bytes.\n", sizeof(float));
        printf("double      has a size of %zu bytes.\n", sizeof(double));
        printf("long double has a size of %zu bytes.\n", sizeof(long double));
}
