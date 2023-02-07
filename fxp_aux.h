/*
 * fxp_aux.h
 *
 * Utility functions built to test and debug
 * the fxp code, together with the tester itself,
 * and the execution time (xtime) auxiliary program.
 *
 * By Raul Saavedra, 2023-01-28, Bonn Germany
 */

#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

long double int_to_frac(long frac_value);

long double lim_frac(long double x, int fbp);

long double dfxp(int fxp);

long double dfxp_undef();

long double dfxp_pos_inf();

long double dfxp_neg_inf();

void print_int_as_bin(int n, int width);

void print_fxp_as_bin(int x, int width);

void print_fxp(int fxp);

void print_fxp_div(int startmask, int nmaskbits, int n, int frac_bits);

void trace_fxp_div( char * msg,
            int iteration, int frac_bits, int bindex, \
            int dividend, int divisor, \
            int quotient, int newqbit, \
            int mc, \
            int difference);

int fxp_from_ldouble(long double x);

