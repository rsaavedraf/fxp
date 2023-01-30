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
            int quotient, int lastqbit, int qbits, \
            int mc, \
            int newqbit, \
            int difference);

// Ones[n] is the lowest +int with n consecutive 1's
// Might be used to speed up some mask operations
static const int ONES[] = {
        0,
        1, 3, 7, 15,
        31, 63, 127, 255,
        511, 1023, 2047, 4095,
        8191, 16383, 32767, 65535,
        131071, 262143, 524287, 1048575,
        2097151, 4194303, 8388607, 16777215,
        33554431, 67108863, 134217727, 268435455,
        536870911, 1073741823, 2147483647 };
