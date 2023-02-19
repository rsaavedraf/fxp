/*
 * fxp_aux.h
 *
 * Utility functions for the fxp testing program,
 * and for the execution time (xtime) auxiliary program.
 *
 * By Raul Saavedra, 2023-01-28, Bonn Germany
 */

#include "fxp.h"
#include "fxp_conv.h"

extern const float FXP_UNDEF_F;
extern const float FXP_NINF_F;
extern const float FXP_PINF_F;
extern const double FXP_UNDEF_D;
extern const double FXP_NINF_D;
extern const double FXP_PINF_D;
extern const long double FXP_UNDEF_LD;
extern const long double FXP_NINF_LD;
extern const long double FXP_PINF_LD;
extern int FXP_frac_bits;
extern int FXP_frac_max;

long double lim_frac(long double x, int fbp);
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
void print_sys_info();
