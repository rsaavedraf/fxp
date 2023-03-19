/* SPDX-License-Identifier: MIT */
/*
 * fxp_aux.h
 *
 * Utility functions for the fxp testing program,
 * and for the execution time (xtime) auxiliary program.
 *
 * By Raul Saavedra, 2023-01-28, Bonn Germany
 */

#include "fxp.h"
#include "fxp_extern.h"
#include "fxp_conv.h"

long double lim_frac(long double x, int fbp);
void print_int_as_bin(int n, int width);
void print_long_as_bin(long n);
void print_uint_as_bin(unsigned int n);
void print_ulong_as_bin(unsigned long n);
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
