/* SPDX-License-Identifier: MIT */
/*
 * print_as_bits.h
 *
 * Utility functions to print ints, longs, and long doubles as binaries
 *
 * By Raul Saavedra, 2023-03-20, Bonn Germany
 */

void print_int_as_bin(int n, int width);
void print_long_as_bin(long n);
void print_uint_as_bin(unsigned int n);
void print_ulong_as_bin(unsigned long n);
void inspect_long_double(long double x);
unsigned long get_ulong_bits_from_ldouble(long double x);
