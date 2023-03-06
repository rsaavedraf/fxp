/* SPDX-License-Identifier: MIT */
/*
 * fxp_conv.h
 *
 * Functions for convertions between fxp and floating-point types
 *
 */
#include "fxp.h"
#include "fxp_extern.h"
#include <float.h>
#include <math.h>

// Convertions between fxp and floating point types
long double fxp2ld(int fxp1);
double fxp2d(int fxp1);
float fxp2f(int fxp1);
int f2fxp(float x);
int d2fxp(double x);
int ld2fxp(long double x);

