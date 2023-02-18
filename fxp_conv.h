/*
 * fxp_conv.h
 *
 * Functions for convertions between fxp and floating-point types
 *
 */
#include "fxp.h"
#include <float.h>
#include <math.h>

extern int FXP_frac_bits;
extern float FXP_max_f;
extern float FXP_min_f;
extern double FXP_max_d;
extern double FXP_min_d;
extern long double FXP_max_ld;
extern long double FXP_min_ld;

// Convertions between fxp and floating point types
long double fxp2ld(int fxp1);
double fxp2d(int fxp1);
float fxp2f(int fxp1);
int f2fxp(float x);
int d2fxp(double x);
int ld2fxp(long double x);

