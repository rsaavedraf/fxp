/* SPDX-License-Identifier: MIT */
/*
 * ulongy.h
 * Emulating an unsigned long by using a couple of uints in a struct
 */
#include "fxp_extern.h"

typedef struct ulongy {
        unsigned int hi;
        unsigned int lo;
} ulongy;

static const ulongy ULONGY_ZERO = { 0u, 0u };

ulongy ulongy_create(unsigned int a, unsigned int b);
unsigned int ulongy_get_hi(ulongy x);
unsigned int ulongy_get_lo(ulongy x);
unsigned int ulongy_hi_to_uint_rounded(ulongy x);
int ulongy_compare(ulongy x, ulongy y);
ulongy ulongy_add(ulongy x, ulongy y);
ulongy ulongy_add_uint(ulongy x, unsigned int y);
ulongy ulongy_rshift(ulongy x, unsigned int rshift);
ulongy ulongy_lshift(ulongy x, unsigned int lshift);
void ulongy_rshift_inplace(ulongy *x, unsigned int rshift);
void ulongy_lshift_inplace(ulongy *x, unsigned int lshift);

ulongy ulongy_from_dmul(unsigned int a, unsigned int b);
unsigned int dmul_into_uint(unsigned int a, unsigned int b);
ulongy dmul_ulongy_x_uint(ulongy x, unsigned int y);
ulongy dmul_ulongys(ulongy x, ulongy y);



