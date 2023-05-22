/* SPDX-License-Identifier: MIT */
/*
 * ulongy.c
 * Emulates an unsigned long by using a couple of uints in a struct (ulongy).
 *
 * For systems where sizeof(longs) is either not >= 2 * sizeof(int)),
 * and/or when simply choosing to strictly use ints and only ints.
 *
 * By Raul Saavedra, Bonn, Germany
 */

#include "ulongy.h"
#include "stdio.h"

const unsigned int ALL_ONES = ~0u;

inline ulongy ulongy_create(unsigned int vhi, unsigned int vlo)
{
        ulongy x = { vhi, vlo };
        return x;
}

inline unsigned int ulongy_get_hi(ulongy x)
{
        return x.hi;
}

inline unsigned int ulongy_get_lo(ulongy x)
{
        return x.lo;
}

inline unsigned int ulongy_hi_to_uint_rounded(ulongy x)
{
        return x.hi + (x.lo >> FXP_INT_BITS_M1);
}

// Returns -1 when x < y, 0 when ==, and 1 when x > y
inline int ulongy_compare(ulongy x, ulongy y)
{
    return (x.hi < y.hi)? -1:
                (x.hi > y.hi)? 1:
                    (x.lo == y.lo)? 0: ((x.lo < y.lo)? -1: 1);
}

// Compares the ulongy x to y as if y was the .hi half of
// another ulongy that had lo == 0
inline int ulongy_compare_to_uint(ulongy x, unsigned int y)
{
    return (x.hi < y)? -1:
                (x.hi > y)? 1: (x.lo > 0);
}

inline ulongy ulongy_add(ulongy x, ulongy y)
{
        unsigned int sumlo = x.lo + y.lo;
        unsigned int co = (sumlo < x.lo);
        unsigned int sumhi = x.hi + y.hi + co;
        //printf("x.lo : x%X,  y.lo: x%X\n", x.lo, y.lo);
        //printf("sumlo: x%X,  co: x%X\n", sumlo, co);
        //printf("x.hi : x%X,  y.hi: x%X\n", x.hi, y.hi);
        //printf("sumhi: x%X\n", sumhi);
        ulongy result = { sumhi, sumlo };
        return result;
}

// Returns the ulongy that would correspond to the representation
// of the argument's value negated into two's complement
inline ulongy ulongy_negate(ulongy x)
{
        x.hi = ~x.hi;
        x.lo = ~x.lo;
        return ulongy_add_uint(x, 1u);
}

inline ulongy ulongy_sub(ulongy x, ulongy y)
{
        return ulongy_add(x, ulongy_negate(y));
}

inline ulongy ulongy_add_uint(ulongy x, unsigned int b)
{
        unsigned int sumlo = x.lo + b;
        ulongy result = { x.hi + (sumlo < x.lo), sumlo };
        return result;
}

inline ulongy lshift_ulongy(ulongy x, unsigned int lshift)
{
        if (lshift < FXP_INT_BITS) {
                if (lshift == 0) return x;
                int d = FXP_INT_BITS - lshift;
                ulongy shifted = { (x.hi << lshift) | \
                                        ((x.lo & (ALL_ONES << d)) >> d),
                                   x.lo << lshift };
                return shifted;
        } else {
                if (lshift >= FXP_LONG_BITS) return ULONGY_ZERO;
                ulongy shifted = { x.lo << (lshift - FXP_INT_BITS), 0u};
                return shifted;
        }
}

inline ulongy rshift_ulongy(ulongy x, unsigned int rshift)
{
        if (rshift < FXP_INT_BITS) {
                if (rshift == 0) return x;
                int d = FXP_INT_BITS - rshift;
                ulongy shifted = {  x.hi >> rshift, \
                                    ((x.hi & (ALL_ONES >> d)) << d) | \
                                            (x.lo >> rshift) };
                return shifted;
        } else {
                if (rshift >= FXP_LONG_BITS) return ULONGY_ZERO;
                ulongy shifted = { 0u, x.hi >> (rshift - FXP_INT_BITS) };
                return shifted;
        }
}

inline ulongy rshift_ulongy_rounding(ulongy x, unsigned int rshift)
{
        if (rshift == 0) return x;
        ulongy r = rshift_ulongy(x, rshift - 1);
        unsigned int rbit = r.lo & 1u;
        r = rshift_ulongy(r, 1);
        return ulongy_add_uint(r, rbit);
}

inline unsigned int rshift_ulongy_into_uint_rounding(ulongy x)
{
        ulongy r = rshift_ulongy(x, FXP_INT_BITS_M1);
        unsigned int rbit = r.lo & 1u;
        r = rshift_ulongy(r, 1);
        r = ulongy_add_uint(r, rbit);
        return r.lo;
}

inline ulongy ulongy_bitwise_and(ulongy x, ulongy y)
{
        ulongy anded = { x.hi & y.hi, x.lo & y.lo };
        return anded;
}

inline ulongy ulongy_hi_bitwise_and(ulongy x, unsigned int y)
{
        ulongy anded = { x.hi & y, x.lo };
        return anded;
}

inline ulongy ulongy_hi_bitwise_or(ulongy x, unsigned int y)
{
        ulongy ored = { x.hi | y, x.lo };
        return ored;
}


/*
 * Distributive multiplication approach for two uints,
 * returning a ulongy.
 * With xa, xb, ya and yb as the inner hi and lo "words"
 * (e.g. half ints) of the unsigned int arguments, so
 *      x == (xa << NWORD_BITS) | xb
 *      y == (ya << NWORD_BITS) | yb
 * NWORD_BITS == (FXP_INT_BITS / 2), then:
 *
 *                        ya |    yb
 *                   *    xa |    xb
 *   ---------------------------------
 *   |          c.o. |  k0hi    k0lo | -> k0 = xb * yb
 *   |          k1hi |  k1lo         | -> k1 = xb * ya
 *   |          k2hi |  k2lo         | -> k2 = xa * yb
 *   |  k3hi    k3lo |               | -> k3 = xa * ya
 *   --------------------------------
 *   [    dmul_hi   ]-[    dmul_lo   ]
 *
 * k# are all intermediate uints resulting from multiplying
 * the words in x and y, as indicated on the right side above.
 * k#hi and k#lo are the high and low words in k#.
 * The result of the calculation will be two uints shown at the
 * bottom: dmul_hi and dmul_lo. Together they make up
 * the resulting pseudo-ulong.
 *
 * The operations needed to build those two final uints:
 * For dmul_lo:
 *    1. rsum = xbyb_hi + xbya_lo + xayb_lo
 *    2. dmul_lo = (rsum << NWORD_BITS) | xbyb_lo
 * For dmul_hi:
 *    3. carry_over = rsum >> NWORD_BITS
 *    4. dmul_hi = xaya + xayb_hi + xbya_hi + carry_over
 */
inline ulongy ulongy_from_dmul(unsigned int x, unsigned int y)
{
        unsigned int xa, xb, ya, yb, xbyb, xbya, xayb, xaya;
        unsigned int xbyb_hi, xbyb_lo, xbya_hi, xbya_lo, xayb_hi, xayb_lo;
        unsigned int rsum, carry_over, dproduct_hi, dproduct_lo;
        xa = x >> FXP_WORD_BITS;
        xb = (x & FXP_RWORD_MASK);
        ya = y >> FXP_WORD_BITS;
        yb = (y & FXP_RWORD_MASK);
        xbyb = xb * yb;   // k0
        xbya = xb * ya;   // k1
        xayb = xa * yb;   // k2
        xaya = xa * ya;   // k3
        xbyb_hi = xbyb >> FXP_WORD_BITS;
        xbyb_lo = xbyb & FXP_RWORD_MASK;
        xbya_hi = xbya >> FXP_WORD_BITS;
        xbya_lo = xbya & FXP_RWORD_MASK;
        xayb_hi = xayb >> FXP_WORD_BITS;
        xayb_lo = xayb & FXP_RWORD_MASK;
        rsum = xbyb_hi + xbya_lo + xayb_lo;
        carry_over = rsum >> FXP_WORD_BITS;
        dproduct_hi = xaya + xayb_hi + xbya_hi + carry_over;
        dproduct_lo = (rsum << FXP_WORD_BITS) | xbyb_lo;
        ulongy dproduct = { dproduct_hi, dproduct_lo };
        return dproduct;
}

inline unsigned int dmul_into_uint(unsigned int x, unsigned int y)
{
        return ulongy_hi_to_uint_rounded(ulongy_from_dmul(x, y));
}

/*
 * Equivalent to dmul_ulongs in fxp_l.c, but here using ulongys
 */
inline ulongy dmul_ulongys(ulongy x, ulongy y)
{
        ulongy xayb, yaxb, xaya, xbyb;
        ulongy qr1, qrsum, product;
        xayb = ulongy_from_dmul(x.hi, y.lo);
        yaxb = ulongy_from_dmul(y.hi, x.lo);
        xaya = ulongy_from_dmul(x.hi, y.hi);
        xbyb = ulongy_from_dmul(x.lo, y.lo);
        //qr1 = (xayb & FXP_RINT_MASK);
        qr1 = ulongy_create(0u, xayb.lo);
        //qr2 = (yaxb & FXP_RINT_MASK);
        //qr3 = (xbyb >> FXP_INT_BITS);
        unsigned int rbit1 = (xbyb.lo >> FXP_INT_BITS_M1);
        //qrsum = qr1 + qr2 + qr3 + rbit1;
        qrsum = ulongy_add_uint(
                    ulongy_add_uint(
                        ulongy_add_uint(qr1, yaxb.lo),  // qr1 + qr2
                        xbyb.hi),                       // + qr3
                    rbit1);                             // + rbit1
        //ql1 = (xayb >> FXP_INT_BITS);
        //ql2 = (yaxb >> FXP_INT_BITS);
        unsigned int rbit2 = (qrsum.lo >> FXP_INT_BITS_M1);
        //ql3 = (qrsum >> FXP_INT_BITS);
        //product = xaya + ql1 + ql2 + ql3 + rbit2;
        product = ulongy_add_uint(
                    ulongy_add_uint(
                        ulongy_add_uint(
                            ulongy_add_uint(xaya, xayb.hi), // xaya + ql1
                            yaxb.hi),                       // + ql2
                        qrsum.hi),                          // + ql3
                    rbit2);                                 // + rbit2
        return product;
}

// Same as above, but as if y.b was 0
inline ulongy dmul_ulongy_x_uint(ulongy x, unsigned int y)
{
        ulongy yaxa, yaxb, product;
        yaxa = ulongy_from_dmul(y, x.hi);
        yaxb = ulongy_from_dmul(y, x.lo);
        unsigned int rbit = (yaxb.lo >> FXP_INT_BITS_M1);
        product = ulongy_add_uint(
                        ulongy_add_uint(
                                yaxa, yaxb.hi), rbit);
        return product;
}

void print_ulongy_as_hex(ulongy x)
{
        printf("x%X%X", x.hi, x.lo);
}

void print_ulongy_as_bin(ulongy x)
{
        int i = FXP_INT_BITS;
        while (i > 0) {
                int bit = (x.hi >> (i - 1)) & 1u;
                printf("%d", bit);
                i--;
        }
        i = FXP_INT_BITS;
        while (i > 0) {
                int bit = (x.lo >> (i - 1)) & 1u;
                printf("%d", bit);
                i--;
        }
}
