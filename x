#!/bin/bash

clear
if [ -f "x.out" ]; then
    rm x.out
fi
#gcc -O2 fxp.h fxp.c fxp_aux.h fxp_aux.c print_as_bits.h print_as_bits.c fxp_conv.h fxp_conv.c fxp_l.h fxp_l.c fxp_xtimes.c -lm -o x.out
gcc fxp.h fxp.c fxp_aux.h fxp_aux.c print_as_bits.h print_as_bits.c fxp_conv.h fxp_conv.c fxp_l.h fxp_l.c fxp_xtimes.c -lm -o x.out
if [ -f "x.out" ]; then
    ./x.out
fi
