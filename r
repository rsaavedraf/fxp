#!/bin/bash

clear
if [ -f "r.out" ]; then
    rm r.out
fi
gcc fxp.h fxp_extern.h fxp_l.h fxp_l.c fxp_aux.h fxp_aux.c print_as_bits.h print_as_bits.c fxp.c fxp_conv.h fxp_conv.c fxp_tester.c -lm -o r.out
if [ -f "r.out" ]; then
    ./r.out
fi
