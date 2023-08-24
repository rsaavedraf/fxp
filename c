#!/bin/bash
# Compiles and runs fxp_tconst.c

clear
if [ -f "c.out" ]; then
    rm c.out
fi
gcc fxp.h fxp.c ulongy.h ulongy.c fxp_extern.h fxp_aux.h fxp_aux.c print_as_bits.h print_as_bits.c fxp_conv.h fxp_conv.c fxp_tconst.c -lm -o c.out
if [ -f "c.out" ]; then
    ./c.out
fi
