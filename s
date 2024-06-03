#!/bin/bash
# Compiles and runs sqrt.c

clear
if [ -f "s.out" ]; then
    rm s.out
fi
gcc print_as_bits.h sqrt.c fxp.h fxp_aux.h fxp.c ulongy.h ulongy.c fxp_extern.h fxp_conv.h fxp_conv.c print_as_bits.c -lm -o s.out
if [ -f "s.out" ]; then
    ./s.out
fi
