#!/bin/bash
# Compiles and runs dmul_tester.c

clear
if [ -f "d.out" ]; then
    rm d.out
fi
gcc fxp.h fxp_extern.h ulongy.h ulongy.c fxp_l.h fxp_l.c fxp_aux.h fxp_aux.c print_as_bits.h print_as_bits.c fxp.c fxp_conv.h fxp_conv.c dmul_tester.c -lm -o d.out
if [ -f "d.out" ]; then
    ./d.out
fi
