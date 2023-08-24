#!/bin/bash
# Compiles and runs bkm.c

clear
if [ -f "b.out" ]; then
    rm b.out
fi
gcc bkm.c fxp.h fxp_aux.h fxp.c ulongy.h ulongy.c fxp_extern.h fxp_conv.h fxp_conv.c -lm -o b.out
if [ -f "b.out" ]; then
    ./b.out
fi
