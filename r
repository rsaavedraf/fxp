#!/bin/bash

clear
if [ -f "r.out" ]; then
    rm r.out
fi
gcc fxp.h fxp.c fxp_l.h fxp_l.c fxp_aux.h fxp_aux.c fxp_conv.h fxp_conv.c fxp_tester.c -lm -o r.out
if [ -f "r.out" ]; then
    ./r.out
fi
