#!/bin/bash

clear
if [ -f "b.out" ]; then
    rm b.out
fi
gcc bkm.c fxp.h fxp_aux.h fxp.c fxp_extern.h fxp_conv.h fxp_conv.c -lm -o b.out
if [ -f "b.out" ]; then
    ./b.out
fi
