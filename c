#!/bin/bash

clear
if [ -f "c.out" ]; then
    rm c.out
fi
gcc fxp.h fxp.c fxp_extern.h fxp_aux.h fxp_aux.c fxp_conv.h fxp_conv.c fxp_tconst.c -lm -o c.out
if [ -f "c.out" ]; then
    ./c.out
fi
