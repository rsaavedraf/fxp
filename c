#!/bin/bash

rm a.out
clear
gcc fxp.h fxp.c fxp_aux.h fxp_aux.c fxp_constants.h fxp_constants.c -lm
./a.out
