#!/bin/bash

rm a.out
clear
gcc fxp.h fxp.c fxp_aux.h fxp_aux.c fxp_tester.c -lm
./a.out
