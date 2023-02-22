#!/bin/bash

clear
if [ -f "b.out" ]; then
    rm b.out
fi
gcc bkm.c -lm -o b.out
if [ -f "b.out" ]; then
    ./b.out
fi
