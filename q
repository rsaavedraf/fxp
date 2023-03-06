#!/bin/bash

rm q.out
clear
#gcc test_quad.c -march=skylake -lquadmath
gcc test_quad.c -lquadmath -o q.out
./q.out

