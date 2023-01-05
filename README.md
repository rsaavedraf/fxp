# fxp
Implementation of Fixed Point numbers in C using integers and their 12 least 
significant bits (configurable in fxp.h) for the fraction part.

For more information about Fixed Point arithmetic:
[Wikipedia: Fixed-point arithmetic](https://en.wikipedia.org/wiki/Fixed-point_arithmetic)

By Raul Saavedra
Nov 15, 2022, Bonn

## Last updates:
- 2023-01-04: bug fix in fxp_sum
- 2023-01-05: multiplication now done not using longs, entirely with ints, so
supporting systems where sizeof(long) == sizeof(int)

To Do: write alternative implementations for division also entirely with ints.
