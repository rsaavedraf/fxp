# fxp
Implementation of Fixed Point numbers in C using integers and their 12 least 
significant bits (configurable in fxp.h) for the fraction part.

For more information about Fixed Point arithmetic:
[Wikipedia: Fixed-point arithmetic](https://en.wikipedia.org/wiki/Fixed-point_arithmetic)

By Raul Saavedra

## Last updates:
- 2022-11-15: first version
- 2023-01-04: bug fix in fxp_sum
- 2023-01-05: multiplication now done not using longs, entirely with ints, so
supporting systems where sizeof(long) == sizeof(int)
- 2023-01-06: bugs fixed in new safe multiplication, also in the safe sum for a
border case.
- 2023-01-08: new version with runtime-modifiable number of bits to use for
the frac part. Fine-tuned all code, now all tests run with 0 warnings for
any chosen number of frac bits (tester runs with 0, 1, 12, 16, 24, and 31)

To Do: write alternative implementations for division also using just ints,
not longs.
