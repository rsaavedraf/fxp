# fxp
Implementation of Fixed Point numbers in C using integers. Number of bits
to use for fraction part now configurable at runtime.

For more information about Fixed Point arithmetic:
[Wikipedia: Fixed-point arithmetic](https://en.wikipedia.org/wiki/Fixed-point_arithmetic)

The code implements overflow-safe operations for +, -, *, and /,
following CMU SEI's INT32-C recommendations. For further details:
[INT32-C. Ensure that operations on signed integers do not result in overflow](https://wiki.sei.cmu.edu/confluence/display/c/INT32-C.+Ensure+that+operations+on+signed+integers+do+not+result+in+overflow)


By Raul Saavedra

## Last updates:
- 2022-11-15: first version
- 2023-01-04: bug fix in fxp_sum
- 2023-01-05: multiplication now done not using longs, entirely with ints, so
supporting systems where sizeof(long) == sizeof(int), but for now at the
cost of performance, and precision loss of up to ~half the frac bits in use.
- 2023-01-06: bugs fixed in new safe multiplication, also in the safe sum for a
border case.
- 2023-01-08: new version with runtime-modifiable number of bits to use for
the frac part. Fine-tuned all code, now all tests run with 0 warnings for
any chosen number of frac bits. Tester checks different nums of frac bits.
- 2023-01-09: Default number of frac bits is set to 14.
Beautified tester output organization.
- 2023-01-11: Using gcc's own built in to count leading zeros.
Tester program now using long-doubles to better compare
different implementations of the operations. Thanks to that, a few bugs
exposed and fixed along the way.
- 2023-01-15: implemented fxp division using only ints (fxp_div)
Basically a software-based implementation of binary division
(tailored to fxp's.) Because of this, it is significantly slower:
average execution time of first version is ~6x that of fxp_div_l.
However it will work for systems in which sizeof(long) is not larger
than sizeof(int).
- 2023-01-29: improved testing framework using long doubles, and
refactored some utility functions related to testing and tracing
in a separate file: fxp_aux.c.
- 2023-01-31: significant precision improvement in fxp_mul(),
now offering pretty much the same precision as fxp_mul_l() while
still only using ints, no longs. The calculation uses the
distributive approach twice, now the second time within the frac
part only.
- 2023-02-01: bug fixed in fxp_xtimes.c. After this fix,
the true relative execution time of fxp_mul() is shown to be
effectively ~6x that of fxp_mul_l(). This makes a lot of sense,
as a similar ratio stands between fxp_div() and fxp_div_l().
Both * and / operations implemented using longs are then
significantly faster (~6x) than the int-only implementations.
