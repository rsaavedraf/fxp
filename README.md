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
supporting systems where sizeof(long) == sizeof(int)
- 2023-01-06: bugs fixed in new safe multiplication, also in the safe sum for a
border case.
- 2023-01-08: new version with runtime-modifiable number of bits to use for
the frac part. Fine-tuned all code, now all tests run with 0 warnings for
any chosen number of frac bits. Tester checks different nums of frac bits.
- 2023-01-09: Default number of frac bits is set to 14, which the tester
runs last. Beautified tester output organization.
- 2023-01-11: Using gcc's own built in to count leading zeros.
Tester program now using long-doubles to better compare
different implementations of the operations. Thanks to that, a few bugs
exposed and fixed along the way.

To Do: write alternative implementations for division also using just ints,
not longs.
