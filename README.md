# fxp
Implementation of Fixed Point numbers in C using integers. Number of bits
to use for fraction part configurable at runtime.

For more information about Fixed Point arithmetic:
[Wikipedia: Fixed-point arithmetic](https://en.wikipedia.org/wiki/Fixed-point_arithmetic)

The code implements overflow-safe operations for +, -, *, and /,
following CMU SEI's INT32-C recommendations. For further details:
[INT32-C. Ensure that operations on signed integers do not result in overflow](https://wiki.sei.cmu.edu/confluence/display/c/INT32-C.+Ensure+that+operations+on+signed+integers+do+not+result+in+overflow)

lg2(), pow2() implemented using the [BKM algorithm (Wikipedia)](https://en.wikipedia.org/wiki/BKM_algorithm)

ln() and lg10() implemented multiplying lg2 by needed factor at full 
int-size precision (not just current fxp frac bits) to avoid precision loss.

exp(), pow(), and sqrt() coming soon.

Once completed with trigonometric functions, planning to rewrite in Rust.

By Raul Saavedra
