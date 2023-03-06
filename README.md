# fxp
Implementation of Fixed Point numbers in C using integers. Number of bits
to use for fraction part configurable at runtime.

For more information about Fixed Point arithmetic:
[Wikipedia: Fixed-point arithmetic](https://en.wikipedia.org/wiki/Fixed-point_arithmetic)

The code implements overflow-safe operations for +, -, *, and /,
following CMU SEI's INT32-C recommendations. For further details:
[INT32-C. Ensure that operations on signed integers do not result in overflow](https://wiki.sei.cmu.edu/confluence/display/c/INT32-C.+Ensure+that+operations+on+signed+integers+do+not+result+in+overflow)

lg2() now also implemented (ln() and lg10() as well by just calling lg2 and converting, but
being reworked to multiply by needed factor at full precision to avoid precision loss)

exp(), pow(), and sqrt() coming soon.

Once completed with trigonometric functions, planning to rewrite in Rust.

By Raul Saavedra
