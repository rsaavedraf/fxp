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

exp() and pow10() now implemented (using longs, int-only versions + powxy() and sqrt() coming very soon.)

Once completed with trigonometric functions, planning to rewrite in Rust.

## To try it
You only need a C compiler like gcc installed. Simply clone this repo and then run the following script on its folder:

    ./r

That will compile locally and run the tester program. It should produce a long testing output like the one in output.txt.

Up to now run and tested on the following CPUs using gcc v11.3.0 with no -O options:

[Little Endian](https://en.wikipedia.org/wiki/Endianness) CPUs:
- Intel Core i7-6700K
- Arm Cortex-A72 (-> Raspberry Pi 4 Model B)

[Big Endian](https://en.wikipedia.org/wiki/Endianness) CPUs:
- (To be added when tried)

To use the Fix Point Numbers yourself, you only need files *fxp.h* and *fxp.c* when strictly
using ints and only ints. If also using longs, then also files *fxp_l.h* and *fxp_l.c*. All the other files
are auxiliary (to convert to and from floating points, print out, test, etc.)

If you encounter any problems or errors please don't hesitate to let me know.
