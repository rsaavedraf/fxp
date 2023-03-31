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
Goals of this implementation are mostly flexibility (hence the configurable frac bits,)
yet ultimate precision within the chosen number of frac bits.

exp() and pow10() now implemented (using longs, int-only versions + powxy() and sqrt() coming very soon.)

Once completed with trigonometric functions, planning to rewrite in Rust.

## How To Try It
You only need a C compiler like gcc installed.
Simply clone this repo and then run the following script on its folder:

    ./r

That will compile locally and run the tester program.
It should produce a long testing output, like the ones in the 
*output.\*.txt* files.
At the end it will ideally show a total of zero warnings.
If an error larger than the warning tolerance is found, then an assert 
will get triggered stopping the program immediately, showing the
function call that was being tested right then.

Up to now I have run and tested it successfully with zero warnings on 
a couple of CPUs listed below, using gcc v11.3.0 both without and 
with some -O options.

[Little Endian](https://en.wikipedia.org/wiki/Endianness) CPUs:
- Intel Core i7-6700K
- Arm Cortex-A72 (-> Raspberry Pi 4 Model B)

[Big Endian](https://en.wikipedia.org/wiki/Endianness) CPUs:
(To be added when tried/receiving output from test runs)

## Your Feedback is Welcome!
If you try it on any CPU not listed above, specially one with an 
instruction set architecture different from the ones represented above,
please get in touch and send me your output. I'll be glad to credit you 
for it here if you'd like! :)
Either email me or send a PR for an output file named along the 
lines of *output.cpu_model.txt*

If you encounter any problems or errors, please don't hesitate to 
let me know as well, glad to also credit troubleshooting input/feedback.
Just please mention compiler, compiler options used, and your hardware details.

To use these Fix Point Numbers yourself, for now you only need 
files *fxp.h* and *fxp.c* if strictly using ints and only ints.
If also using longs, then also files *fxp_l.h* and *fxp_l.c*. 
All the other files are auxiliary (to convert to and from floating points, 
print out, test, etc.)
