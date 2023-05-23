# fxp
Implementation of Fixed Point numbers in C using integers. Number of 
bits to use for fraction part configurable at runtime.

For more information about Fixed Point arithmetic:
[Wikipedia: Fixed-point arithmetic](https://en.wikipedia.org/wiki/Fixed-point_arithmetic)

The code implements overflow-safe operations for +, -, *, and /,
following CMU SEI's INT32-C recommendations. For further details:
[INT32-C. Ensure that operations on signed integers do not result in overflow](https://wiki.sei.cmu.edu/confluence/display/c/INT32-C.+Ensure+that+operations+on+signed+integers+do+not+result+in+overflow)

lg2(), pow2() implemented using the [BKM algorithm (Wikipedia)](https://en.wikipedia.org/wiki/BKM_algorithm)

Additional logarithm and power functions (ln, lg10, exp, pow10, sqrt, and 
powxy,) all implemented through lg2 and pow2. Behind the scenes, these 
functions work with longs or emulated longs, so as to avoid calculation 
inaccuracies that would inevitably appear if working just with the chosen 
number of frac bits, or just with int-size precision. Goals of this 
implementation are mostly flexibility (hence the configurable frac bits,) yet 
also ultimate precision if desired. All tests, including those for the power 
functions, should ideally run with no inaccuracy warnings, or as few as 
possible, regardless of frac bits in use.

Later:
- Trigonometric functions
- Saturated mode
- Rewrite in Rust maybe?

&nbsp;
## How To Try It
You only need a C compiler like gcc installed.
Simply clone this repo and then run the following script on its folder:

    ./r

That will compile with gcc locally, and will run the tester program. 
It should produce a long testing output, like the one in the 
*output.\*.txt* file(s). At the end it will ideally show a list with 
the observed number of warnings for each of the frac bit 
configurations tested. If an error larger than the maximum allowed 
is found, an assert will get triggered stopping the program 
immediately, showing the function call that was being tested right 
then.

Up to now I have run and tested it successfully (zero asserts)
on the CPUs listed below, using gcc v11.3.0 both without and with some 
-O options:

[Little Endian](https://en.wikipedia.org/wiki/Endianness) CPUs:
- Intel Core i7-6700K
- Arm Cortex-A72 (-> Raspberry Pi 4 Model B)

[Big Endian](https://en.wikipedia.org/wiki/Endianness) CPUs:
- (To be added when tried/receiving output from test runs)

&nbsp;
## Your Feedback/Help is Welcome!
If you try it on any CPU not listed above, specially one with an 
instruction set architecture different from the ones represented above,
please get in touch and send me your output. I'll be glad to credit you 
for it here if you'd like! :)
Either email me or send a PR for an output file named along the 
lines of *output.cpu_model.txt*

If you encounter any problems or errors, please don't hesitate to 
let me know as well, glad to also credit troubleshooting 
input/feedback. Just please mention compiler, compiler options used, 
and your hardware details.

To use these Fix Point Numbers yourself, for now you only need files 
*fxp.c* and its dependencies if strictly using ints and only ints. 
If allowing longs, then also file *fxp_l.c*. Most other files are 
auxiliary (to convert to and from floating points, to print out, 
test, etc.)
