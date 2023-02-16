/*
 * fxp_constants.c
 * By Raul Saavedra, Bonn, Germany
 *
 * v1: 2023-02-05
 */

// Important constants as fxps

//            ----*----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8----*----9----*----0
// e:       2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274
// bin e:  10.1011011111100001010100010110001010001010111011010010101001101010101111110111000101011000100000001001
// e using 32 bits (1 sign + 2 whole + 29 frac bits), and 64 bits (1+2+61)
#define FXP_E_I32 0x56fc2a2c
#define FXP_E_I64 0x56fc2a2c515da54d
// Values of e at different precision levels:
// Bits Expansion
//         7    6    5    4    3    2    1    0
// 32   0101 0110 1111 1100 0010 1010 0010 1100
// 64   0101 0110 1111 1100 0010 1010 0010 1100 0101 0001 0101 1101 1010 0101 0100 1101
// For more details, or even more precise versions of e:
//  https://apod.nasa.gov/htmltest/gifcity/e.2mil
//  https://www.math.utah.edu/~pa/math/e.html
//  https://rosettacode.org/wiki/Calculating_the_value_of_e#C
//  https://www.exploringbinary.com/pi-and-e-in-binary/

//            ----*----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8----*----9----*----0
// pi:      3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679
// bin pi: 11.0010010000111111011010101000100010000101101000110000100011010011000100110001100110001010001011100000
#define FXP_PI_I32 0x6487ed51
#define FXP_PI_I64 0x6487ed5110b4611a
// Values of pi at different precision levels:
// Bits Expansion
//         7    6    5    4    3    2    1    0
// 32   0110 0100 1000 0111 1110 1101 0101 0001
// 64   0110 0100 1000 0111 1110 1101 0101 0001 0001 0000 1011 0100 0110 0001 0001 1010
// For more details, or even more precise versions of Pi:
//  https://oeis.org/A004601/constant
//  https://oeis.org/A000796/constant
//  http://www.befria.nu/elias/pi/binpi.html
//  https://www.exploringbinary.com/pi-and-e-in-binary/

// ln(2) = 0.693147180559945309417232121458176568075500134360255254120680009493393621969694715605863326996418687
#define FXP_LN_2_I32 0x58b90bfc
#define FXP_LN_2_I64 0x58b90bfbe8e7be39
// https://oeis.org/A002162

// lg10(2) = .301029995663981195213738894724493026768189881462108541310427461127108189274424509486927252118186172040684
#define FXP_LG10_2_I32 0x268826a1
#define FXP_LG10_2_I64 0x268826a13ef3fe80
//  https://oeis.org/A007524


// lg2(10) = 3.32192809488736234787031942948939017586483139302458061205475639581593477660862521585013974335937015
//#define FXP_LG2_10_I32 0x0
//#define FXP_LG2_10_I64 0x0
// Bits Expansion
//  https://oeis.org/A020862

// lg2(e) = 1.4426950408889634073599246810019
//#define FXP_LOG2E_I32 0x5c551d95
//#define FXP_LOG2E_I64 0x5c551d94ae0bf8cf
// Bits Expansion
// 32   0-101 1100 0101 0101 0001 1101 1001 0101  (<== last bit here rounded)
// 64   0-101 1100 0101 0101 0001 1101 1001 0100 1010 1110 0000 1011 1111 1000 1100 1111
// https://oeis.org/A007525/constant

// sqrt(2) = 1.414213562373095048801688724209698078569671875376948073176679737990732478462107038...
//#define FXP_SQRT2_I32 0x5a82799a
//#define FXP_SQRT2_I64 0x5a827999fcef3242
// In binary: 1.01101010000010011110011001100111111100111011110011001001000010001011001011111011...
//  https://community.wolfram.com/groups/-/m/t/1063480
//  https://oeis.org/A002193
//  https://oeis.org/A004539
