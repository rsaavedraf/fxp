/*
 * fxp_constants.c
 * By Raul Saavedra, Bonn, Germany
 *
 * v1: 2023-02-05
 */

// Important constants as fxps

// pi = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803...
// Pi using 32 bits (1 sign + 2 whole + 29 frac bits), and 64 bits (1+2+61)
#define FXP_PI_I32 0x6487ed51
#define FXP_PI_I64 0x6487ed5110b4611a
// Rounded values of pi at different precision levels:
// Bits Expansion
// 32   0-110 0100 1000 0111 1110 1101 0101 0001
// 64   0-110 0100 1000 0111 1110 1101 0101 0001 0001 0000 1011 0100 0110 0001 0001 1010
// For more details, or even more precise versions of Pi:
//  https://www.exploringbinary.com/pi-and-e-in-binary/
//  https://oeis.org/A004601/constant
//  http://www.befria.nu/elias/pi/binpi.html

// e = 2.718281828459045235360287471352662497757247093699959574966967627724076630353547594571382...
#define FXP_E_I32 0x56fc2a2c
#define FXP_E_I64 0x56fc2a2c515da54d
// Rounded values of e at different precision levels:
// Bits Expansion
// 32   0-101 0110 1111 1100 0010 1010 0010 1100
// 64   0-101 0110 1111 1100 0010 1010 0010 1100 0101 0001 0101 1101 1010 0101 0100 1101
// For more details, or even more precise versions of e:
//  https://www.exploringbinary.com/pi-and-e-in-binary/
//  https://rosettacode.org/wiki/Calculating_the_value_of_e#C

// ln(2) = 0.69314718055994530941723212145818
#define FXP_LN2_I32 0x58b90bfc
#define FXP_LN2_I64 0x58b90bfbe8e7be3a
// Bits Expansion
// 32   0-101 1000 1011 1001 0000 1011 1111 1011  (<== last bit here rounded)
// 64   0-101 1000 1011 1001 0000 1011 1111 1010 1110 1000 1110 0111 1011 1110 0011 1010

// log2(e) = 1.4426950408889634073599246810019
#define FXP_LOG2E_I32 0x5c551d95
#define FXP_LOG2E_I64 0x5c551d94ae0bf8cf
// Bits Expansion
// 32   0-101 1100 0101 0101 0001 1101 1001 0101  (<== last bit here rounded)
// 64   0-101 1100 0101 0101 0001 1101 1001 0100 1010 1110 0000 1011 1111 1000 1100 1111

// sqrt(2) = 1.414213562373095048801688724209698078569671875376948073176679737990732478462107038...
#define FXP_SQRT2_I32 0x5a82799a
#define FXP_SQRT2_I64 0x5a827999fcef3242
// In binary: 1.01101010000010011110011001100111111100111011110011001001000010001011001011111011...
// For details on sqrt calculation of binary numbers:
//  https://community.wolfram.com/groups/-/m/t/1063480
//  https://oeis.org/A002193
//  https://oeis.org/A004539
