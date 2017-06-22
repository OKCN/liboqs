#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "local.h"

#include <oqs/aes.h>

#define min(x, y) (((x) < (y)) ? (x) : (y))

// Pack the input uint16 vector into a char output vector, copying lsb bits
// from each input element. If inlen * lsb / 8 > outlen, only outlen * 8 bits
// are copied.
void oqs_kex_lwe_okcn_pack(unsigned char *out, const size_t outlen, const uint16_t *in, const size_t inlen, const unsigned char lsb) {
	memset(out, 0, outlen);

	size_t i = 0;            // whole bytes already filled in
	size_t j = 0;            // whole uint16_t already copied
	uint16_t w = 0;          // the leftover, not yet copied
	unsigned char bits = 0;  // the number of lsb in w
	while (i < outlen && (j < inlen || ((j == inlen) && (bits > 0)))) {
		/*
		in: |        |        |********|********|
		                      ^
		                      j
		w : |   ****|
		        ^
		       bits
		out:|**|**|**|**|**|**|**|**|* |
		                            ^^
		                            ib
		*/
		unsigned char b = 0;  // bits in out[i] already filled in
		while (b < 8) {
			int nbits = min(8 - b, bits);
			uint16_t mask = (1 << nbits) - 1;
			unsigned char t = (w >> (bits - nbits)) & mask;  // the bits to copy from w to out
			out[i] = out[i] + (t << (8 - b - nbits));
			b += nbits;
			bits -= nbits;
			w &= ~(mask << bits);  // not strictly necessary; mostly for debugging

			if (bits == 0) {
				if (j < inlen) {
					w = in[j];
					bits = lsb;
					j++;
				} else {
					break;  // the input vector is exhausted
				}
			}
		}
		if (b == 8) {  // out[i] is filled in
			i++;
		}
	}
}

// Unpack the input char vector into a uint16_t output vector, copying lsb bits
// for each output element from input. outlen must be at least ceil(inlen * 8 /
// lsb).
void oqs_kex_lwe_okcn_unpack(uint16_t *out, const size_t outlen, const unsigned char *in, const size_t inlen, const unsigned char lsb) {
	memset(out, 0, outlen * sizeof(uint16_t));

	size_t i = 0;            // whole uint16_t already filled in
	size_t j = 0;            // whole bytes already copied
	unsigned char w = 0;     // the leftover, not yet copied
	unsigned char bits = 0;  // the number of lsb bits of w
	while (i < outlen && (j < inlen || ((j == inlen) && (bits > 0)))) {
		/*
		in: |  |  |  |  |  |  |**|**|...
		                      ^
		                      j
		w : | *|
		      ^
		      bits
		out:|   *****|   *****|   ***  |        |...
		                      ^   ^
		                      i   b
		*/
		unsigned char b = 0;  // bits in out[i] already filled in
		while (b < lsb) {
			int nbits = min(lsb - b, bits);
			uint16_t mask = (1 << nbits) - 1;
			unsigned char t = (w >> (bits - nbits)) & mask;  // the bits to copy from w to out
			out[i] = out[i] + (t << (lsb - b - nbits));
			b += nbits;
			bits -= nbits;
			w &= ~(mask << bits);  // not strictly necessary; mostly for debugging

			if (bits == 0) {
				if (j < inlen) {
					w = in[j];
					bits = 8;
					j++;
				} else {
					break;  // the input vector is exhausted
				}
			}
		}
		if (b == lsb) {  // out[i] is filled in
			i++;
		}
	}
}

#define LWE_DIV_ROUNDUP(x, y) (((x) + (y)-1) / y)

// define parameters for "recommended" parameter set
#define PARAMS_N 712
#define PARAMS_NBAR 8
#define PARAMS_LOG2Q 14
#define PARAMS_Q (1 << PARAMS_LOG2Q)
#define PARAMS_EXTRACTED_BITS 4
#define PARAMS_KEY_BITS 256
#define PARAMS_STRIPE_STEP 8
#define PARAMS_SINGLE_HINT_BITS 8
#define PARAMS_REC_HINT_LENGTH LWE_DIV_ROUNDUP(PARAMS_NBAR * PARAMS_NBAR * PARAMS_SINGLE_HINT_BITS, 8)
// pre-process code to obtain "recommended" functions
#define MACRIFY(NAME) NAME ## _recommended
#include "lwe_macrify.c"
// undefine macros to avoid any confusion later
#undef PARAMS_SINGLE_HINT_BITS
#undef PARAMS_N
#undef PARAMS_NBAR
#undef PARAMS_LOG2Q
#undef PARAMS_Q
#undef PARAMS_EXTRACTED_BITS
#undef PARAMS_KEY_BITS
#undef PARAMS_STRIPE_STEP
#undef PARAMS_REC_HINT_LENGTH
#undef MACRIFY
