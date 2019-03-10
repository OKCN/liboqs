#include "params.h"
#include <oqs/rand.h>
#include <oqs/sha3.h>

 #define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

typedef struct {
	uint16_t coeffs[PARAM_N];
#if defined(WINDOWS)
} poly;
#else
} poly __attribute__((aligned(32)));
#endif

static const uint32_t qinv = 12287; // -inverse_mod(p,2^18)
static const uint32_t rlog = 18;

static uint16_t montgomery_reduce(uint32_t a) {
	uint32_t u;

	u = (a * qinv);
	u &= ((1 << rlog) - 1);
	u *= PARAM_Q;
	a = a + u;
	return a >> 18;
}

static uint16_t barrett_reduce(uint16_t a) {
	uint32_t u;

	u = ((uint32_t) a * 5) >> 16;
	u *= PARAM_Q;
	a -= u;
	return a;
}

static void bitrev_vector(uint16_t *poly) {
	unsigned int i, r;
	uint16_t tmp;

	for (i = 0; i < PARAM_N; i++) {
		r = E8_bitrev_table[i];
		if (i < r) {
			tmp = poly[i];
			poly[i] = poly[r];
			poly[r] = tmp;
		}
	}
}

static void mul_coefficients(uint16_t *poly, const uint16_t *factors) {
	unsigned int i;

	for (i = 0; i < PARAM_N; i++) {
		poly[i] = montgomery_reduce((poly[i] * factors[i]));
	}
}

/* GS_bo_to_no; omegas need to be in Montgomery domain */
static void ntt(uint16_t *a, const uint16_t *omega) {
	int i, start, j, jTwiddle, distance;
	uint16_t temp, W;

	for (i = 0; i < 10; i += 2) {
		// Even level
		distance = (1 << i);
		for (start = 0; start < distance; start++) {
			jTwiddle = 0;
			for (j = start; j < PARAM_N - 1; j += 2 * distance) {
				W = omega[jTwiddle++];
				temp = a[j];
				a[j] = (temp + a[j + distance]); // Omit reduction (be lazy)
				a[j + distance] = montgomery_reduce((W * ((uint32_t) temp + 3 * PARAM_Q - a[j + distance])));
			}
		}

		// Odd level
		distance <<= 1;
		for (start = 0; start < distance; start++) {
			jTwiddle = 0;
			for (j = start; j < PARAM_N - 1; j += 2 * distance) {
				W = omega[jTwiddle++];
				temp = a[j];
				a[j] = barrett_reduce((temp + a[j + distance]));
				a[j + distance] = montgomery_reduce((W * ((uint32_t) temp + 3 * PARAM_Q - a[j + distance])));
			}
		}
	}
}

// Pack the input uint16 vector into a char output vector, copying lsb bits
// from each input element. If inlen * lsb / 8 > outlen, only outlen * 8 bits
// are copied.
static void inline frodo_pack(unsigned char *out, const size_t outlen, const uint16_t *in, const size_t inlen, const unsigned char lsb) {
	memset(out, 0, outlen);

	size_t i = 0;           // whole bytes already filled in
	size_t j = 0;           // whole uint16_t already copied
	uint16_t w = 0;         // the leftover, not yet copied
	unsigned char bits = 0; // the number of lsb in w
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
		unsigned char b = 0; // bits in out[i] already filled in
		while (b < 8) {
			int nbits = min(8 - b, bits);
			uint16_t mask = (1 << nbits) - 1;
			unsigned char t = (w >> (bits - nbits)) & mask; // the bits to copy from w to out
			out[i] += t << (8 - b - nbits);
			b += nbits;
			bits -= nbits;
			w &= ~(mask << bits); // not strictly necessary; mostly for debugging

			if (bits == 0) {
				if (j < inlen) {
					w = in[j];
					bits = lsb;
					j++;
				} else {
					break; // the input vector is exhausted
				}
			}
		}
		if (b == 8) { // out[i] is filled in
			i++;
		}
	}
}

// Unpack the input char vector into a uint16_t output vector, copying lsb bits
// for each output element from input. outlen must be at least ceil(inlen * 8 /
// lsb).
static inline void frodo_unpack(uint16_t *out, const size_t outlen, const unsigned char *in, const size_t inlen, const unsigned char lsb) {
	memset(out, 0, outlen * sizeof(uint16_t));

	size_t i = 0;           // whole uint16_t already filled in
	size_t j = 0;           // whole bytes already copied
	unsigned char w = 0;    // the leftover, not yet copied
	unsigned char bits = 0; // the number of lsb bits of w
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
		unsigned char b = 0; // bits in out[i] already filled in
		while (b < lsb) {
			int nbits = min(lsb - b, bits);
			uint16_t mask = (1 << nbits) - 1;
			unsigned char t = (w >> (bits - nbits)) & mask; // the bits to copy from w to out
			out[i] += t << (lsb - b - nbits);
			b += nbits;
			bits -= nbits;
			w &= ~(mask << bits); // not strictly necessary; mostly for debugging

			if (bits == 0) {
				if (j < inlen) {
					w = in[j];
					bits = 8;
					j++;
				} else {
					break; // the input vector is exhausted
				}
			}
		}
		if (b == lsb) { // out[i] is filled in
			i++;
		}
	}
}

static void poly_frombytes(poly *r, const unsigned char *a, const int lenPerCoeff) {
    frodo_unpack(r->coeffs, PARAM_N, a, POLY_BYTES, lenPerCoeff);
}

static inline void poly_reduce(poly *p) {
    int i, t, m, c;

    for (i = 0;i < PARAM_N;i++) {
        t = barrett_reduce(p->coeffs[i]); //Make sure that coefficients have only 14 bits
		m = t - PARAM_Q;
		c = m;
		c >>= 15;
		t = m ^ ((t ^ m) & c); // <Make sure that coefficients are in [0,q]
        p->coeffs[i] = t;
    }
}

static inline void poly_cutTail(poly *p, const int cutBitNum) {
    int i;
    
    for (i = 0;i < PARAM_N;i++) {
        p->coeffs[i] >>= cutBitNum;
    }
}

static inline void poly_extTail(poly *p, const int cutBitNum) {
    int i;

    for (i = 0;i < PARAM_N;i++) {
        p->coeffs[i] = (p->coeffs[i] << 1) | 0x1;
        p->coeffs[i] <<= cutBitNum - 1;
    }
}

static void poly_tobytes(unsigned char *r, const poly *p, const int lenPerCoeff) {
    frodo_pack(r, POLY_BYTES, p->coeffs, PARAM_N, lenPerCoeff);
}

static void poly_uniform(poly *a, const unsigned char *seed) {
	unsigned int pos = 0, ctr = 0;
	uint16_t val;
	uint64_t state[OQS_SHA3_STATESIZE];
	unsigned int nblocks = 8;
	uint8_t buf[OQS_SHA3_SHAKE128_RATE * nblocks];

	OQS_SHA3_shake128_absorb(state, seed, E8_SEEDBYTES);

	OQS_SHA3_shake128_squeezeblocks((unsigned char *) buf, nblocks, state);

	while (ctr < PARAM_N) {
		val = (buf[pos] | ((uint16_t) buf[pos + 1] << 8)) & 0x3fff; // Specialized for q = 12889
		if (val < PARAM_Q) {
			a->coeffs[ctr++] = val;
		}
		pos += 2;
		if (pos > OQS_SHA3_SHAKE128_RATE * nblocks - 2) {
			nblocks = 1;
			OQS_SHA3_shake128_squeezeblocks((unsigned char *) buf, nblocks, state);
			pos = 0;
		}
	}
}

#define MASK_32_1  0x55555555
#define MASK_32_2  0x33333333
#define MASK_32_4  0x0f0f0f0f
#define MASK_32_8  0x00ff00ff
#define MASK_32_16 0x0000ffff

static uint16_t popcnt32(uint32_t x)
{
    x = (x & MASK_32_1 ) + ((x >>  1) & MASK_32_1 ); 
    x = (x & MASK_32_2 ) + ((x >>  2) & MASK_32_2 );
    x = (x & MASK_32_4 ) + ((x >>  4) & MASK_32_4 );
    x = (x & MASK_32_8 ) + ((x >>  8) & MASK_32_8 );
    x = (x & MASK_32_16) + ((x >> 16) & MASK_32_16);
    return x;
}

#define MASK_16_1  0x5555
#define MASK_16_2  0x3333
#define MASK_16_4  0x0f0f
#define MASK_16_8  0x00ff

static uint16_t popcnt16(uint16_t x)
{
    x = (x & MASK_16_1 ) + ((x >>  1) & MASK_16_1 ); 
    x = (x & MASK_16_2 ) + ((x >>  2) & MASK_16_2 );
    x = (x & MASK_16_4 ) + ((x >>  4) & MASK_16_4 );
    x = (x & MASK_16_8 ) + ((x >>  8) & MASK_16_8 );
    return x;
}

#define MASK_8_1  0x55
#define MASK_8_2  0x33
#define MASK_8_4  0x0f

static uint8_t popcnt8(uint8_t x)
{
    x = (x & MASK_8_1 ) + ((x >>  1) & MASK_8_1 ); 
    x = (x & MASK_8_2 ) + ((x >>  2) & MASK_8_2 );
    x = (x & MASK_8_4 ) + ((x >>  4) & MASK_8_4 );
    return x;
}

static void poly_getnoise(poly *r, OQS_RAND *rand)
{
    unsigned char buf[3*PARAM_N];
    uint16_t *tp_a = (uint32_t *) buf;
    uint16_t ones, twos;
    uint8_t *tp_b = buf + PARAM_N*2;
    int i;

	rand->rand_n(rand, buf, 3*PARAM_N);

    for(i=0;i<PARAM_N;i++)
    {
        ones = popcnt16(tp_a[i]); // 16-bits a = 16
        twos = popcnt8(tp_b[i]);  // 8-bits  b = 8
        r->coeffs[i] = (twos<<1) + ones - 16 + PARAM_Q;
    }
}

static void poly_pointwise(poly *r, const poly *a, const poly *b) {
	int i;
	uint16_t t;
	for (i = 0; i < PARAM_N; i++) {
		t = montgomery_reduce(3186 * b->coeffs[i]);         /* t is now in Montgomery domain */
		r->coeffs[i] = montgomery_reduce(a->coeffs[i] * t); /* r->coeffs[i] is back in normal domain */
	}
}

static void poly_add(poly *r, const poly *a, const poly *b) {
	int i;
	for (i = 0; i < PARAM_N; i++) {
		r->coeffs[i] = barrett_reduce(a->coeffs[i] + b->coeffs[i]);
	}
}

static void poly_ntt(poly *r) {
	bitrev_vector(r->coeffs);
	mul_coefficients(r->coeffs, E8_psis_bitrev_montgomery);
	ntt((uint16_t *) r->coeffs, E8_omegas_montgomery);
}

static void poly_invntt(poly *r) {
	bitrev_vector(r->coeffs);
	ntt((uint16_t *) r->coeffs, E8_omegas_inv_montgomery);
	mul_coefficients(r->coeffs, E8_psis_inv_montgomery);
}

//Error Correction:

static uint8_t encode_hamming(uint8_t m)
{
    uint8_t a0 = m & 1;
    uint8_t a1 = (m >> 1) & 1;
    uint8_t a2 = (m >> 2) & 1;
    uint8_t a3 = (m >> 3) & 1;
    return ((-a1) & 0x0F) ^ ((-a2) & 0x3c) ^ ((-a3) & 0xF0) ^ ((-a0) & 0x55);
}

static int32_t sqr(int32_t x) {
    return x*x;
}

static uint32_t abs_q(uint32_t x) {
    uint32_t mask = -(((PARAM_Q - (x<<1))>>31) & 1);
    return (mask & (PARAM_Q - x)) | ((~mask) & x); 
}

static uint8_t decode_D8_00(uint32_t *cost, uint32_t tmp_cost[8][2])
{
    uint8_t m[2] = {0, 0}, xor_sum = 0, res;
    uint32_t min_diff = ~0U>>2, r;
    uint32_t min_i = 0;
    int i;

    *cost = 0;
    for (i = 0;i < 4;i++) {
        uint32_t c[2];
        c[0] = tmp_cost[i<<1][0] + tmp_cost[i<<1|1][0];
        c[1] = tmp_cost[i<<1][1] + tmp_cost[i<<1|1][1];
        r = ((c[1] - c[0]) >> 31) & 1;
        m[0] |= r << i;
        xor_sum ^= r;
        *cost += c[r];

        uint32_t diff = c[r^1] - c[r];
        r = ((diff - min_diff) >> 31) & 1;
        min_diff = ((-r) & diff) | ((-(r^1)) & min_diff);
        min_i = ((-r) & i) | ((-(r^1)) & min_i);
    }
    m[1] = m[0] ^ (1<<min_i);
    res = m[xor_sum];
    *cost += (-xor_sum) & min_diff;
    return res;
}

static uint8_t decode_D8_10(uint32_t *cost, uint32_t tmp_cost[8][2])
{
    uint8_t m[2] = {0, 0}, xor_sum = 0, res;
    uint32_t min_diff = ~0U>>2, r;
    uint32_t min_i = 0;
    int i;

    *cost = 0;
    for (i = 0;i < 4;i++) {
        uint32_t c[2];
        c[0] = tmp_cost[i<<1][1] + tmp_cost[i<<1|1][0];
        c[1] = tmp_cost[i<<1][0] + tmp_cost[i<<1|1][1];
        r = ((c[1] - c[0]) >> 31) & 1;
        m[0] |= r << i;
        xor_sum ^= r;
        *cost += c[r];

        uint32_t diff = c[r^1] - c[r];
        r = ((diff - min_diff) >> 31) & 1;
        min_diff = ((-r) & diff) | ((-(r^1)) & min_diff);
        min_i = ((-r) & i) | ((-(r^1)) & min_i);
    }
    m[1] = m[0] ^ (1<<min_i);
    res = m[xor_sum];
    *cost += (-xor_sum) & min_diff;
    return res;
}

static uint32_t const_abs(int32_t x) {
    uint32_t mask = x >> 31;
    return (x ^ mask) - mask;
}

static uint8_t decode_E8(uint32_t vec[8])
{
    uint32_t cost[2], r;
    uint8_t m[2], res;
    uint32_t tmp_cost[8][2];
    int i;

    for (i = 0;i < 8;i++) {
        tmp_cost[i][0] = sqr(abs_q(vec[i]));
        tmp_cost[i][1] = sqr(const_abs(vec[i] - (PARAM_Q>>1)));
    }

    m[0] = decode_D8_00(cost+0, tmp_cost);
    m[1] = decode_D8_10(cost+1, tmp_cost);
    r = ((cost[1] - cost[0]) >> 31) & 1;

    res = m[r];
    res = ((((res ^ (res << 1)) & 0x3) | ((res >> 1) & 4)) << 1) | r;

    return res;
}

static void con(poly *v, unsigned char *key, const poly *sigma, OQS_RAND *oqs_rand)
{
    int i;
    uint16_t rbit;
    uint8_t exkey[128];

	oqs_rand->rand_n(oqs_rand, key, 64);

    for (i = 0;i < 128;i++) {
        uint8_t m = (key[i>>1] >> ((i&1)<<2)) & 0xF;
        exkey[i] = encode_hamming(m);
    }
    for (i = 0;i < PARAM_N;i++) {
        rbit = (exkey[i>>3] >> (i & 0x7)) & 1;
        v->coeffs[i] = (((((uint32_t)(sigma->coeffs[i] + ((PARAM_Q>>1) & (-rbit)))) << 5) + (PARAM_Q>>1)) / PARAM_Q) & 0x1F;
    }
}

static void rec(unsigned char *key, const poly *v, const poly *sigma)
{
    int i, j, k;

    memset(key, 0, sizeof(unsigned char)*64);
    for (i = 0;i < 128;i++) {
        uint32_t tmp_v[8];
        for (j = 0;j < 8;j++) {
            k = (i<<3) | j;
            tmp_v[j] = ((((((uint32_t)v->coeffs[k])*PARAM_Q) + 16) >> 5) - sigma->coeffs[k] + PARAM_Q) % PARAM_Q;
        }
        key[i>>1] |= decode_E8(tmp_v) << ((i & 1) << 2);
    }
}
