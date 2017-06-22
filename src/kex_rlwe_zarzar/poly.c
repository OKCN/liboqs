#include "params.h"
#include <oqs/rand.h>
#include <oqs/sha3.h>

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
		r = zarzar_bitrev_table[i];
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

	for (i = 0; i < 8; i += 2) {
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
    // last level
    distance <<= 1;
    for(start = 0; start < distance;start++)
    {
        jTwiddle = 0;
        for(j=start;j<PARAM_N-1;j+=2*distance)
        {
            W = omega[jTwiddle++];
            temp = a[j];
            a[j] = barrett_reduce((temp + a[j + distance]));
            a[j + distance] = montgomery_reduce((W * ((uint32_t)temp + 3*PARAM_Q - a[j + distance])));
        }
    }

}

static void poly_frombytes(poly *r, const unsigned char *a) {
	int i;
	for (i = 0; i < PARAM_N / 4; i++) {
		r->coeffs[4 * i + 0] = a[7 * i + 0] | (((uint16_t) a[7 * i + 1] & 0x3f) << 8);
		r->coeffs[4 * i + 1] = (a[7 * i + 1] >> 6) | (((uint16_t) a[7 * i + 2]) << 2) | (((uint16_t) a[7 * i + 3] & 0x0f) << 10);
		r->coeffs[4 * i + 2] = (a[7 * i + 3] >> 4) | (((uint16_t) a[7 * i + 4]) << 4) | (((uint16_t) a[7 * i + 5] & 0x03) << 12);
		r->coeffs[4 * i + 3] = (a[7 * i + 5] >> 2) | (((uint16_t) a[7 * i + 6]) << 6);
	}
}

static void poly_tobytes(unsigned char *r, const poly *p) {
	int i;
	uint16_t t0, t1, t2, t3, m;
	int16_t c;
	for (i = 0; i < PARAM_N / 4; i++) {
		t0 = barrett_reduce(p->coeffs[4 * i + 0]); //Make sure that coefficients have only 14 bits
		t1 = barrett_reduce(p->coeffs[4 * i + 1]);
		t2 = barrett_reduce(p->coeffs[4 * i + 2]);
		t3 = barrett_reduce(p->coeffs[4 * i + 3]);

		m = t0 - PARAM_Q;
		c = m;
		c >>= 15;
		t0 = m ^ ((t0 ^ m) & c); // <Make sure that coefficients are in [0,q]

		m = t1 - PARAM_Q;
		c = m;
		c >>= 15;
		t1 = m ^ ((t1 ^ m) & c); // <Make sure that coefficients are in [0,q]

		m = t2 - PARAM_Q;
		c = m;
		c >>= 15;
		t2 = m ^ ((t2 ^ m) & c); // <Make sure that coefficients are in [0,q]

		m = t3 - PARAM_Q;
		c = m;
		c >>= 15;
		t3 = m ^ ((t3 ^ m) & c); // <Make sure that coefficients are in [0,q]

		r[7 * i + 0] = t0 & 0xff;
		r[7 * i + 1] = (t0 >> 8) | (t1 << 6);
		r[7 * i + 2] = (t1 >> 2);
		r[7 * i + 3] = (t1 >> 10) | (t2 << 4);
		r[7 * i + 4] = (t2 >> 4);
		r[7 * i + 5] = (t2 >> 12) | (t3 << 2);
		r[7 * i + 6] = (t3 >> 6);
	}
}

static void poly_uniform(poly *a, const unsigned char *seed) {
	unsigned int pos = 0, ctr = 0;
	uint16_t val;
	uint64_t state[OQS_SHA3_STATESIZE];
	unsigned int nblocks = 8;
	uint8_t buf[OQS_SHA3_SHAKE128_RATE * nblocks];

	OQS_SHA3_shake128_absorb(state, seed, ZARZAR_SEEDBYTES);

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

static void poly_getnoise(poly *r, OQS_RAND *rand)
{
    unsigned char buf[5*PARAM_N];
    uint32_t *tp_a = (uint32_t *) buf;
    uint16_t ones, twos;
    uint8_t *tp_b = buf + PARAM_N*4;
    int i;

	rand->rand_n(rand, buf, 5 * PARAM_N);

    for(i=0;i<PARAM_N;i++)
    {
        uint32_t rnd32 = tp_a[i];
        ones = popcnt32(rnd32 >> 8); // 24 bit
        twos = popcnt16((tp_b[i] << 8) | (rnd32 & 0xFF)); // 16 bit
        r->coeffs[i] = (twos<<1) + ones - 28 + PARAM_Q;
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
	mul_coefficients(r->coeffs, zarzar_psis_bitrev_montgomery);
	ntt((uint16_t *) r->coeffs, zarzar_omegas_montgomery);
}

static void poly_invntt(poly *r) {
	bitrev_vector(r->coeffs);
	ntt((uint16_t *) r->coeffs, zarzar_omegas_inv_montgomery);
	mul_coefficients(r->coeffs, zarzar_psis_inv_montgomery);
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
    uint8_t exkey[64];

	oqs_rand->rand_n(oqs_rand, key, 32);
    for (i = 0;i < 64;i++) {
        uint8_t m = (key[i>>1] >> ((i&1)<<2)) & 0xF;
        exkey[i] = encode_hamming(m);
    }
    for (i = 0;i < PARAM_N;i++) {
        rbit = (exkey[i>>3] >> (i & 0x7)) & 1;
        v->coeffs[i] = (((((uint32_t)(sigma->coeffs[i] + ((PARAM_Q>>1) & (-rbit)))) << 6) + (PARAM_Q>>1)) / PARAM_Q) & 0x3F;
    }
}

static void rec(unsigned char *key, const poly *v, const poly *sigma)
{
    int i, j, k;

    for (i = 0;i < 32;i++) {
        key[i] = 0;
    }
    for (i = 0;i < 64;i++) {
        uint32_t tmp_v[8];
        for (j = 0;j < 8;j++) {
            k = (i<<3) | j;
            tmp_v[j] = ((((((uint32_t)v->coeffs[k])*PARAM_Q) + 32) >> 6) - sigma->coeffs[k] + PARAM_Q) % PARAM_Q;
        }
        key[i>>1] |= decode_E8(tmp_v) << ((i & 1) << 2);
    }
}
