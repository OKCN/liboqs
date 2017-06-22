#include <stdint.h>

#include <oqs/sha3.h>

// clang-format off
// (order of include matters)
#include "precomp.c"
#include "poly.c"
// clang-format on

static void encode_a(unsigned char *r, const poly *pk, const unsigned char *seed) {
	int i;
	poly_tobytes(r, pk);
	for (i = 0; i < ZARZAR_SEEDBYTES; i++) {
		r[POLY_BYTES + i] = seed[i];
	}
}

static void decode_a(poly *pk, unsigned char *seed, const unsigned char *r) {
	int i;
	poly_frombytes(pk, r);
	for (i = 0; i < ZARZAR_SEEDBYTES; i++) {
		seed[i] = r[POLY_BYTES + i];
	}
}

static void encode_b(unsigned char *r, const poly *b, const poly *c)
{
    int i;
    uint16_t *r16 = (uint16_t *) r;
    for(i = 0;i < PARAM_N;i++) {
        r16[i] = (b->coeffs[i] << 2) | (c->coeffs[i] & 0x3); 
    }
    r = (unsigned char *)(r16 + PARAM_N);
    for(i = 0;i < PARAM_N/2;i++) {
        r[i] = (c->coeffs[i<<1]>>2) | ((c->coeffs[(i<<1)|1]>>2) << 4);
    }
}

static void decode_b(poly *b, poly *c, const unsigned char *r)
{
    int i;
    uint16_t *r16 = (uint16_t *) r;
    for(i = 0;i < PARAM_N;i++) {
        b->coeffs[i] = r16[i] >> 2;
        c->coeffs[i] = r16[i] & 0x3;
    }
    for(i = 0;i < PARAM_N/2;i++) {
        c->coeffs[i<<1] |= (r[PARAM_N*2 + i] & 0xF) << 2;
        c->coeffs[i<<1|1] |= (r[PARAM_N*2 + i] >> 4) << 2;
    }
}

static void gen_a(poly *a, const unsigned char *seed) {
	poly_uniform(a, seed);
}

// API FUNCTIONS

static void keygen(unsigned char *send, poly *sk, OQS_RAND *rand) {
	poly a, e, r, pk;
	unsigned char seed[ZARZAR_SEEDBYTES];

	rand->rand_n(rand, seed, ZARZAR_SEEDBYTES);

	gen_a(&a, seed);

	poly_getnoise(sk, rand);
	poly_ntt(sk);

	poly_getnoise(&e, rand);
	poly_ntt(&e);

	poly_pointwise(&r, sk, &a);
	poly_add(&pk, &e, &r);

	encode_a(send, &pk, seed);
}

static void sharedb(unsigned char *sharedkey, unsigned char *send, const unsigned char *received, OQS_RAND *rand) {
	poly sp, ep, v, a, pka, c, epp, bp;
	unsigned char seed[ZARZAR_SEEDBYTES];

	decode_a(&pka, seed, received);
	gen_a(&a, seed);

	poly_getnoise(&sp, rand);
	poly_ntt(&sp);
	poly_getnoise(&ep, rand);
	poly_ntt(&ep);

	poly_pointwise(&bp, &a, &sp);
	poly_add(&bp, &bp, &ep);

	poly_pointwise(&v, &pka, &sp);
	poly_invntt(&v);

	poly_getnoise(&epp, rand);
	poly_add(&v, &v, &epp);

    con(&c, sharedkey, &v, rand);

	encode_b(send, &bp, &c);

#ifndef STATISTICAL_TEST
	OQS_SHA3_sha3256(sharedkey, sharedkey, 32);
#endif
}

static void shareda(unsigned char *sharedkey, const poly *sk, const unsigned char *received) {
	poly v, bp, c;

	decode_b(&bp, &c, received);

	poly_pointwise(&v, sk, &bp);
	poly_invntt(&v);

	rec(sharedkey, &c, &v);

#ifndef STATISTICAL_TEST
	OQS_SHA3_sha3256(sharedkey, sharedkey, 32);
#endif
}
