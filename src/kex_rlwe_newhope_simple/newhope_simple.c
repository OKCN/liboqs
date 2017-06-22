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
	for (i = 0; i < NEWHOPE_SIMPLE_SEEDBYTES; i++) {
		r[POLY_BYTES + i] = seed[i];
	}
}

static void decode_a(poly *pk, unsigned char *seed, const unsigned char *r) {
	int i;
	poly_frombytes(pk, r);
	for (i = 0; i < NEWHOPE_SIMPLE_SEEDBYTES; i++) {
		seed[i] = r[POLY_BYTES + i];
	}
}

static void encode_b(unsigned char *r, const poly *b, const poly *c) {
	int i, ptr = POLY_BYTES;
	poly_tobytes(r, b);
    uint16_t stack = 0;
    int top = 0;
	for (i = 0; i < PARAM_N; i++) {
        stack |= c->coeffs[i] << top;
        top += 3;
        if (top >= 8) {
            r[ptr] = stack & 0xFF;
            ptr ++;
            stack >>= 8;
            top -= 8;
        }
	}
}

static void decode_b(poly *b, poly *c, const unsigned char *r) {
	int i, ptr = 0, top = 0;
    uint16_t stack = 0;
	poly_frombytes(b, r);
	for (i = 0; i < NEWHOPE_SIMPLE_RECBYTES; i++) {
        stack |= r[POLY_BYTES + i] << top;
        top += 8;
        while (top >= 3) {
            c->coeffs[ptr] = stack & 7;
            ptr ++;
            stack >>= 3;
            top -= 3;
        }
	}
}

static void gen_a(poly *a, const unsigned char *seed) {
	poly_uniform(a, seed);
}

// API FUNCTIONS

static void keygen(unsigned char *send, poly *sk, OQS_RAND *rand) {
	poly a, e, r, pk;
	unsigned char seed[NEWHOPE_SIMPLE_SEEDBYTES];

	rand->rand_n(rand, seed, NEWHOPE_SIMPLE_SEEDBYTES);

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
	poly sp, ep, v, a, pka, c, epp, bp, enc_key;
	unsigned char seed[NEWHOPE_SIMPLE_SEEDBYTES];

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

	rand->rand_n(rand, sharedkey, 32);
    NHSEncode(&enc_key, sharedkey);
    poly_add(&v, &enc_key, &v);

    NHSCompress(&c, &v);
	encode_b(send, &bp, &c);

#ifndef STATISTICAL_TEST
	OQS_SHA3_sha3256(sharedkey, sharedkey, 32);
#endif
}

static void shareda(unsigned char *sharedkey, const poly *sk, const unsigned char *received) {
	poly v, bp, c, cp;

	decode_b(&bp, &c, received);

	poly_pointwise(&v, sk, &bp);
	poly_invntt(&v);

    NHSDecompress(&cp, &c);
    poly_substract(&cp, &cp, &v);
    NHSDecode(sharedkey, &cp);
#ifndef STATISTICAL_TEST
	OQS_SHA3_sha3256(sharedkey, sharedkey, 32);
#endif
}
