#include <stdint.h>
#include <stdio.h>

#include <oqs/sha3.h>

// clang-format off
// (order of include matters)
#include "precomp.c"
#include "poly.c"
// clang-format on

static void encode_a(unsigned char *r, poly *pk, const unsigned char *seed) {
	int i;

    poly_reduce(pk);
    //for (i = 0;i <  10;i++) {
    //    printf("%d ", pk->coeffs[i]);
    //}
    poly_cutTail(pk, PARAM_T);
    frodo_pack(r, POLY_BYTES, pk->coeffs, PARAM_N, PARAM_LOGQ - PARAM_T);
    //puts("@@@@");
	for (i = 0; i < E8_SEEDBYTES; i++) {
		r[POLY_BYTES + i] = seed[i];
	}
}

static void decode_a(poly *pk, unsigned char *seed, const unsigned char *r) {
	int i;

    frodo_unpack(pk->coeffs, PARAM_N, r, POLY_BYTES, PARAM_LOGQ - PARAM_T);
    poly_extTail(pk, PARAM_T);

    //for (i = 0;i <  10;i++) {
    //    printf("%d ", pk->coeffs[i]);
    //}
    //puts("@@@@");
	for (i = 0; i < E8_SEEDBYTES; i++) {
		seed[i] = r[POLY_BYTES + i];
	}
}

static void encode_b(unsigned char *r, poly *b, const poly *c)
{
    //int i;
    poly_reduce(b);
    //for (i = 0;i < 20;i++) {
    //    printf("%d ", b->coeffs[i]);
    //}
    //puts("QQQQ");
    poly_cutTail(b, PARAM_T);
    frodo_pack(r, POLY_BYTES, b->coeffs, PARAM_N, PARAM_LOGQ - PARAM_T);
    frodo_pack(r + POLY_BYTES, E8_RECBYTES, c->coeffs, PARAM_N, PARAM_REC_BITNUM);
}

static void decode_b(poly *b, poly *c, const unsigned char *r)
{
    frodo_unpack(b->coeffs, PARAM_N, r, POLY_BYTES, PARAM_LOGQ - PARAM_T);
    frodo_unpack(c->coeffs, PARAM_N, r + POLY_BYTES, E8_RECBYTES, PARAM_REC_BITNUM);
    poly_extTail(b, PARAM_T);

    //int i;
    //for (i = 0;i < 20;i++) {
    //    printf("%d ", b->coeffs[i]);
    //}
    //puts("QQQQ");
}

static void gen_a(poly *a, const unsigned char *seed) {
	poly_uniform(a, seed);
}

// API FUNCTIONS

static void keygen(unsigned char *send, poly *sk, OQS_RAND *rand) {
	poly a, e, r, pk;
	unsigned char seed[E8_SEEDBYTES];

	rand->rand_n(rand, seed, E8_SEEDBYTES);

	gen_a(&a, seed);

	poly_getnoise(sk, rand);
    //int i;
    //for (i = 0;i < PARAM_N;i++) {
    //    printf("%d ", sk->coeffs[i]);
    //}
    //puts("sk");
	poly_ntt(sk);
    //for (i = 0;i < 20;i++) {
    //    printf("%d ", sk->coeffs[i]);
    //}
    //puts("sk");

	poly_getnoise(&e, rand);

	poly_pointwise(&r, sk, &a);
    poly_invntt(&r);

	poly_add(&pk, &e, &r);
    
	encode_a(send, &pk, seed);
}

static void sharedb(unsigned char *sharedkey, unsigned char *send, const unsigned char *received, OQS_RAND *rand) {
	poly sp, ep, v, a, pka, c, epp, bp;
	unsigned char seed[E8_SEEDBYTES];

	decode_a(&pka, seed, received);
    poly_ntt(&pka);
	gen_a(&a, seed);

	poly_getnoise(&sp, rand);
	poly_ntt(&sp);
	poly_getnoise(&ep, rand);

	poly_pointwise(&bp, &a, &sp);
    poly_invntt(&bp);
	poly_add(&bp, &bp, &ep);

	poly_pointwise(&v, &pka, &sp);
	poly_invntt(&v);

	poly_getnoise(&epp, rand);
	poly_add(&v, &v, &epp);

    con(&c, sharedkey, &v, rand);

	encode_b(send, &bp, &c);

#ifndef STATISTICAL_TEST
	//OQS_SHA3_sha3256(sharedkey, sharedkey, 64);
#endif
}

static void shareda(unsigned char *sharedkey, const poly *sk, const unsigned char *received) {
	poly v, bp, c;

	decode_b(&bp, &c, received);
    poly_ntt(&bp);

	poly_pointwise(&v, sk, &bp);
	poly_invntt(&v);

	rec(sharedkey, &c, &v);

#ifndef STATISTICAL_TEST
	//OQS_SHA3_sha3256(sharedkey, sharedkey, 64);
#endif
}
