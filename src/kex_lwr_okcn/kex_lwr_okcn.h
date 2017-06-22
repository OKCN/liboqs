/**
 * \file kex_lwr_okcn.h
 * \brief Header for LWR key exchange protocol OKCN
 */

#ifndef __OQS_KEX_LWR_OKCN_H
#define __OQS_KEX_LWR_OKCN_H

#include <stddef.h>
#include <stdint.h>

#include <oqs/kex.h>
#include <oqs/rand.h>

OQS_KEX *OQS_KEX_lwr_okcn_new(OQS_RAND *rand, const uint8_t *seed, const size_t seed_len, const char *named_parameters);

int OQS_KEX_lwr_okcn_alice_0_recommended(OQS_KEX *k, void **alice_priv, uint8_t **alice_msg, size_t *alice_msg_len);
int OQS_KEX_lwr_okcn_bob_recommended(OQS_KEX *k, const uint8_t *alice_msg, const size_t alice_msg_len, uint8_t **bob_msg, size_t *bob_msg_len, uint8_t **key, size_t *key_len);
int OQS_KEX_lwr_okcn_alice_1_recommended(OQS_KEX *k, const void *alice_priv, const uint8_t *bob_msg, const size_t bob_msg_len, uint8_t **key, size_t *key_len);

void OQS_KEX_lwr_okcn_alice_priv_free(OQS_KEX *k, void *alice_priv);
void OQS_KEX_lwr_okcn_free(OQS_KEX *k);

#endif
