int MACRIFY(OQS_KEX_lwr_okcn_alice_0)(
        OQS_KEX *k, 
        void **alice_priv, 
        uint8_t **alice_msg, 
        size_t *alice_msg_len) 
{
	int ret;

	struct oqs_kex_lwr_okcn_params *params = (struct oqs_kex_lwr_okcn_params *) k->params;
	uint16_t *b = NULL;

	*alice_priv = NULL;
	*alice_msg = NULL;

	/* allocate private key and outgoing message */
	*alice_priv = malloc(params->n * params->nbar * sizeof(uint16_t));
	if (*alice_priv == NULL) {
		goto err;
	}
	b = (uint16_t *) malloc(params->n * params->nbar * sizeof(uint16_t));
	if (b == NULL) {
		goto err;
	}
	*alice_msg = malloc(params->pub_len);
	if (*alice_msg == NULL) {
		goto err;
	}

	/* generate S */
	ret = oqs_kex_lwr_okcn_sample_n(*alice_priv, params->n * params->nbar, params, k->rand);
	if (ret != 1) {
		goto err;
	}

	/* compute B = round(AS)_p */
	ret = MACRIFY(oqs_kex_lwr_okcn_mul_round_as_on_the_fly)(b, *alice_priv, params);
	if (ret != 1) {
		goto err;
	}
	oqs_kex_lwr_okcn_pack(*alice_msg, params->pub_len, b, params->n * params->nbar, params->log2_p);

	*alice_msg_len = params->pub_len;

	ret = 1;
	goto cleanup;

err:
	ret = 0;
	free(*alice_msg);
	*alice_msg = NULL;
	free(*alice_priv);
	*alice_priv = NULL;

cleanup:
	free(b);
	return ret;

}

int MACRIFY(OQS_KEX_lwr_okcn_bob)(
        OQS_KEX *k, 
        const uint8_t *alice_msg, 
        const size_t alice_msg_len, 
        uint8_t **bob_msg, 
        size_t *bob_msg_len, 
        uint8_t **key, 
        size_t *key_len) 
{
	int ret;

	struct oqs_kex_lwr_okcn_params *params = (struct oqs_kex_lwr_okcn_params *) k->params;

	uint16_t *bob_priv = NULL;
	uint8_t *bob_rec = NULL;
	uint16_t *b = NULL, *bprime = NULL;
	uint16_t *v = NULL;
	*bob_msg = NULL;
	*key = NULL;

	/* check length of other party's public key */
	if (alice_msg_len != params->pub_len) {
		goto err;
	}

	/* allocate private key outgoing message, and key */
	bob_priv = malloc(params->n * params->nbar * sizeof(uint16_t));
	if (bob_priv == NULL) {
		goto err;
	}
	bprime = (uint16_t *) malloc(params->n * params->nbar * sizeof(uint16_t));
	if (bprime == NULL) {
		goto err;
	}
	b = (uint16_t *) malloc(params->n * params->nbar * sizeof(uint16_t));
	if (b == NULL) {
		goto err;
	}
	v = (uint16_t *) malloc(params->nbar * params->nbar * sizeof(uint16_t));
	if (v == NULL) {
		goto err;
	}
    uint16_t bob_pub_len = LWR_DIV_ROUNDUP(params->n * params->nbar * params->log2_p, 8);
	*bob_msg = malloc(bob_pub_len + params->rec_hint_len);
	if (*bob_msg == NULL) {
		goto err;
	}
	bob_rec = *bob_msg + bob_pub_len;
	*key = malloc(params->key_bits >> 3);
	if (*key == NULL) {
		goto err;
	}

	/* generate S' */
	ret = oqs_kex_lwr_okcn_sample_n(bob_priv, params->n * params->nbar, params, k->rand);
	if (ret != 1) {
		goto err;
	}
	/* compute B' = round(S'A)_p */
	ret = MACRIFY(oqs_kex_lwr_okcn_mul_round_sa_on_the_fly)(bprime, bob_priv, params);
	if (ret != 1) {
		goto err;
	}
	oqs_kex_lwr_okcn_pack(*bob_msg, bob_pub_len, bprime, params->n * params->nbar, params->log2_p);

	/* unpack B */
	oqs_kex_lwr_okcn_unpack(b, params->n * params->nbar, alice_msg, alice_msg_len, params->log2_p);

	/* compute V = round(S'(B + Eps))_p */
	MACRIFY(oqs_kex_lwr_okcn_mul_round_sb)(v, b, bob_priv, k->rand);

    MACRIFY(oqs_kex_lwr_okcn_con)(bob_rec, (unsigned char *)(*key), v);

	*bob_msg_len = bob_pub_len + params->rec_hint_len;
	*key_len = params->key_bits >> 3;

	ret = 1;
	goto cleanup;

err:
	ret = 0;
	free(*bob_msg);
	*bob_msg = NULL;
	if (*key != NULL) {
		memset(*key, 0, params->key_bits >> 3);
	}
	free(*key);
	*key = NULL;

cleanup:
	free(bob_priv);
	free(bprime);
	free(b);
	if (v != NULL) {
		memset(v, 0, params->nbar * params->nbar * sizeof(uint16_t));
	}
	free(v);

	return ret;

}

int MACRIFY(OQS_KEX_lwr_okcn_alice_1)(OQS_KEX *k, const void *alice_priv, const uint8_t *bob_msg, const size_t bob_msg_len, uint8_t **key, size_t *key_len) {

	int ret;

	struct oqs_kex_lwr_okcn_params *params = (struct oqs_kex_lwr_okcn_params *) k->params;

	uint16_t *bprime = NULL, *w = NULL;
	*key = NULL;

	/* check length of other party's public key */
    uint16_t bob_pub_len = LWR_DIV_ROUNDUP(params->n * params->nbar * params->log2_p, 8);
	if (bob_msg_len != bob_pub_len + params->rec_hint_len) {
		goto err;
	}

	/* allocate working values and session key */
	bprime = malloc(params->n * params->nbar * sizeof(uint16_t));
	if (bprime == NULL) {
		goto err;
	}
	w = malloc(params->nbar * params->nbar * sizeof(uint16_t));
	if (w == NULL) {
		goto err;
	}
	*key = malloc(params->key_bits >> 3);
	if (*key == NULL) {
		goto err;
	}

	/* unpack B' */
	oqs_kex_lwr_okcn_unpack(bprime, params->n * params->nbar, bob_msg, bob_pub_len, params->log2_p);

	/* compute W = B'S */
	MACRIFY(oqs_kex_lwr_okcn_mul_bs)(w, bprime, (uint16_t *) alice_priv);

	/* compute K = rec(B'S, C) */
	const uint8_t *bob_rec = bob_msg + bob_pub_len;
	MACRIFY(oqs_kex_lwr_okcn_rec)(*key, w, bob_rec);

	*key_len = params->key_bits >> 3;

	ret = 1;
	goto cleanup;

err:
	ret = 0;
	memset(key, 0, params->key_bits >> 3);
	free(*key);
	*key = NULL;

cleanup:
	free(w);
	free(bprime);
	return ret;
}
