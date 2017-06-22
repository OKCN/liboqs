#ifndef _OQS_KEX_LWE_OKCN_LOCAL_H_
#define _OQS_KEX_LWE_OKCN_LOCAL_H_

#include <stdint.h>

#include <oqs/rand.h>

struct oqs_kex_lwe_okcn_params {
	uint8_t *seed;
	size_t seed_len;
	char *param_name;
	uint16_t log2_q;
	uint16_t q;
	uint16_t n;
	uint16_t extracted_bits;
	uint16_t nbar;
	uint16_t key_bits;
	uint16_t rec_hint_len;
	uint16_t single_hint_len;
	uint32_t pub_len;
	uint16_t stripe_step;
	int sampler_num;
	uint16_t *cdf_table;
	size_t cdf_table_len;
};

void oqs_kex_lwe_okcn_con_recommended(unsigned char *bob_rec, unsigned char *key, uint16_t *in);
void oqs_kex_lwe_okcn_rec_recommended(unsigned char *out, uint16_t *w, const unsigned char *hint);

void oqs_kex_lwe_okcn_pack(unsigned char *out, const size_t outlen, const uint16_t *in, const size_t inlen, const unsigned char lsb);
void oqs_kex_lwe_okcn_unpack(uint16_t *out, const size_t outlen, const unsigned char *in, const size_t inlen, const unsigned char lsb);

int oqs_kex_lwe_okcn_sample_n(uint16_t *s, const size_t n, struct oqs_kex_lwe_okcn_params *params, OQS_RAND *rand);

int oqs_kex_lwe_okcn_mul_add_as_plus_e_on_the_fly_recommended(uint16_t *b, const uint16_t *s, const uint16_t *e, struct oqs_kex_lwe_okcn_params *params);
int oqs_kex_lwe_okcn_mul_add_sa_plus_e_on_the_fly_recommended(uint16_t *b, const uint16_t *s, const uint16_t *e, struct oqs_kex_lwe_okcn_params *params);
void oqs_kex_lwe_okcn_mul_add_sb_plus_e_recommended(uint16_t *out, const uint16_t *b, const uint16_t *s, const uint16_t *e);
void oqs_kex_lwe_okcn_mul_bs_recommended(uint16_t *out, const uint16_t *b, const uint16_t *s);

#endif /* _OQS_KEX_RLWE_BCNS15_LOCAL_H_ */
