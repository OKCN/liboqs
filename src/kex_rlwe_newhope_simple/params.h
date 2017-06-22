#ifndef PARAMS_H
#define PARAMS_H

#include <stdint.h>

#define PARAM_N 1024

#define PARAM_K 16 /* used in sampler */
#define PARAM_Q 12289

#define POLY_BYTES 1792
#define NEWHOPE_SIMPLE_SEEDBYTES 32
#define NEWHOPE_SIMPLE_RECBYTES 384

#define NEWHOPE_SIMPLE_SENDABYTES (POLY_BYTES + NEWHOPE_SIMPLE_SEEDBYTES)
#define NEWHOPE_SIMPLE_SENDBBYTES (POLY_BYTES + NEWHOPE_SIMPLE_RECBYTES)

extern uint16_t simple_bitrev_table[];
extern uint16_t simple_omegas_montgomery[];
extern uint16_t simple_omegas_inv_montgomery[];
extern uint16_t simple_psis_inv_montgomery[];
extern uint16_t simple_psis_bitrev_montgomery[];

#if defined(WINDOWS)
typedef unsigned __int16 uint16_t;
#endif

#endif
