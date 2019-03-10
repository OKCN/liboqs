#ifndef PARAMS_H
#define PARAMS_H

#include <stdint.h>

#define PARAM_N 1024

#define PARAM_Q 12289
#define PARAM_LOGQ 14
#define PARAM_T 3
#define PARAM_REC_BITNUM 5

#define POLY_BYTES 1408
#define E8_SEEDBYTES 32
#define E8_RECBYTES 640

#define E8_SENDABYTES (POLY_BYTES + E8_SEEDBYTES)
#define E8_SENDBBYTES (POLY_BYTES + E8_RECBYTES)

extern uint16_t E8_bitrev_table[];
extern uint16_t E8_omegas_montgomery[];
extern uint16_t E8_omegas_inv_montgomery[];
extern uint16_t E8_psis_inv_montgomery[];
extern uint16_t E8_psis_bitrev_montgomery[];

#if defined(WINDOWS)
typedef unsigned __int16 uint16_t;
#endif

#endif
