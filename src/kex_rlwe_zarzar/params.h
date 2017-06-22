#ifndef PARAMS_H
#define PARAMS_H

#include <stdint.h>

#define PARAM_N 512 

#define PARAM_Q 12289

#define POLY_BYTES 896
#define ZARZAR_SEEDBYTES 32
#define ZARZAR_RECBYTES 384

#define ZARZAR_SENDABYTES (POLY_BYTES + ZARZAR_SEEDBYTES)
#define ZARZAR_SENDBBYTES (POLY_BYTES + ZARZAR_RECBYTES)

extern uint16_t zarzar_bitrev_table[];
extern uint16_t zarzar_omegas_montgomery[];
extern uint16_t zarzar_omegas_inv_montgomery[];
extern uint16_t zarzar_psis_inv_montgomery[];
extern uint16_t zarzar_psis_bitrev_montgomery[];

#if defined(WINDOWS)
typedef unsigned __int16 uint16_t;
#endif

#endif
