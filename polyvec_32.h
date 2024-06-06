#ifndef POLYVEC_32_H
#define POLYVEC_32_H

#include <stdint.h>
#include "params.h"
#include "poly_32.h"

typedef struct{
  poly vec[KYBER_K];
} polyvec;
typedef struct{
  poly_32 vec[KYBER_K];
} polyvec_32;

#define polyvec_compress KYBER_NAMESPACE(polyvec_compress)
void polyvec_compress(uint8_t r[KYBER_POLYVECCOMPRESSEDBYTES+2], const polyvec_32 *a);
#define polyvec_decompress KYBER_NAMESPACE(polyvec_decompress)
void polyvec_decompress(polyvec_32 *r, const uint8_t a[KYBER_POLYVECCOMPRESSEDBYTES]);

#define polyvec_tobytes KYBER_NAMESPACE(polyvec_tobytes)
void polyvec_tobytes(uint8_t r[KYBER_POLYVECBYTES*32], const polyvec_32 *a);
#define polyvec_frombytes KYBER_NAMESPACE(polyvec_frombytes)
void polyvec_frombytes(polyvec_32 *r, const uint8_t a[KYBER_POLYVECBYTES*32]);

#define polyvec_ntt KYBER_NAMESPACE(polyvec_ntt)
void polyvec_ntt(polyvec_32 *r);
#define polyvec_invntt_tomont KYBER_NAMESPACE(polyvec_invntt_tomont)
void polyvec_invntt_tomont(polyvec_32 *r);

#define polyvec_basemul_acc_montgomery KYBER_NAMESPACE(polyvec_basemul_acc_montgomery)
void polyvec_basemul_acc_montgomery(poly_32 *r, const polyvec_32 *a, const polyvec_32 *b);

#define polyvec_reduce KYBER_NAMESPACE(polyvec_reduce)
void polyvec_reduce(polyvec_32 *r);

#define polyvec_add KYBER_NAMESPACE(polyvec_add)
void polyvec_add(polyvec_32 *r, const polyvec_32 *a, const polyvec_32 *b);

#endif
