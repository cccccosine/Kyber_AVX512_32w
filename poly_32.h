#ifndef POLY_32_H
#define POLY_32_H

#include <stdint.h>
#include "align.h"
#include "params.h"

typedef ALIGNED_INT16(KYBER_N) poly;
typedef ALIGNED_INT16_512(KYBER_N*32) poly_32;

#define poly_compress KYBER_NAMESPACE(poly_compress)
void poly_compress(uint8_t r[KYBER_POLYCOMPRESSEDBYTES], const poly_32 *a);
#define poly_decompress KYBER_NAMESPACE(poly_decompress)
void poly_decompress(poly_32 *r, const uint8_t a[KYBER_POLYCOMPRESSEDBYTES]);

#define poly_tobytes KYBER_NAMESPACE(poly_tobytes)
void poly_tobytes(uint8_t r[KYBER_POLYBYTES*32], const poly_32 *a);
#define poly_frombytes KYBER_NAMESPACE(poly_frombytes)
void poly_frombytes(poly_32 *r, const uint8_t a[KYBER_POLYBYTES*32]);

#define poly_frommsg_32 KYBER_NAMESPACE(poly_frommsg_32)
void poly_frommsg_32(poly_32 *r, const uint8_t msg[KYBER_INDCPA_MSGBYTES*32]);
#define poly_tomsg_32 KYBER_NAMESPACE(poly_tomsg_32)
void poly_tomsg_32(uint8_t msg[KYBER_INDCPA_MSGBYTES*32], const poly_32 *r);

#define poly_getnoise_eta1 KYBER_NAMESPACE(poly_getnoise_eta1)
void poly_getnoise_eta1(poly *r, const uint8_t seed[KYBER_SYMBYTES], uint8_t nonce);

#define poly_getnoise_eta2 KYBER_NAMESPACE(poly_getnoise_eta2)
void poly_getnoise_eta2(poly *r, const uint8_t seed[KYBER_SYMBYTES], uint8_t nonce);

#ifndef KYBER_90S
#define poly_getnoise_eta1_8x KYBER_NAMESPACE(poly_getnoise_eta1_8x)
void poly_getnoise_eta1_8x(poly_32 *r0,
                           poly_32 *r1,
                           poly_32 *r2,
                           poly_32 *r3,
                           poly_32 *r4,
                           poly_32 *r5,
                           poly_32 *r6,
                           poly_32 *r7,
                           uint8_t seed[32*(32*2-1)],
                           uint8_t nonce0,
                           uint8_t nonce1,
                           uint8_t nonce2,
                           uint8_t nonce3,
                           uint8_t nonce4,
                           uint8_t nonce5,
                           uint8_t nonce6,
                           uint8_t nonce7);

#define poly_getnoise_eta2_4x KYBER_NAMESPACE(poly_getnoise_eta2_4x)
void poly_getnoise_eta2_4x(poly_32 *r, const uint8_t seed[KYBER_SYMBYTES*(16*2-1)], uint8_t nonce);

#if KYBER_K == 2
#define poly_getnoise_eta1122_4x KYBER_NAMESPACE(poly_getnoise_eta1122_4x)
void poly_getnoise_eta1122_4x(poly_16 *r0,
                              poly_16 *r1,
                              poly_16 *r2,
                              poly_16 *r3,
                              const uint8_t seed[32*(16*2-1)],
                              uint8_t nonce0,
                              uint8_t nonce1,
                              uint8_t nonce2,
                              uint8_t nonce3);

#endif
#endif


#define poly_ntt KYBER_NAMESPACE(poly_ntt)
void poly_ntt(poly_32 *r);
#define poly_invntt_tomont KYBER_NAMESPACE(poly_invntt_tomont)
void poly_invntt_tomont(poly_32 *r);
#define poly_nttunpack KYBER_NAMESPACE(poly_nttunpack)
void poly_nttunpack(poly_32 *r);
#define poly_basemul_montgomery KYBER_NAMESPACE(poly_basemul_montgomery)
void poly_basemul_montgomery(poly_32 *r, const poly_32 *a, const poly_32 *b);
#define poly_tomont KYBER_NAMESPACE(poly_tomont)
void poly_tomont(poly_32 *r);

#define poly_reduce KYBER_NAMESPACE(poly_reduce)
void poly_reduce(poly_32 *r);

#define poly_add KYBER_NAMESPACE(poly_add)
void poly_add(poly_32 *r, const poly_32 *a, const poly_32 *b);
#define poly_sub KYBER_NAMESPACE(poly_sub)
void poly_sub(poly_32 *r, const poly_32 *a, const poly_32 *b);

#endif
