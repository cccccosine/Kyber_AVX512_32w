#ifndef INDCPA_32_H
#define INDCPA_32_H

#include <stdint.h>
#include "params.h"
#include "polyvec_32.h"

#define matrix_formseqto32 KYBER_NAMESPACE(matrix_formseqto32)
void matrix_formseqto32(polyvec_32 *a, polyvec_32 *t, polyvec_32 *aseq);
#define polyvec_formseqto32 KYBER_NAMESPACE(polyvec_formseqto32)
void polyvec_formseqto32(polyvec_32 *pv, polyvec_32 *t, polyvec_32 *pvseq);
#define poly_formseqto32 KYBER_NAMESPACE(poly_formseqto32)
void poly_formseqto32(poly_32 *p, poly_32 *t, poly_32 *pseq);
#define pk_formseqfrom32 KYBER_NAMESPACE(pk_formseqfrom32)
void pk_formseqfrom32(uint8_t *keyseq, uint8_t *t, uint8_t *key);
#define pk_formseqto32 KYBER_NAMESPACE(pk_formseqto32)
void pk_formseqto32(uint8_t *key, uint8_t *t, uint8_t *keyseq);
#define sk_formseqfrom32 KYBER_NAMESPACE(sk_formseqfrom32)
void sk_formseqfrom32(uint8_t *keyseq, uint8_t *t, uint8_t *key);
#define sk_formseqto32 KYBER_NAMESPACE(sk_formseqto32)
void sk_formseqto32(uint8_t *key, uint8_t *t, uint8_t *keyseq);
#define msg_formseqto32 KYBER_NAMESPACE(msg_formseqto32)
void msg_formseqto32(const uint8_t *m, uint8_t *mseq);
#define msg_formseqfrom32 KYBER_NAMESPACE(msg_formseqfrom32)
void msg_formseqfrom32(uint8_t *mseq, uint8_t *m);
#define cipher_formseqfrom32 KYBER_NAMESPACE(cipher_formseqfrom32)
void cipher_formseqfrom32(uint8_t *cseq, uint8_t *t, uint8_t *c);
#define cipher_formseqto32 KYBER_NAMESPACE(cipher_formseqto32)
void cipher_formseqto32(uint8_t *c, uint8_t *t, uint8_t *cseq);

#define gen_matrix KYBER_NAMESPACE(gen_matrix)
void gen_matrix(polyvec_32 *a, const uint8_t seed[KYBER_SYMBYTES*2*32], int transposed);
#define indcpa_keypair KYBER_NAMESPACE(indcpa_keypair)
void indcpa_keypair(uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                    uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES]
                    );

#define indcpa_enc KYBER_NAMESPACE(indcpa_enc)
void indcpa_enc(uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t m[KYBER_INDCPA_MSGBYTES*32*2],
                uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                uint8_t coins[KYBER_SYMBYTES*(32*2-1)]
                );

#define indcpa_dec KYBER_NAMESPACE(indcpa_dec)
void indcpa_dec(uint8_t m[KYBER_INDCPA_MSGBYTES*32*2],
                uint8_t c[KYBER_INDCPA_BYTES],
                uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES]
                );

#endif
