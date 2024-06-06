#ifndef NTT_32_H
#define NTT_32_H

#include <stdint.h>
#include <immintrin.h>

#define ntt_avx_32 KYBER_NAMESPACE(ntt_avx_32)
void ntt_avx_32(__m512i *r, const __m512i *qdata_32);

#define basemul_avx_32 KYBER_NAMESPACE(basemul_avx_32)
void basemul_avx_32(__m512i *r,
                 const __m512i *a,
                 const __m512i *b,
                 const __m512i *qdata_32);

#define invntt_avx_32 KYBER_NAMESPACE(invntt_avx_32)
void invntt_avx_32(__m512i *r, const __m512i *qdata_32);


#define ntttobytes_avx512_32 KYBER_NAMESPACE(ntttobytes_avx512_32)
void ntttobytes_avx512_32(uint8_t *r, const __m512i *a, const __m512i *qdata_32);
#define nttfrombytes_avx512_32 KYBER_NAMESPACE(nttfrombytes_avx512_32)
void nttfrombytes_avx512_32(__m512i *r, const uint8_t *a, const __m512i *qdata_32);

#define poly_formseqto32_AVX512 KYBER_NAMESPACE(poly_formseqto32_AVX512)
void poly_formseqto32_AVX512(__m512i *a, __m512i *t, __m512i *aseq, const __m512i *qdata_32);

#define pk_formseqfrom32_AVX512 KYBER_NAMESPACE(pk_formseqfrom32_AVX512)
void pk_formseqfrom32_AVX512(uint8_t *k, uint8_t *kseq, uint8_t *t, const __m512i *qdata_32);
#define pk_formseqto32_AVX512 KYBER_NAMESPACE(pk_formseqto32_AVX512)
void pk_formseqto32_AVX512(uint8_t *k, uint8_t *t, uint8_t *kseq, const __m512i *qdata_32);
#define sk_formseqfrom32_AVX512 KYBER_NAMESPACE(sk_formseqfrom32_AVX512)
void sk_formseqfrom32_AVX512(uint8_t *k, uint8_t *kseq, uint8_t *t, const __m512i *qdata_32);
#define sk_formseqto32_AVX512 KYBER_NAMESPACE(sk_formseqto32_AVX512)
void sk_formseqto32_AVX512(uint8_t *k, uint8_t *t, uint8_t *kseq, const __m512i *qdata_32);

#define cipher_formseqfrom32_AVX512 KYBER_NAMESPACE(cipher_formseqfrom32_AVX512)
void cipher_formseqfrom32_AVX512(uint8_t *c, uint8_t *cseq, uint8_t *t,  const __m512i *qdata_32);
#define cipher_formseqto32_AVX512 KYBER_NAMESPACE(cipher_formseqto32_AVX512)
void cipher_formseqto32_AVX512(uint8_t *c, uint8_t *t, uint8_t *cseq,  const __m512i *qdata_32);

#endif