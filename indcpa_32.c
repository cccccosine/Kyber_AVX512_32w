#include <stddef.h>
#include <stdint.h>
#include <immintrin.h>
#include <string.h>
#include <malloc.h>
#include "align.h"
#include "params.h"
#include "consts_32.h"
#include "indcpa_32.h"
#include "polyvec_32.h"
#include "poly_32.h"
#include "ntt_32.h"
#include "cbd.h"
#include "rejsample.h"
#include "symmetric.h"
#include "randombytes.h"

void matrix_formseqto32(polyvec_32 *a, polyvec_32 *t, polyvec_32 *aseq) {
  for(int i = 0; i < KYBER_K; i++) {
    for(int j = 0; j < KYBER_K; j++) {
      poly_formseqto32_AVX512(a[i].vec[j].vec, t[i].vec[j].vec, aseq[i].vec[j].vec, qdata_32.vec);
    }
  }
}

void polyvec_formseqto32(polyvec_32 *pv, polyvec_32 *t, polyvec_32 *pvseq) {
  for(int i = 0; i < KYBER_K; i++) {
    poly_formseqto32_AVX512(pv->vec[i].vec, t->vec[i].vec, pvseq->vec[i].vec, qdata_32.vec);
  }
}

void poly_formseqto32(poly_32 *p, poly_32 *t, poly_32 *pseq) {
  poly_formseqto32_AVX512(p->vec, t->vec, pseq->vec, qdata_32.vec);
}

void pk_formseqfrom32(uint8_t *keyseq, uint8_t *t, uint8_t *key) {
  pk_formseqfrom32_AVX512(key, keyseq, t, qdata_32.vec);  //这里参数位置不同是为了适应汇编中的宏函数，想和to16共用formseq_supp宏函数
}

void pk_formseqto32(uint8_t *key, uint8_t *t, uint8_t *keyseq) {
  pk_formseqto32_AVX512(key, t, keyseq, qdata_32.vec);
}

void sk_formseqfrom32(uint8_t *keyseq, uint8_t *t, uint8_t *key) {
  sk_formseqfrom32_AVX512(key, keyseq, t, qdata_32.vec);
}

void sk_formseqto32(uint8_t *key, uint8_t *t, uint8_t *keyseq) {
  sk_formseqto32_AVX512(key, t, keyseq, qdata_32.vec);
}

void msg_formseqto32(const uint8_t *m, uint8_t *mseq) {  //目前不考虑msg的from/to32变换，因为msg本身是uint8_t类型，不是很适配AVX的16bit运算
  for(int i = 0; i < 32; i++) {
    for(int j = 0; j < 32; j++) {
      mseq[i*32+j] = m[j*64+i];  //kem中的每个单路msg后面还包括了H(pk),所以总长度是32*32*2
    }
  }
}

void msg_formseqfrom32(uint8_t *mseq, uint8_t *m) {
  for(int i = 0; i < 32; i++) {
    for(int j = 0; j < 32; j++) {
      m[j*64+i] = mseq[i*32+j];  //要间隔留出空间来连接kem后续的H(pk)
    }
  }
}

void cipher_formseqfrom32(uint8_t *cseq, uint8_t *t, uint8_t *c) {
  // for(int k = 0; k < 16; k++) {
  //   for(int i = 0; i < KYBER_K; i++) {
  //     for(int j = 0; j < 160; j++) { 
  //       c[k*(KYBER_K*320+128)+i*320+j*2] = cseq[i*16*320+j*16*2+k*2];
  //       c[k*(KYBER_K*320+128)+i*320+j*2+1] = cseq[i*16*320+j*16*2+k*2+1];
  //     }
  //   }

  //   for(int i = 0; i < 64; i++) {
  //     c[k*(KYBER_K*320+128)+KYBER_K*320+i*2] = cseq[KYBER_K*320*16+k*2+i*32];
  //     c[k*(KYBER_K*320+128)+KYBER_K*320+i*2+1] = cseq[KYBER_K*320*16+k*2+i*32+1];
  //   }
  // }

  cipher_formseqfrom32_AVX512(c, cseq, t, qdata_32.vec);
}

void cipher_formseqto32(uint8_t *c, uint8_t *t, uint8_t *cseq) {
  // for(int k = 0; k < 16; k++) {
  //   for(int i = 0; i < KYBER_K; i++) {
  //     for(int j = 0; j < 160; j++) { 
  //       cseq[i*16*320+j*16*2+k*2] = c[k*(KYBER_K*320+128)+i*320+j*2];
  //       cseq[i*16*320+j*16*2+k*2+1] = c[k*(KYBER_K*320+128)+i*320+j*2+1];
  //     }
  //   }
  //   for(int i = 0; i < 64; i++) {
  //     cseq[KYBER_K*320*16+k*2+i*32] = c[k*(KYBER_K*320+128)+KYBER_K*320+i*2];
  //     cseq[KYBER_K*320*16+k*2+i*32+1] = c[k*(KYBER_K*320+128)+KYBER_K*320+i*2+1];
  //   }
  // }

  cipher_formseqto32_AVX512(c, t, cseq, qdata_32.vec);
}

static void pack_pk(uint8_t r[KYBER_INDCPA_PUBLICKEYBYTES],
                    polyvec_32 *pk
                    // const uint8_t seed[KYBER_SYMBYTES]
                    )
{
  polyvec_tobytes(r, pk);
  // memcpy(r+KYBER_POLYVECBYTES*16, seed, KYBER_SYMBYTES);
}


static void unpack_pk(polyvec_32 *pk,
                      // uint8_t seed[KYBER_SYMBYTES],
                      const uint8_t packedpk[KYBER_INDCPA_PUBLICKEYBYTES]
                      )
{
  polyvec_frombytes(pk, packedpk);
  // memcpy(seed, packedpk+KYBER_POLYVECBYTES*16, KYBER_SYMBYTES);
}


static void pack_sk(uint8_t r[KYBER_INDCPA_SECRETKEYBYTES], polyvec_32 *sk)
{
  polyvec_tobytes(r, sk);
}


static void unpack_sk(polyvec_32 *sk, const uint8_t packedsk[KYBER_INDCPA_SECRETKEYBYTES])
{
  polyvec_frombytes(sk, packedsk);
}


static void pack_ciphertext(uint8_t r[KYBER_INDCPA_BYTES], polyvec_32 *b, poly_32 *v)
{
  polyvec_compress(r, b);
  poly_compress(r+KYBER_POLYVECCOMPRESSEDBYTES, v);
}


static void unpack_ciphertext(polyvec_32 *b, poly_32 *v, const uint8_t c[KYBER_INDCPA_BYTES])
{
  polyvec_decompress(b, c);
  poly_decompress(v, c+KYBER_POLYVECCOMPRESSEDBYTES);
}


static unsigned int rej_uniform(int16_t *r,
                                unsigned int len,
                                const uint8_t *buf,
                                unsigned int buflen)
{
  unsigned int ctr, pos;
  uint16_t val0, val1;

  ctr = pos = 0;
  while(ctr < len && pos <= buflen - 3) {  // buflen is always at least 3
    val0 = ((buf[pos+0] >> 0) | ((uint16_t)buf[pos+1] << 8)) & 0xFFF;
    val1 = ((buf[pos+1] >> 4) | ((uint16_t)buf[pos+2] << 4)) & 0xFFF;
    pos += 3;

    if(val0 < KYBER_Q)
      r[ctr++] = val0;
    if(ctr < len && val1 < KYBER_Q)
      r[ctr++] = val1;
  }

  return ctr;
}

#define gen_a(A,B)  gen_matrix(A,B,0)
#define gen_at(A,B) gen_matrix(A,B,1)


#ifdef KYBER_90S
void gen_matrix(polyvec *a, const uint8_t seed[32], int transposed)
{
  unsigned int ctr, i, j, k;
  unsigned int buflen, off;
  uint64_t nonce = 0;
  ALIGNED_UINT8(REJ_UNIFORM_AVX_NBLOCKS*AES256CTR_BLOCKBYTES) buf;
  aes256ctr_ctx state;

  aes256ctr_init(&state, seed, 0);

  for(i=0;i<KYBER_K;i++) {
    for(j=0;j<KYBER_K;j++) {
      if(transposed)
        nonce = (j << 8) | i;
      else
        nonce = (i << 8) | j;

      state.n = _mm_loadl_epi64((__m128i *)&nonce);
      aes256ctr_squeezeblocks(buf.coeffs, REJ_UNIFORM_AVX_NBLOCKS, &state);
      buflen = REJ_UNIFORM_AVX_NBLOCKS*AES256CTR_BLOCKBYTES;
      ctr = rej_uniform_avx(a[i].vec[j].coeffs, buf.coeffs);

      while(ctr < KYBER_N) {
        off = buflen % 3;
        for(k = 0; k < off; k++)
          buf.coeffs[k] = buf.coeffs[buflen - off + k];
        aes256ctr_squeezeblocks(buf.coeffs + off, 1, &state);
        buflen = off + AES256CTR_BLOCKBYTES;
        ctr += rej_uniform(a[i].vec[j].coeffs + ctr, KYBER_N - ctr, buf.coeffs, buflen);
      }

      poly_nttunpack(&a[i].vec[j]);
    }
  }
}
#else
#if KYBER_K == 2
void gen_matrix(polyvec_16 *a, const uint8_t seed[32*(2*16-1)], int transposed)
{
  unsigned int ctr0, ctr1, ctr2, ctr3;
  ALIGNED_UINT8(REJ_UNIFORM_AVX_NBLOCKS*SHAKE128_RATE) buf[4];
  __m256i f;
  keccakx4_state state;

  for(int i = 0; i < 16; i++) {
    f = _mm256_loadu_si256((__m256i *)(seed+i*32*2));
    _mm256_store_si256(buf[0].vec, f);
    _mm256_store_si256(buf[1].vec, f);
    _mm256_store_si256(buf[2].vec, f);
    _mm256_store_si256(buf[3].vec, f);

    if(transposed) {
      buf[0].coeffs[32] = 0;
      buf[0].coeffs[33] = 0;
      buf[1].coeffs[32] = 0;
      buf[1].coeffs[33] = 1;
      buf[2].coeffs[32] = 1;
      buf[2].coeffs[33] = 0;
      buf[3].coeffs[32] = 1;
      buf[3].coeffs[33] = 1;
    }
    else {
      buf[0].coeffs[32] = 0;
      buf[0].coeffs[33] = 0;
      buf[1].coeffs[32] = 1;
      buf[1].coeffs[33] = 0;
      buf[2].coeffs[32] = 0;
      buf[2].coeffs[33] = 1;
      buf[3].coeffs[32] = 1;
      buf[3].coeffs[33] = 1;
    }

    shake128x4_absorb_once(&state, buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 34);
    shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, REJ_UNIFORM_AVX_NBLOCKS, &state);

    ctr0 = rej_uniform_avx(a[0].vec[0].coeffs + i*KYBER_N, buf[0].coeffs);
    ctr1 = rej_uniform_avx(a[0].vec[1].coeffs + i*KYBER_N, buf[1].coeffs);
    ctr2 = rej_uniform_avx(a[1].vec[0].coeffs + i*KYBER_N, buf[2].coeffs);
    ctr3 = rej_uniform_avx(a[1].vec[1].coeffs + i*KYBER_N, buf[3].coeffs);

    while(ctr0 < KYBER_N || ctr1 < KYBER_N || ctr2 < KYBER_N || ctr3 < KYBER_N) {
      shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 1, &state);

      ctr0 += rej_uniform(a[0].vec[0].coeffs + i*KYBER_N + ctr0, KYBER_N - ctr0, buf[0].coeffs, SHAKE128_RATE);
      ctr1 += rej_uniform(a[0].vec[1].coeffs + i*KYBER_N + ctr1, KYBER_N - ctr1, buf[1].coeffs, SHAKE128_RATE);
      ctr2 += rej_uniform(a[1].vec[0].coeffs + i*KYBER_N + ctr2, KYBER_N - ctr2, buf[2].coeffs, SHAKE128_RATE);
      ctr3 += rej_uniform(a[1].vec[1].coeffs + i*KYBER_N + ctr3, KYBER_N - ctr3, buf[3].coeffs, SHAKE128_RATE);
    }
  }

}
#elif KYBER_K == 3
void gen_matrix(polyvec_32 *a, const uint8_t seed[KYBER_SYMBYTES*2*32], int transposed)
{
  unsigned int ctr0, ctr1, ctr2, ctr3, ctr4, ctr5, ctr6, ctr7;
  ALIGNED_UINT8(REJ_UNIFORM_AVX_NBLOCKS*SHAKE128_RATE) buf[8];   //3*168
  __m256i f0, f1, f2, f3, f4, f5, f6, f7;
  keccakx8_state state;

  for(int i = 0; i < 9; i++) {
    for(int j = 0; j < 4; j++) {
      f0 = _mm256_loadu_si256((__m256i *)(seed+(j*8+0)*32*2));
      f1 = _mm256_loadu_si256((__m256i *)(seed+(j*8+1)*32*2));
      f2 = _mm256_loadu_si256((__m256i *)(seed+(j*8+2)*32*2));
      f3 = _mm256_loadu_si256((__m256i *)(seed+(j*8+3)*32*2));
      f4 = _mm256_loadu_si256((__m256i *)(seed+(j*8+4)*32*2));
      f5 = _mm256_loadu_si256((__m256i *)(seed+(j*8+5)*32*2));
      f6 = _mm256_loadu_si256((__m256i *)(seed+(j*8+6)*32*2));
      f7 = _mm256_loadu_si256((__m256i *)(seed+(j*8+7)*32*2));

      _mm256_store_si256(buf[0].vec, f0);
      _mm256_store_si256(buf[1].vec, f1);
      _mm256_store_si256(buf[2].vec, f2);
      _mm256_store_si256(buf[3].vec, f3);
      _mm256_store_si256(buf[4].vec, f4);
      _mm256_store_si256(buf[5].vec, f5);
      _mm256_store_si256(buf[6].vec, f6);
      _mm256_store_si256(buf[7].vec, f7);

      switch (i)
      {
      case 0:
        for(int k = 0; k < 8; k++) {
          if(transposed) {
            buf[k].coeffs[32] = 0;
            buf[k].coeffs[33] = 0;
          }
          else {
            buf[k].coeffs[32] = 0;
            buf[k].coeffs[33] = 0;
          }
        }
        break;
      case 1:
        for(int k = 0; k < 8; k++) {
          if(transposed) {
            buf[k].coeffs[32] = 0;
            buf[k].coeffs[33] = 1;
          }
          else {
            buf[k].coeffs[32] = 1;
            buf[k].coeffs[33] = 0;
          }
        }
        break;
      case 2:
        for(int k = 0; k < 8; k++) {
          if(transposed) {
            buf[k].coeffs[32] = 0;
            buf[k].coeffs[33] = 2;
          }
          else {
            buf[k].coeffs[32] = 2;
            buf[k].coeffs[33] = 0;
          }
        }
        break;
      case 3:
        for(int k = 0; k < 8; k++) {
          if(transposed) {
            buf[k].coeffs[32] = 1;
            buf[k].coeffs[33] = 0;
          }
          else {
            buf[k].coeffs[32] = 0;
            buf[k].coeffs[33] = 1;
          }
        }
        break;
      case 4:
        for(int k = 0; k < 8; k++) {
          if(transposed) {
            buf[k].coeffs[32] = 1;
            buf[k].coeffs[33] = 1;
          }
          else {
            buf[k].coeffs[32] = 1;
            buf[k].coeffs[33] = 1;
          }
        }
        break;
      case 5:
        for(int k = 0; k < 8; k++) {
          if(transposed) {
            buf[k].coeffs[32] = 1;
            buf[k].coeffs[33] = 2;
          }
          else {
            buf[k].coeffs[32] = 2;
            buf[k].coeffs[33] = 1;
          }
        }
        break;
      case 6:
        for(int k = 0; k < 8; k++) {
          if(transposed) {
            buf[k].coeffs[32] = 2;
            buf[k].coeffs[33] = 0;
          }
          else {
            buf[k].coeffs[32] = 0;
            buf[k].coeffs[33] = 2;
          }
        }
        break;
      case 7:
        for(int k = 0; k < 8; k++) {
          if(transposed) {
            buf[k].coeffs[32] = 2;
            buf[k].coeffs[33] = 1;
          }
          else {
            buf[k].coeffs[32] = 1;
            buf[k].coeffs[33] = 2;
          }
        }
        break;
      case 8:
        for(int k = 0; k < 8; k++) {
          if(transposed) {
            buf[k].coeffs[32] = 2;
            buf[k].coeffs[33] = 2;
          }
          else {
            buf[k].coeffs[32] = 2;
            buf[k].coeffs[33] = 2;
          }
        }
        break;
      default:
        printf("The i in loop is invalid!");
        break;
      }

      shake128x8_absorb_once(&state, buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, buf[4].coeffs, buf[5].coeffs, buf[6].coeffs, buf[7].coeffs, 34);
      shake128x8_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, buf[4].coeffs, buf[5].coeffs, buf[6].coeffs, buf[7].coeffs, REJ_UNIFORM_AVX_NBLOCKS, &state);

      ctr0 = rej_uniform_avx(a[i/3].vec[i%3].coeffs + (j*8+0)*KYBER_N, buf[0].coeffs); 
      ctr1 = rej_uniform_avx(a[i/3].vec[i%3].coeffs + (j*8+1)*KYBER_N, buf[1].coeffs); 
      ctr2 = rej_uniform_avx(a[i/3].vec[i%3].coeffs + (j*8+2)*KYBER_N, buf[2].coeffs); 
      ctr3 = rej_uniform_avx(a[i/3].vec[i%3].coeffs + (j*8+3)*KYBER_N, buf[3].coeffs); 
      ctr4 = rej_uniform_avx(a[i/3].vec[i%3].coeffs + (j*8+4)*KYBER_N, buf[4].coeffs); 
      ctr5 = rej_uniform_avx(a[i/3].vec[i%3].coeffs + (j*8+5)*KYBER_N, buf[5].coeffs); 
      ctr6 = rej_uniform_avx(a[i/3].vec[i%3].coeffs + (j*8+6)*KYBER_N, buf[6].coeffs); 
      ctr7 = rej_uniform_avx(a[i/3].vec[i%3].coeffs + (j*8+7)*KYBER_N, buf[7].coeffs); 

      while (ctr0 < KYBER_N || ctr1 < KYBER_N || ctr2 < KYBER_N || ctr3 < KYBER_N || ctr4 < KYBER_N || ctr5 < KYBER_N || ctr6 < KYBER_N || ctr7 < KYBER_N)
      {
        shake128x8_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, buf[4].coeffs, buf[5].coeffs, buf[6].coeffs, buf[7].coeffs, 1, &state); 

        ctr0 += rej_uniform(a[i/3].vec[i%3].coeffs + (j*8+0)*KYBER_N + ctr0, KYBER_N - ctr0, buf[0].coeffs, SHAKE128_RATE); 
        ctr1 += rej_uniform(a[i/3].vec[i%3].coeffs + (j*8+1)*KYBER_N + ctr1, KYBER_N - ctr1, buf[1].coeffs, SHAKE128_RATE); 
        ctr2 += rej_uniform(a[i/3].vec[i%3].coeffs + (j*8+2)*KYBER_N + ctr2, KYBER_N - ctr2, buf[2].coeffs, SHAKE128_RATE); 
        ctr3 += rej_uniform(a[i/3].vec[i%3].coeffs + (j*8+3)*KYBER_N + ctr3, KYBER_N - ctr3, buf[3].coeffs, SHAKE128_RATE); 
        ctr4 += rej_uniform(a[i/3].vec[i%3].coeffs + (j*8+4)*KYBER_N + ctr4, KYBER_N - ctr4, buf[4].coeffs, SHAKE128_RATE); 
        ctr5 += rej_uniform(a[i/3].vec[i%3].coeffs + (j*8+5)*KYBER_N + ctr5, KYBER_N - ctr5, buf[5].coeffs, SHAKE128_RATE); 
        ctr6 += rej_uniform(a[i/3].vec[i%3].coeffs + (j*8+6)*KYBER_N + ctr6, KYBER_N - ctr6, buf[6].coeffs, SHAKE128_RATE); 
        ctr7 += rej_uniform(a[i/3].vec[i%3].coeffs + (j*8+7)*KYBER_N + ctr7, KYBER_N - ctr7, buf[7].coeffs, SHAKE128_RATE);

      }
    }
  }

}
#elif KYBER_K == 4
void gen_matrix(polyvec_16 *a, const uint8_t seed[32*(2*16-1)], int transposed)
{
  unsigned int i, ctr0, ctr1, ctr2, ctr3;
  ALIGNED_UINT8(REJ_UNIFORM_AVX_NBLOCKS*SHAKE128_RATE) buf[4];
  __m256i f;
  keccakx4_state state;

  for(i=0;i<4;i++) {
    for(int j = 0; j < 16; j++) {
      f = _mm256_loadu_si256((__m256i *)(seed+j*32*2));
      _mm256_store_si256(buf[0].vec, f);
      _mm256_store_si256(buf[1].vec, f);
      _mm256_store_si256(buf[2].vec, f);
      _mm256_store_si256(buf[3].vec, f);

      if(transposed) {
        buf[0].coeffs[32] = i;
        buf[0].coeffs[33] = 0;
        buf[1].coeffs[32] = i;
        buf[1].coeffs[33] = 1;
        buf[2].coeffs[32] = i;
        buf[2].coeffs[33] = 2;
        buf[3].coeffs[32] = i;
        buf[3].coeffs[33] = 3;
      }
      else {
        buf[0].coeffs[32] = 0;
        buf[0].coeffs[33] = i;
        buf[1].coeffs[32] = 1;
        buf[1].coeffs[33] = i;
        buf[2].coeffs[32] = 2;
        buf[2].coeffs[33] = i;
        buf[3].coeffs[32] = 3;
        buf[3].coeffs[33] = i;
      }

      shake128x4_absorb_once(&state, buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 34);
      shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, REJ_UNIFORM_AVX_NBLOCKS, &state);

      ctr0 = rej_uniform_avx(a[i].vec[0].coeffs + j*KYBER_N, buf[0].coeffs);
      ctr1 = rej_uniform_avx(a[i].vec[1].coeffs + j*KYBER_N, buf[1].coeffs);
      ctr2 = rej_uniform_avx(a[i].vec[2].coeffs + j*KYBER_N, buf[2].coeffs);
      ctr3 = rej_uniform_avx(a[i].vec[3].coeffs + j*KYBER_N, buf[3].coeffs);

      while(ctr0 < KYBER_N || ctr1 < KYBER_N || ctr2 < KYBER_N || ctr3 < KYBER_N) {
        shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 1, &state);

        ctr0 += rej_uniform(a[i].vec[0].coeffs + j*KYBER_N + ctr0, KYBER_N - ctr0, buf[0].coeffs, SHAKE128_RATE);
        ctr1 += rej_uniform(a[i].vec[1].coeffs + j*KYBER_N + ctr1, KYBER_N - ctr1, buf[1].coeffs, SHAKE128_RATE);
        ctr2 += rej_uniform(a[i].vec[2].coeffs + j*KYBER_N + ctr2, KYBER_N - ctr2, buf[2].coeffs, SHAKE128_RATE);
        ctr3 += rej_uniform(a[i].vec[3].coeffs + j*KYBER_N + ctr3, KYBER_N - ctr3, buf[3].coeffs, SHAKE128_RATE);
      }
    }
    
  }
}
#endif
#endif


void indcpa_keypair(uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],     
                    uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES]
                    // uint16_t pkpvprint[KYBER_INDCPA_PUBLICKEYBYTES]
                    )
{
  unsigned int i;
  uint8_t buf[2*KYBER_SYMBYTES*32] = {0};
  uint8_t *pkseq = (uint8_t *)malloc(KYBER_INDCPA_PUBLICKEYBYTES);
  uint8_t *skseq = (uint8_t *)malloc(KYBER_INDCPA_SECRETKEYBYTES);
  uint8_t *tkp = (uint8_t *)malloc(KYBER_INDCPA_PUBLICKEYBYTES);
  const uint8_t *publicseed = buf;
  const uint8_t *noiseseed = buf + KYBER_SYMBYTES;
  polyvec_32 a[KYBER_K], aseq[KYBER_K], t[KYBER_K], skpv, skpvseq, tpv, e, eseq, pkpv, pkpvseq;

  randombytes(buf+KYBER_SYMBYTES*32, KYBER_SYMBYTES*32);  //只需要生成一半的随机数并放在后半部分的内存中，之后hash_gx8从后半部分开始取，生成的数从头开始存
  // printf("indcpa_keygen--randombytes1:");
  // for(int i = 0; i < KYBER_SYMBYTES*16; i++) {
  //   printf("%d,", *(buf + KYBER_SYMBYTES*16 + i));
  // }
  // printf("\n");

  // for(i = 0; i < KYBER_SYMBYTES*32; i++) {
  //   buf[i+KYBER_SYMBYTES*32] = 1;
  // }

  for(i = 0; i < 4; i++) {
    hash_gx8(buf+16*i*KYBER_SYMBYTES, 
             buf+(16*i+2)*KYBER_SYMBYTES, 
             buf+(16*i+4)*KYBER_SYMBYTES, 
             buf+(16*i+6)*KYBER_SYMBYTES, 
             buf+(16*i+8)*KYBER_SYMBYTES, 
             buf+(16*i+10)*KYBER_SYMBYTES, 
             buf+(16*i+12)*KYBER_SYMBYTES, 
             buf+(16*i+14)*KYBER_SYMBYTES,
             buf+KYBER_SYMBYTES*32+KYBER_SYMBYTES*i*8, 
             buf+KYBER_SYMBYTES*32+KYBER_SYMBYTES*(i*8+1), 
             buf+KYBER_SYMBYTES*32+KYBER_SYMBYTES*(i*8+2), 
             buf+KYBER_SYMBYTES*32+KYBER_SYMBYTES*(i*8+3), 
             buf+KYBER_SYMBYTES*32+KYBER_SYMBYTES*(i*8+4), 
             buf+KYBER_SYMBYTES*32+KYBER_SYMBYTES*(i*8+5), 
             buf+KYBER_SYMBYTES*32+KYBER_SYMBYTES*(i*8+6), 
             buf+KYBER_SYMBYTES*32+KYBER_SYMBYTES*(i*8+7),
             KYBER_SYMBYTES);
  }
  // printf("indcpa_keygen--hashgx4-1:");
  // for(int i = 0; i < KYBER_SYMBYTES*32*2; i++) {
  //   printf("%d,", *(buf + i));
  //   if (i % 16 == 15) printf("\n");
  // }
  // printf("\n");

  // for (int i = 0; i < KYBER_K; i++) {
  //   for (j = 0; j < KYBER_K; j++) {
  //     for(k = 0; k < KYBER_N; k++){
  //       for(p = 0; p < 16; p++) {
  //         a[i].vec[j].coeffs[k*16+p] = 19;
  //       }
  //     }
  //   }
  // }

  // for(int i = 0; i < 2*KYBER_SYMBYTES*16; i++) {
  //   buf[i] = 1;
  // }

  gen_a(a, publicseed);
  matrix_formseqto32(a, t, aseq);
  // printf("indcpa_keygen--gen_a:");
  // for (int i = 0; i < KYBER_K; i++) {
  //   for (int j = 0; j < KYBER_K; j++) {
  //     for(int k = 0; k < KYBER_N; k++){
  //       for(int p = 0; p < 32; p++) {
  //         printf("%5d,", aseq[i].vec[j].coeffs[k*32+p]);
  //         if (p % 16 == 15) printf("\n");
  //       }
  //       // printf("\n");
  //     }
  //     printf("\n////////////////////////////////\n");
  //   }
  //   printf("\n*******************************\n");
  // }
  // printf("\n");

#ifdef KYBER_90S  //not changed
#define NOISE_NBLOCKS ((KYBER_ETA1*KYBER_N/4)/AES256CTR_BLOCKBYTES) /* Assumes divisibility */
  uint64_t nonce = 0;
  ALIGNED_UINT8(NOISE_NBLOCKS*AES256CTR_BLOCKBYTES+32) coins; // +32 bytes as required by poly_cbd_eta1
  aes256ctr_ctx state;
  aes256ctr_init(&state, noiseseed, nonce++);
  for(i=0;i<KYBER_K;i++) {
    aes256ctr_squeezeblocks(coins.coeffs, NOISE_NBLOCKS, &state);
    state.n = _mm_loadl_epi64((__m128i *)&nonce);
    nonce += 1;
    poly_cbd_eta1(&skpv.vec[i], coins.vec);
  }
  for(i=0;i<KYBER_K;i++) {
    aes256ctr_squeezeblocks(coins.coeffs, NOISE_NBLOCKS, &state);
    state.n = _mm_loadl_epi64((__m128i *)&nonce);
    nonce += 1;
    poly_cbd_eta1(&e.vec[i], coins.vec);
  }
#else
#if KYBER_K == 2
  poly_getnoise_eta1_4x(skpv.vec+0, skpv.vec+1, e.vec+0, e.vec+1, noiseseed, 0, 1, 2, 3);
  polyvec_formseqto16(&skpv, &tpv, &skpvseq);
  polyvec_formseqto16(&e, &tpv, &eseq);  
#elif KYBER_K == 3 
  // for (int j = 0; j < KYBER_K; j++) {
  //   for(int k = 0; k < 256; k++) {
  //     for(int p = 0; p < 16; p++) {
  //       skpv.vec[j].coeffs[k*16+p] = 19;
  //       e.vec[j].coeffs[k*16+p] = 19;
  //     }
  //   }
  // }
  poly_getnoise_eta1_8x(skpv.vec+0, skpv.vec+1, skpv.vec+2, e.vec+0, e.vec+1, e.vec+2, pkpv.vec+0, pkpv.vec+1, noiseseed, 0, 1, 2, 3, 4, 5, 6, 7);
  polyvec_formseqto32(&skpv, &tpv, &skpvseq);
  polyvec_formseqto32(&e, &tpv, &eseq);
  // for (int j = 0; j < KYBER_K; j++) {
  //   for(int k = 0; k < KYBER_N; k++){
  //     for (int p = 0; p < 32; p++) {
  //       printf("%5d,", skpvseq.vec[j].coeffs[k*32+p]);
  //       if (p % 16 == 15) printf("\n");
  //     }
  //   }
  //   printf("\n\n");
  // }
#elif KYBER_K == 4    //not changed
  poly_getnoise_eta1_4x(skpv.vec+0, skpv.vec+1, skpv.vec+2, skpv.vec+3, noiseseed,  0, 1, 2, 3);
  poly_getnoise_eta1_4x(e.vec+0, e.vec+1, e.vec+2, e.vec+3, noiseseed, 4, 5, 6, 7);
#endif
#endif

  polyvec_ntt(&skpvseq);
  polyvec_reduce(&skpvseq);
  polyvec_ntt(&eseq);

  // for (int j = 0; j < KYBER_K; j++) {
  //   for(int k = 0; k < 256; k++) {
  //     for(int p = 0; p < 32; p++) {
  //       printf("%6d,", eseq.vec[j].coeffs[k*32+p]);
  //       if (p % 16 == 15) printf("\n");
  //     }
  //   }
  //   printf("\n/////////////////////////////\n");
  // }
  // printf("\n");

  // matrix-vector multiplication
  for(i=0;i<KYBER_K;i++) {
    polyvec_basemul_acc_montgomery(&pkpvseq.vec[i], &aseq[i], &skpvseq);
    poly_tomont(&pkpvseq.vec[i]);
  }

  polyvec_add(&pkpvseq, &pkpvseq, &eseq);
  polyvec_reduce(&pkpvseq);
  // for (int j = 0; j < 3; j++) {
  //   for(int k = 0; k < 256; k++) {
  //     for(int p = 0; p < 32; p++) {
  //       printf("%6d,", pkpvseq.vec[j].coeffs[k*32+p]);
  //       if (p % 16 == 15) printf("\n");
  //     }
  //   }
  //   printf("\n/////////////////////////////\n");
  // }
  // printf("\n");

  pack_sk(skseq, &skpvseq);
  pack_pk(pkseq, &pkpvseq);
  sk_formseqfrom32(skseq, tkp, sk);
  pk_formseqfrom32(pkseq, tkp, pk);

  for(i = 0; i < 32; i++) {
    memcpy(pk+KYBER_POLYVECBYTES+KYBER_INDCPA_PUBLICKEYBYTES/32*i, publicseed+KYBER_SYMBYTES*2*i, KYBER_SYMBYTES);
  }
  
  // for (int i = 0; i < KYBER_INDCPA_SECRETKEYBYTES; i++) {
  //   printf("%5d, ", sk[i]);
  //   if (i % 16 == 15) printf("\n");
  // }

  // for(i = 0; i < KYBER_K; i++) {
  //   // for(j = 0; j < KYBER_N; j++) {
  //   for(j = 0; j < 384; j++) {
  //     for(k = 0; k < 16; k++) {
  //       // skpvprint[i*KYBER_N+j] = skpv.vec[i].coeffs[j];
  //       // pkpvprint[(i*KYBER_N+j)*16+k] = skpvseq.vec[i].coeffs[j*16+k];
  //       // pkpvprint[(i*KYBER_N+j)*16+k] = pkpvseq.vec[i].coeffs[j*16+k];
  //       // pkpvprint[(i*384+j)*16+k] = sk[(i*384+j)*16+k];
  //       pkpvprint[(i*384+j)*16+k] = pk[(i*384+j)*16+k];
  //       // pkpvprint[i*KYBER_N+j] = e.vec[i].coeffs[j];
  //       // pkpvprint[i*KYBER_N+j] = a[1].vec[i].coeffs[j];
  //     }
  //   }
  // }


  free(pkseq);
  free(skseq);
  free(tkp);

}


void indcpa_enc(uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t m[KYBER_INDCPA_MSGBYTES*32*2],
                uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                uint8_t coins[KYBER_SYMBYTES*(32*2-1)]
                // int16_t pkpvprint[KYBER_K*KYBER_N*16],
                // int16_t vprint[KYBER_N*16]
                )
{
  unsigned int i;
  uint8_t seed[KYBER_SYMBYTES*32*2], mseq[KYBER_INDCPA_MSGBYTES*32];
  uint8_t *cseq = (uint8_t *)malloc(KYBER_INDCPA_BYTES);
  uint8_t *tc = (uint8_t *)malloc(KYBER_INDCPA_BYTES);
  uint8_t *pkseq = (uint8_t *)malloc(KYBER_INDCPA_PUBLICKEYBYTES);
  uint8_t *tpk = (uint8_t *)malloc(KYBER_INDCPA_PUBLICKEYBYTES);
  polyvec_32 sp, spseq, tpv, pkpvseq, ep, epseq, at[KYBER_K], t[KYBER_K], atseq[KYBER_K], b;
  poly_32 v, k, epp, tp, eppseq;

  // for(int i = 0; i < 3*192; i++) {
  //   for(int j = 0; j < 16; j++) {
  //     pkseq[i*32+j*2] = pk[j*3*384+i*2];
  //     pkseq[i*32+j*2+1] = pk[j*3*384+1+i*2];
  //   }
  // }

  for (i = 0; i < 32; i++) {
    memcpy(seed+2*KYBER_SYMBYTES*i, pk+KYBER_INDCPA_PUBLICKEYBYTES/32*i+KYBER_POLYVECBYTES, KYBER_SYMBYTES);
  }

  pk_formseqto32(pk, tpk, pkseq);
  unpack_pk(&pkpvseq, pkseq);
  // for (int j = 0; j < 3; j++) {
  //   for(int k = 0; k < 256; k++) {
  //     for(int p = 0; p < 32; p++) {
  //       printf("%6d,", pkpvseq.vec[j].coeffs[k*32+p]);
  //       if (p % 16 == 15) printf("\n");
  //     }
  //   }
  //   printf("\n/////////////////////////////\n");
  // }
  // printf("\n");

  msg_formseqto32(m, mseq);
  // for (int i = 0; i < 1024; i++) {
  //   printf("%5d", mseq[i]);
  //   if (i%16 == 15) printf("\n");
  // }

  poly_frommsg_32(&k, mseq);

  gen_at(at, seed);
  matrix_formseqto32(at, t, atseq);
  // for (int i = 0; i < KYBER_K; i++) {
  //   for (int j = 0; j < KYBER_K; j++) {
  //     for(int k = 0; k < KYBER_N; k++){
  //       for(int p = 0; p < 32; p++) {
  //         printf("%5d,", atseq[i].vec[j].coeffs[k*32+p]);
  //         if (p % 16 == 15) printf("\n");
  //       }
  //       // printf("\n");
  //     }
  //     printf("\n////////////////////////////////\n");
  //   }
  //   printf("\n*******************************\n");
  // }
  // printf("\n");

#ifdef KYBER_90S
#define NOISE_NBLOCKS ((KYBER_ETA1*KYBER_N/4)/AES256CTR_BLOCKBYTES) /* Assumes divisibility */
#define CIPHERTEXTNOISE_NBLOCKS ((KYBER_ETA2*KYBER_N/4)/AES256CTR_BLOCKBYTES) /* Assumes divisibility */
  uint64_t nonce = 0;
  ALIGNED_UINT8(NOISE_NBLOCKS*AES256CTR_BLOCKBYTES+32) buf; /* +32 bytes as required by poly_cbd_eta1 */
  aes256ctr_ctx state;
  aes256ctr_init(&state, coins, nonce++);
  for(i=0;i<KYBER_K;i++) {
    aes256ctr_squeezeblocks(buf.coeffs, NOISE_NBLOCKS, &state);
    state.n = _mm_loadl_epi64((__m128i *)&nonce);
    nonce += 1;
    poly_cbd_eta1(&sp.vec[i], buf.vec);
  }
  for(i=0;i<KYBER_K;i++) {
    aes256ctr_squeezeblocks(buf.coeffs, CIPHERTEXTNOISE_NBLOCKS, &state);
    state.n = _mm_loadl_epi64((__m128i *)&nonce);
    nonce += 1;
    poly_cbd_eta2(&ep.vec[i], buf.vec);
  }
  aes256ctr_squeezeblocks(buf.coeffs, CIPHERTEXTNOISE_NBLOCKS, &state);
  poly_cbd_eta2(&epp, buf.vec);
#else
#if KYBER_K == 2
  poly_getnoise_eta1122_4x(sp.vec+0, sp.vec+1, ep.vec+0, ep.vec+1, coins, 0, 1, 2, 3);
  poly_getnoise_eta2_4x(&epp, coins, 4);
#elif KYBER_K == 3
  // for (j = 0; j < KYBER_K; j++) {
  //   for(l = 0; l < KYBER_N; l++) {
  //     for(p = 0; p < 16; p++) {
  //       sp.vec[j].coeffs[l*16+p] = 21;
  //       ep.vec[j].coeffs[l*16+p] = 21;
  //       epp.coeffs[l*16+p] = 21;
  //     }
  //   }
  // }
  poly_getnoise_eta1_8x(sp.vec+0, sp.vec+1, sp.vec+2, ep.vec+0, ep.vec+1, ep.vec+2, &epp, b.vec+0, coins, 0, 1, 2 ,3, 4, 5, 6, 7);
  polyvec_formseqto32(&sp, &tpv, &spseq);
  polyvec_formseqto32(&ep, &tpv, &epseq);
  poly_formseqto32(&epp, &tp, &eppseq);
  // for (int j = 0; j < KYBER_K; j++) {
  //   for(int k = 0; k < KYBER_N; k++){
  //     for (int p = 0; p < 32; p++) {
  //       printf("%5d,", epseq.vec[j].coeffs[k*32+p]);
  //       if (p % 32 == 31) printf("\n");
  //     }
  //   }
  //   printf("////////////////////////////\n\n");
  // }

  // for(int k = 0; k < KYBER_N; k++){
  //   for (int p = 0; p < 32; p++) {
  //     printf("%5d,", eppseq.coeffs[k*32+p]);
  //     if (p % 16 == 15) printf("\n");
  //   }
  // }
  // printf("////////////////////////////\n\n");
#elif KYBER_K == 4
  poly_getnoise_eta1_4x(sp.vec+0, sp.vec+1, sp.vec+2, sp.vec+3, coins, 0, 1, 2, 3);
  poly_getnoise_eta1_4x(ep.vec+0, ep.vec+1, ep.vec+2, ep.vec+3, coins, 4, 5, 6, 7);
  poly_getnoise_eta2_4x(&epp, coins, 8);
#endif
#endif

  polyvec_ntt(&spseq);
  // for (int j = 0; j < KYBER_K; j++) {
  //   for(int k = 0; k < KYBER_N; k++){
  //     for (int p = 0; p < 32; p++) {
  //       printf("%5d,", spseq.vec[j].coeffs[k*32+p]);
  //       if (p % 16 == 15) printf("\n");
  //     }
  //   }
  //   printf("////////////////////////////\n\n");
  // }

  // matrix-vector multiplication
  for(i=0;i<KYBER_K;i++)
    polyvec_basemul_acc_montgomery(&b.vec[i], &atseq[i], &spseq);
  polyvec_basemul_acc_montgomery(&v, &pkpvseq, &spseq);

  polyvec_invntt_tomont(&b);
  poly_invntt_tomont(&v);

  polyvec_add(&b, &b, &epseq);
  poly_add(&v, &v, &eppseq);
  poly_add(&v, &v, &k);
  polyvec_reduce(&b);
  poly_reduce(&v);

  // for (int j = 0; j < 3; j++) {
  //   for(int k = 0; k < 256; k++){
  //     for (int p = 0; p < 32; p++) {
  //       printf("%6d,", b.vec[j].coeffs[k*32+p]);
  //       if (p % 16 == 15) printf("\n");
  //     }
  //   }
  //   printf("////////////////////////////\n\n");
  // }

  // for(int l = 0; l < KYBER_N; l++){
  //   for (int p = 0; p < 32; p++) {
  //     // printf("%6d,", v.coeffs[l*32+p]);
  //     printf("%6d,", k.coeffs[l*32+p]);
  //     if (p % 16 == 15) printf("\n");
  //   }
  // }
  // printf("////////////////////////////\n\n");

  pack_ciphertext(cseq, &b, &v);
  cipher_formseqfrom32(cseq, tc, c);

  // for (int i = 0; i < 128*32; i++) {
  //   printf("%5d, ", cseq[3*320*32+i]);
  //   if (i % 16 == 15) printf("\n");
  // }

  free(cseq);
  free(tc);
  free(pkseq);
  free(tpk);

}


void indcpa_dec(uint8_t m[KYBER_INDCPA_MSGBYTES*32*2],
                uint8_t c[KYBER_INDCPA_BYTES],
                uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES]
                )
{
  polyvec_32 b, skpvseq;
  poly_32 v, mp;
  uint8_t mseq[KYBER_INDCPA_MSGBYTES*32*2];
  uint8_t *cseq = (uint8_t *)malloc(KYBER_INDCPA_BYTES);
  uint8_t *tc = (uint8_t *)malloc(KYBER_INDCPA_BYTES);
  uint8_t *skseq = (uint8_t *)malloc(KYBER_INDCPA_SECRETKEYBYTES);
  uint8_t *tsk = (uint8_t *)malloc(KYBER_INDCPA_SECRETKEYBYTES);

  cipher_formseqto32(c, tc, cseq);

  // for (int i = 0; i < 128*32; i++) {
  //   printf("%5d, ", cseq[3*320*32+i]);
  //   if (i % 16 == 15) printf("\n");
  // }

  unpack_ciphertext(&b, &v, cseq);

  sk_formseqto32(sk, tsk, skseq);
  unpack_sk(&skpvseq, skseq);

  polyvec_ntt(&b);
  polyvec_basemul_acc_montgomery(&mp, &skpvseq, &b);
  poly_invntt_tomont(&mp);

  poly_sub(&mp, &v, &mp);
  poly_reduce(&mp);

  // for (int j = 0; j < KYBER_K; j++) {
  //   for(int k = 0; k < KYBER_N; k++){
  //     for (int p = 0; p < 32; p++) {
  //       printf("%6d,", b.vec[j].coeffs[k*32+p]);
  //       // printf("%6d,", skpvseq.vec[j].coeffs[k*32+p]);
  //       if (p % 16 == 15) printf("\n");
  //     }
  //   }
  //   printf("////////////////////////////\n\n");
  // }

  // for(int l = 0; l < KYBER_N; l++){
  //   for (int p = 0; p < 32; p++) {
  //     printf("%6d,", v.coeffs[l*32+p]);
  //     // printf("%6d,", mp.coeffs[l*32+p]);
  //     if (p % 16 == 15) printf("\n");
  //   }
  // }
  // printf("////////////////////////////\n\n");

  // for(int j = 0; j < KYBER_N*16; j++) {
  //    vprint[j] = mp.coeffs[j];
  // }

  poly_tomsg_32(mseq, &mp);
  msg_formseqfrom32(mseq, m);

  // for (int i = 0; i < 1024; i++) {
  //   printf("%5d", mseq[i]);
  //   if (i%16 == 15) printf("\n");
  // }

  free(cseq);
  free(tc);
  free(skseq);
  free(tsk);
  
}
