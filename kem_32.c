#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include "align.h"
#include "params.h"
#include "kem_32.h"
#include "indcpa_32.h"
#include "verify_32.h"
#include "symmetric.h"
#include "randombytes.h"
#include "rejsample.h"


int crypto_kem_keypair(uint8_t *pk,
                       uint8_t *sk)
{
  uint8_t buf[KYBER_SYMBYTES*32];
  uint8_t *sk_32 = (uint8_t *)malloc(KYBER_INDCPA_SECRETKEYBYTES);

  indcpa_keypair(pk, sk_32);
  // 32个分离的sk||pk||publicseed
  for(int i = 0; i < 32; i++) {
    memcpy(sk+(KYBER_SECRETKEYBYTES/32)*i, sk_32+KYBER_POLYVECBYTES*i, KYBER_POLYVECBYTES);
    memcpy(sk+(KYBER_SECRETKEYBYTES/32)*i+KYBER_POLYVECBYTES, pk+KYBER_INDCPA_PUBLICKEYBYTES/32*i, KYBER_INDCPA_PUBLICKEYBYTES/32);
  }
  
  for(int i = 0; i < 4; i++) {
    hash_hx8(sk+(8*i)*(KYBER_SECRETKEYBYTES/32)+2*KYBER_POLYVECBYTES+KYBER_SYMBYTES,
             sk+(8*i+1)*(KYBER_SECRETKEYBYTES/32)+2*KYBER_POLYVECBYTES+KYBER_SYMBYTES,
             sk+(8*i+2)*(KYBER_SECRETKEYBYTES/32)+2*KYBER_POLYVECBYTES+KYBER_SYMBYTES,
             sk+(8*i+3)*(KYBER_SECRETKEYBYTES/32)+2*KYBER_POLYVECBYTES+KYBER_SYMBYTES,
             sk+(8*i+4)*(KYBER_SECRETKEYBYTES/32)+2*KYBER_POLYVECBYTES+KYBER_SYMBYTES,
             sk+(8*i+5)*(KYBER_SECRETKEYBYTES/32)+2*KYBER_POLYVECBYTES+KYBER_SYMBYTES,
             sk+(8*i+6)*(KYBER_SECRETKEYBYTES/32)+2*KYBER_POLYVECBYTES+KYBER_SYMBYTES,
             sk+(8*i+7)*(KYBER_SECRETKEYBYTES/32)+2*KYBER_POLYVECBYTES+KYBER_SYMBYTES,
             sk+(8*i)*(KYBER_SECRETKEYBYTES/32)+KYBER_POLYVECBYTES,
             sk+(8*i+1)*(KYBER_SECRETKEYBYTES/32)+KYBER_POLYVECBYTES,
             sk+(8*i+2)*(KYBER_SECRETKEYBYTES/32)+KYBER_POLYVECBYTES,
             sk+(8*i+3)*(KYBER_SECRETKEYBYTES/32)+KYBER_POLYVECBYTES,
             sk+(8*i+4)*(KYBER_SECRETKEYBYTES/32)+KYBER_POLYVECBYTES,
             sk+(8*i+5)*(KYBER_SECRETKEYBYTES/32)+KYBER_POLYVECBYTES,
             sk+(8*i+6)*(KYBER_SECRETKEYBYTES/32)+KYBER_POLYVECBYTES,
             sk+(8*i+7)*(KYBER_SECRETKEYBYTES/32)+KYBER_POLYVECBYTES,
             KYBER_INDCPA_PUBLICKEYBYTES/32);
  }

  // /* Value z for pseudo-random output on reject */
  randombytes(buf, KYBER_SYMBYTES*32);

  for(int i = 0; i < 32; i++) {
    memcpy(sk+(i+1)*(KYBER_SECRETKEYBYTES/32)-KYBER_SYMBYTES, buf+i*KYBER_SYMBYTES, KYBER_SYMBYTES);
  }

  free(sk_32);

  return 0;
}


int crypto_kem_enc(uint8_t *ct,
                   uint8_t *ss,
                   uint8_t *pk)
{
  uint8_t buf[KYBER_INDCPA_MSGBYTES*32*2];
  /* Will contain key, coins */
  uint8_t kr[2*KYBER_SYMBYTES*32];

  randombytes(buf + KYBER_SYMBYTES*32, KYBER_SYMBYTES*32);

  /* Don't release system RNG output */
  //产生16个H(m)，并且留出下一步进行单路连接H(pk)需要的空间
  for(int i = 0; i < 4; i++) {
    hash_hx8(buf+16*i*KYBER_SYMBYTES, 
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

  /* Multitarget countermeasure for coins + contributory KEM */
  for(int i = 0; i < 4; i++) {
    hash_hx8(buf+(16*i+1)*KYBER_SYMBYTES, 
             buf+(16*i+3)*KYBER_SYMBYTES, 
             buf+(16*i+5)*KYBER_SYMBYTES, 
             buf+(16*i+7)*KYBER_SYMBYTES,
             buf+(16*i+9)*KYBER_SYMBYTES, 
             buf+(16*i+11)*KYBER_SYMBYTES, 
             buf+(16*i+13)*KYBER_SYMBYTES, 
             buf+(16*i+15)*KYBER_SYMBYTES,
             pk+KYBER_PUBLICKEYBYTES/32*i*8, 
             pk+KYBER_PUBLICKEYBYTES/32*(i*8+1), 
             pk+KYBER_PUBLICKEYBYTES/32*(i*8+2), 
             pk+KYBER_PUBLICKEYBYTES/32*(i*8+3), 
             pk+KYBER_PUBLICKEYBYTES/32*(i*8+4), 
             pk+KYBER_PUBLICKEYBYTES/32*(i*8+5), 
             pk+KYBER_PUBLICKEYBYTES/32*(i*8+6), 
             pk+KYBER_PUBLICKEYBYTES/32*(i*8+7),
             KYBER_PUBLICKEYBYTES/32);
  }

  //buf = (m||H(pk)) * 32
  //kr = (K||r) * 32 = G(m||H(pk))
  for(int i = 0; i < 4; i++) {
    hash_gx8(kr+16*i*KYBER_SYMBYTES, 
             kr+(16*i+2)*KYBER_SYMBYTES, 
             kr+(16*i+4)*KYBER_SYMBYTES, 
             kr+(16*i+6)*KYBER_SYMBYTES, 
             kr+(16*i+8)*KYBER_SYMBYTES, 
             kr+(16*i+10)*KYBER_SYMBYTES, 
             kr+(16*i+12)*KYBER_SYMBYTES, 
             kr+(16*i+14)*KYBER_SYMBYTES,
             buf+16*i*KYBER_SYMBYTES, 
             buf+(16*i+2)*KYBER_SYMBYTES, 
             buf+(16*i+4)*KYBER_SYMBYTES, 
             buf+(16*i+6)*KYBER_SYMBYTES,
             buf+(16*i+8)*KYBER_SYMBYTES, 
             buf+(16*i+10)*KYBER_SYMBYTES, 
             buf+(16*i+12)*KYBER_SYMBYTES, 
             buf+(16*i+14)*KYBER_SYMBYTES, 
             KYBER_SYMBYTES*2);
  }

  /* coins are in kr+KYBER_SYMBYTES */
  indcpa_enc(ct, buf, pk, kr+KYBER_SYMBYTES);

  /* overwrite coins in kr with H(c) */
  for(int i = 0; i < 4; i++) {
    hash_hx8(kr+(16*i+1)*KYBER_SYMBYTES, 
             kr+(16*i+3)*KYBER_SYMBYTES, 
             kr+(16*i+5)*KYBER_SYMBYTES, 
             kr+(16*i+7)*KYBER_SYMBYTES,
             kr+(16*i+9)*KYBER_SYMBYTES, 
             kr+(16*i+11)*KYBER_SYMBYTES, 
             kr+(16*i+13)*KYBER_SYMBYTES, 
             kr+(16*i+15)*KYBER_SYMBYTES, 
             ct+KYBER_CIPHERTEXTBYTES/32*i*8, 
             ct+KYBER_CIPHERTEXTBYTES/32*(i*8+1), 
             ct+KYBER_CIPHERTEXTBYTES/32*(i*8+2), 
             ct+KYBER_CIPHERTEXTBYTES/32*(i*8+3),
             ct+KYBER_CIPHERTEXTBYTES/32*(i*8+4), 
             ct+KYBER_CIPHERTEXTBYTES/32*(i*8+5), 
             ct+KYBER_CIPHERTEXTBYTES/32*(i*8+6), 
             ct+KYBER_CIPHERTEXTBYTES/32*(i*8+7),
             KYBER_CIPHERTEXTBYTES/32);
  }
  
  /* hash concatenation of pre-k and H(c) to k */
  for(int i = 0; i < 4; i++) {
    kdfx8(ss + 8*i*KYBER_SSBYTES, 
          ss + (8*i+1)*KYBER_SSBYTES, 
          ss + (8*i+2)*KYBER_SSBYTES, 
          ss + (8*i+3)*KYBER_SSBYTES, 
          ss + (8*i+4)*KYBER_SSBYTES, 
          ss + (8*i+5)*KYBER_SSBYTES, 
          ss + (8*i+6)*KYBER_SSBYTES, 
          ss + (8*i+7)*KYBER_SSBYTES,
          kr+16*i*KYBER_SYMBYTES, 
          kr+(16*i+2)*KYBER_SYMBYTES, 
          kr+(16*i+4)*KYBER_SYMBYTES, 
          kr+(16*i+6)*KYBER_SYMBYTES, 
          kr+(16*i+8)*KYBER_SYMBYTES, 
          kr+(16*i+10)*KYBER_SYMBYTES, 
          kr+(16*i+12)*KYBER_SYMBYTES, 
          kr+(16*i+14)*KYBER_SYMBYTES,
          2*KYBER_SYMBYTES);
  }
  
  return 0;
}


int crypto_kem_dec(uint8_t *ss,
                   uint8_t *ct,
                   const uint8_t *sk)
{
  //现在传进来的sk参数是16个分离的sk||pk||publicseed||H(pk)||random
  int fail;
  uint8_t buf[32*KYBER_SYMBYTES*2];
  /* Will contain key, coins */
  uint8_t kr[2*32*KYBER_SYMBYTES];
  uint8_t *sk_sepa_32 = (uint8_t *)malloc(KYBER_INDCPA_SECRETKEYBYTES);
  uint8_t *pk = (uint8_t *)malloc(KYBER_INDCPA_PUBLICKEYBYTES);
  // uint8_t skseq[KYBER_INDCPA_SECRETKEYBYTES], sk_sepa_16[KYBER_INDCPA_SECRETKEYBYTES], pk[KYBER_INDCPA_PUBLICKEYBYTES];  //sk_sepa_16和pk变量是将sk从整体sk中分离出来
  ALIGNED_UINT8(KYBER_CIPHERTEXTBYTES) cmp;

  for(int i = 0; i < 32; i++) {
    memcpy(sk_sepa_32+i*KYBER_INDCPA_SECRETKEYBYTES/32, sk+KYBER_SECRETKEYBYTES*i/32, KYBER_POLYVECBYTES);
    memcpy(pk+i*KYBER_INDCPA_PUBLICKEYBYTES/32, sk+KYBER_SECRETKEYBYTES*i/32+KYBER_POLYVECBYTES, KYBER_INDCPA_PUBLICKEYBYTES/32);
  }
  indcpa_dec(buf, ct, sk_sepa_32);

  /* Multitarget countermeasure for coins + contributory KEM */
  for(int i = 0; i < 32; i++) {
    memcpy(buf+KYBER_SYMBYTES*(2*i+1), sk+KYBER_SECRETKEYBYTES*i/32+2*KYBER_POLYVECBYTES+KYBER_SYMBYTES, KYBER_SYMBYTES);
  }
  
  //After the step above, buf = (m||H(pk)) * 32
  for(int i = 0; i < 4; i++) {
    hash_gx8(kr+16*i*KYBER_SYMBYTES, 
             kr+(16*i+2)*KYBER_SYMBYTES, 
             kr+(16*i+4)*KYBER_SYMBYTES, 
             kr+(16*i+6)*KYBER_SYMBYTES, 
             kr+(16*i+8)*KYBER_SYMBYTES, 
             kr+(16*i+10)*KYBER_SYMBYTES, 
             kr+(16*i+12)*KYBER_SYMBYTES, 
             kr+(16*i+14)*KYBER_SYMBYTES, 
             buf+16*i*KYBER_SYMBYTES, 
             buf+(16*i+2)*KYBER_SYMBYTES, 
             buf+(16*i+4)*KYBER_SYMBYTES, 
             buf+(16*i+6)*KYBER_SYMBYTES, 
             buf+(16*i+8)*KYBER_SYMBYTES, 
             buf+(16*i+10)*KYBER_SYMBYTES, 
             buf+(16*i+12)*KYBER_SYMBYTES, 
             buf+(16*i+14)*KYBER_SYMBYTES, 
             KYBER_SYMBYTES*2);
  }
  //After the step above, kr = (K||r) * 32

  indcpa_enc(cmp.coeffs, buf, pk, kr+KYBER_SYMBYTES);

  fail = verify(ct, cmp.coeffs, KYBER_CIPHERTEXTBYTES);

  /* overwrite coins in kr with H(c) */
  for(int i = 0; i < 4; i++) {
    hash_hx8(kr+(16*i+1)*KYBER_SYMBYTES, 
             kr+(16*i+3)*KYBER_SYMBYTES, 
             kr+(16*i+5)*KYBER_SYMBYTES, 
             kr+(16*i+7)*KYBER_SYMBYTES, 
             kr+(16*i+9)*KYBER_SYMBYTES, 
             kr+(16*i+11)*KYBER_SYMBYTES, 
             kr+(16*i+13)*KYBER_SYMBYTES, 
             kr+(16*i+15)*KYBER_SYMBYTES, 
             ct+KYBER_CIPHERTEXTBYTES/32*i*8, 
             ct+KYBER_CIPHERTEXTBYTES/32*(i*8+1), 
             ct+KYBER_CIPHERTEXTBYTES/32*(i*8+2), 
             ct+KYBER_CIPHERTEXTBYTES/32*(i*8+3), 
             ct+KYBER_CIPHERTEXTBYTES/32*(i*8+4), 
             ct+KYBER_CIPHERTEXTBYTES/32*(i*8+5), 
             ct+KYBER_CIPHERTEXTBYTES/32*(i*8+6), 
             ct+KYBER_CIPHERTEXTBYTES/32*(i*8+7),
             KYBER_CIPHERTEXTBYTES/32);
  }

  /* Overwrite pre-k with z on re-encryption failure */
  for(int i = 0; i < 32; i++) {
    cmov(kr+KYBER_SYMBYTES*2*i, sk+KYBER_SECRETKEYBYTES/32*i-KYBER_SYMBYTES, KYBER_SYMBYTES, fail);
  }
  
  /* hash concatenation of pre-k and H(c) to k */
  for(int i = 0; i < 4; i++) {
    kdfx8(ss + 8*i*KYBER_SSBYTES, 
          ss + (8*i+1)*KYBER_SSBYTES, 
          ss + (8*i+2)*KYBER_SSBYTES, 
          ss + (8*i+3)*KYBER_SSBYTES, 
          ss + (8*i+4)*KYBER_SSBYTES, 
          ss + (8*i+5)*KYBER_SSBYTES, 
          ss + (8*i+6)*KYBER_SSBYTES, 
          ss + (8*i+7)*KYBER_SSBYTES, 
          kr+16*i*KYBER_SYMBYTES, 
          kr+(16*i+2)*KYBER_SYMBYTES, 
          kr+(16*i+4)*KYBER_SYMBYTES, 
          kr+(16*i+6)*KYBER_SYMBYTES, 
          kr+(16*i+8)*KYBER_SYMBYTES, 
          kr+(16*i+10)*KYBER_SYMBYTES, 
          kr+(16*i+12)*KYBER_SYMBYTES, 
          kr+(16*i+14)*KYBER_SYMBYTES, 
          2*KYBER_SYMBYTES);
  }

  free(sk_sepa_32);
  free(pk);

  return 0;
}
