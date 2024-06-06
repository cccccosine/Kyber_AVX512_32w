#ifndef PARAMS_32_H
#define PARAMS_32_H

#ifndef KYBER_K
#define KYBER_K 3	/* Change this for different security strengths */
#endif

//#define KYBER_90S	/* Uncomment this if you want the 90S variant */

/* Don't change parameters below this line */
#if   (KYBER_K == 2)
#ifdef KYBER_90S
#define KYBER_NAMESPACE(s) pqcrystals_kyber512_90s_avx2_##s
#else
#define KYBER_NAMESPACE(s) pqcrystals_kyber512_avx2_##s
#endif
#elif (KYBER_K == 3)
#ifdef KYBER_90S
#define KYBER_NAMESPACE(s) pqcrystals_kyber768_90s_avx2_##s
#else
#define KYBER_NAMESPACE(s) pqcrystals_kyber768_avx2_##s
#endif
#elif (KYBER_K == 4)
#ifdef KYBER_90S
#define KYBER_NAMESPACE(s) pqcrystals_kyber1024_90s_avx2_##s
#else
#define KYBER_NAMESPACE(s) pqcrystals_kyber1024_avx2_##s
#endif
#else
#error "KYBER_K must be in {2,3,4}"
#endif

#define KYBER_N 256
#define KYBER_Q 3329

#define KYBER_SYMBYTES 32   /* size in bytes of hashes, and seeds */
#define KYBER_SSBYTES  32   /* size in bytes of shared key */

#define KYBER_POLYBYTES		384
#define KYBER_POLYVECBYTES	(KYBER_K * KYBER_POLYBYTES)  //3*384

#if KYBER_K == 2
#define KYBER_ETA1 3
#define KYBER_POLYCOMPRESSEDBYTES    128 * 16
#define KYBER_POLYVECCOMPRESSEDBYTES (KYBER_K * 320 * 16)
#elif KYBER_K == 3
#define KYBER_ETA1 2
#define KYBER_POLYCOMPRESSEDBYTES    128 * 32
#define KYBER_POLYVECCOMPRESSEDBYTES (KYBER_K * 320 * 32)  //3*320*32
#elif KYBER_K == 4
#define KYBER_ETA1 2
#define KYBER_POLYCOMPRESSEDBYTES    160 * 16
#define KYBER_POLYVECCOMPRESSEDBYTES (KYBER_K * 352 * 16)
#endif

#define KYBER_ETA2 2

#define KYBER_INDCPA_MSGBYTES       (KYBER_SYMBYTES)  //32
#define KYBER_INDCPA_PUBLICKEYBYTES (KYBER_POLYVECBYTES + KYBER_SYMBYTES)*32  //(3*384+32)*32, 现在直接在formseqfrom32时将每way的pk后面预留publicseed的位置，所以长度就不需要多一个32byte了
#define KYBER_INDCPA_SECRETKEYBYTES (KYBER_POLYVECBYTES*32)  //3*384*32
#define KYBER_INDCPA_BYTES          (KYBER_POLYVECCOMPRESSEDBYTES + KYBER_POLYCOMPRESSEDBYTES)  //3*320*32 + 128*32  应该是(3*320+128)*32

#define KYBER_PUBLICKEYBYTES  (KYBER_INDCPA_PUBLICKEYBYTES)  //(3*384+32)*32
/* 32 bytes of additional space to save H(pk) */
#define KYBER_SECRETKEYBYTES  (KYBER_INDCPA_SECRETKEYBYTES + KYBER_INDCPA_PUBLICKEYBYTES + 2*KYBER_SYMBYTES*32)  //3*384*32 + (3*384+32)*32 + 2*32*32
#define KYBER_CIPHERTEXTBYTES (KYBER_INDCPA_BYTES)  //3*320*32 + 128*32

#endif
