#ifndef FIPS202X8_H
#define FIPS202X8_H

#include <stddef.h>
#include <stdint.h>
#include <immintrin.h>

#define FIPS202X8_NAMESPACE(s) pqcrystals_kyber_fips202x8_avx512_##s

typedef struct {
  __m512i s[25];
} keccakx8_state;

#define shake128x8_absorb_once FIPS202X8_NAMESPACE(shake128x8_absorb_once)
void shake128x8_absorb_once(keccakx8_state *state,
                            const uint8_t *in0,
                            const uint8_t *in1,
                            const uint8_t *in2,
                            const uint8_t *in3,
                            const uint8_t *in4,
                            const uint8_t *in5,
                            const uint8_t *in6,
                            const uint8_t *in7,
                            size_t inlen);


#define shake128x8_squeezeblocks FIPS202X8_NAMESPACE(shake128x8_squeezeblocks)
void shake128x8_squeezeblocks(uint8_t *out0,
                              uint8_t *out1,
                              uint8_t *out2,
                              uint8_t *out3,
                              uint8_t *out4,
                              uint8_t *out5,
                              uint8_t *out6,
                              uint8_t *out7,
                              size_t nblocks,
                              keccakx8_state *state);

#define shake256x8_absorb_once FIPS202X8_NAMESPACE(shake256x8_absorb_once)
void shake256x8_absorb_once(keccakx8_state *state,
                            const uint8_t *in0,
                            const uint8_t *in1,
                            const uint8_t *in2,
                            const uint8_t *in3,
                            const uint8_t *in4,
                            const uint8_t *in5,
                            const uint8_t *in6,
                            const uint8_t *in7,
                            size_t inlen);

#define shake256x8_squeezeblocks FIPS202X8_NAMESPACE(shake256x8_squeezeblocks)
void shake256x8_squeezeblocks(uint8_t *out0,
                              uint8_t *out1,
                              uint8_t *out2,
                              uint8_t *out3,
                              uint8_t *out4,
                              uint8_t *out5,
                              uint8_t *out6,
                              uint8_t *out7,
                              size_t nblocks,
                              keccakx8_state *state);

#define shake128x8 FIPS202X8_NAMESPACE(shake128x8)
void shake128x8(uint8_t *out0,
                uint8_t *out1,
                uint8_t *out2,
                uint8_t *out3,
                uint8_t *out4,
                uint8_t *out5,
                uint8_t *out6,
                uint8_t *out7,
                size_t outlen,
                const uint8_t *in0,
                const uint8_t *in1,
                const uint8_t *in2,
                const uint8_t *in3,
                const uint8_t *in4,
                const uint8_t *in5,
                const uint8_t *in6,
                const uint8_t *in7,
                size_t inlen);

#define shake256x8 FIPS202X8_NAMESPACE(shake256x8)
void shake256x8(uint8_t *out0,
                uint8_t *out1,
                uint8_t *out2,
                uint8_t *out3,
                uint8_t *out4,
                uint8_t *out5,
                uint8_t *out6,
                uint8_t *out7,
                size_t outlen,
                const uint8_t *in0,
                const uint8_t *in1,
                const uint8_t *in2,
                const uint8_t *in3,
                const uint8_t *in4,
                const uint8_t *in5,
                const uint8_t *in6,
                const uint8_t *in7,
                size_t inlen);

#endif

#define sha3x8_256 FIPS202X8_NAMESPACE(sha3x8_256)
void sha3x8_256(uint8_t *out0,
                uint8_t *out1,
                uint8_t *out2,
                uint8_t *out3,
                uint8_t *out4,
                uint8_t *out5,
                uint8_t *out6,
                uint8_t *out7,
                const uint8_t *in0,
                const uint8_t *in1,
                const uint8_t *in2,
                const uint8_t *in3,
                const uint8_t *in4,
                const uint8_t *in5,
                const uint8_t *in6,
                const uint8_t *in7,
                size_t inlen
               );

#define sha3x8_512 FIPS202X8_NAMESPACE(sha3x8_512)
void sha3x8_512(uint8_t *out0,
                uint8_t *out1,
                uint8_t *out2,
                uint8_t *out3,
                uint8_t *out4,
                uint8_t *out5,
                uint8_t *out6,
                uint8_t *out7,
                const uint8_t *in0,
                const uint8_t *in1,
                const uint8_t *in2,
                const uint8_t *in3,
                const uint8_t *in4,
                const uint8_t *in5,
                const uint8_t *in6,
                const uint8_t *in7,
                size_t inlen
               );
