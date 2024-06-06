#ifndef REDUCE_32_H
#define REDUCE_32_H

#include "params.h"
#include <immintrin.h>

#define reduce_avx_32 KYBER_NAMESPACE(reduce_avx_32)
void reduce_avx_32(__m512i *r, const __m512i *qdata_32);
#define tomont_avx_32 KYBER_NAMESPACE(tomont_avx_32)
void tomont_avx_32(__m512i *r, const __m512i *qdata_32);

#endif
