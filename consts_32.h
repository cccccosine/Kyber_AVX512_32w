#ifndef CONSTS_32_H
#define CONSTS_32_H

#include "params.h"

#define _32XQ              0
#define _ZETAS_EXP_32      32
#define _32XQINV           8160
#define _ZETAS_BASEMUL     8192
#define _32XFLO            12288
#define _32XFHI            12320
#define _32XV              12352
#define _32XMONTSQLO       12384
#define _32XMONTSQHI       12416
#define _32XMASK           12448
#define SEQ_SHUFIDX_32     12480
#define SEQ_PERMIDX_32     12512
#define SEQ_INVSHUFIDX_32  12544
#define SEQ_INVPERMIDX_32  12576
#define K1IDX              12608
#define K2IDX              12609

/* The C ABI on MacOS exports all symbols with a leading
 * underscore. This means that any symbols we refer to from
 * C files (functions) can't be found, and all symbols we
 * refer to from ASM also can't be found.
 *
 * This define helps us get around this
 */
#ifdef __ASSEMBLER__
#if defined(__WIN32__) || defined(__APPLE__)
#define decorate(s) _##s
#define cdecl2(s) decorate(s)
#define cdecl(s) cdecl2(KYBER_NAMESPACE(##s))
#else
#define cdecl(s) KYBER_NAMESPACE(##s)
#endif
#endif

#ifndef __ASSEMBLER__
#include "align.h"
typedef ALIGNED_INT16_512(12620) qdata_t_32;     //大小未定
#define qdata_32 KYBER_NAMESPACE(qdata_32)
extern const qdata_t_32 qdata_32;
#endif

#endif