#include "consts_32.h"
.include "fq.inc"


.macro tobyte256b pos1, pos2
#load
vmovdqa32		(\pos1)*2(%rsi),%zmm5
vmovdqa32		(\pos1 +  32)*2(%rsi),%zmm6
vmovdqa32		(\pos1 +  64)*2(%rsi),%zmm7
vmovdqa32		(\pos1 +  96)*2(%rsi),%zmm8
vmovdqa32		(\pos1 + 128)*2(%rsi),%zmm9
vmovdqa32		(\pos1 + 160)*2(%rsi),%zmm10
vmovdqa32		(\pos1 + 192)*2(%rsi),%zmm11
vmovdqa32		(\pos1 + 224)*2(%rsi),%zmm12
vmovdqa32		(\pos1 + 256)*2(%rsi),%zmm13
vmovdqa32		(\pos1 + 288)*2(%rsi),%zmm14
vmovdqa32		(\pos1 + 320)*2(%rsi),%zmm15
vmovdqa32		(\pos1 + 352)*2(%rsi),%zmm16
vmovdqa32		(\pos1 + 384)*2(%rsi),%zmm17
vmovdqa32		(\pos1 + 416)*2(%rsi),%zmm18
vmovdqa32		(\pos1 + 448)*2(%rsi),%zmm19
vmovdqa32		(\pos1 + 480)*2(%rsi),%zmm20

#csubq
csubq		5,21
csubq		6,21
csubq		7,21
csubq		8,21
csubq		9,21
csubq		10,21
csubq		11,21
csubq		12,21
csubq		13,21
csubq		14,21
csubq		15,21
csubq		16,21
csubq		17,21
csubq		18,21
csubq		19,21
csubq		20,21

#bitpack
vpsllw		$12,%zmm6,%zmm4
vpord		%zmm4,%zmm5,%zmm4

vpsrlw		$4,%zmm6,%zmm5
vpsllw		$8,%zmm7,%zmm6
vpord		%zmm5,%zmm6,%zmm5

vpsrlw		$8,%zmm7,%zmm6
vpsllw		$4,%zmm8,%zmm7
vpord		%zmm6,%zmm7,%zmm6

vpsllw		$12,%zmm10,%zmm7
vpord		%zmm7,%zmm9,%zmm7

vpsrlw		$4,%zmm10,%zmm8
vpsllw		$8,%zmm11,%zmm9
vpord		%zmm8,%zmm9,%zmm8

vpsrlw		$8,%zmm11,%zmm9
vpsllw		$4,%zmm12,%zmm10
vpord		%zmm9,%zmm10,%zmm9

vpsllw		$12,%zmm14,%zmm10
vpord		%zmm10,%zmm13,%zmm10

vpsrlw		$4,%zmm14,%zmm11
vpsllw		$8,%zmm15,%zmm12
vpord		%zmm11,%zmm12,%zmm11

vpsrlw		$8,%zmm15,%zmm12
vpsllw		$4,%zmm16,%zmm13
vpord		%zmm12,%zmm13,%zmm12

vpsllw		$12,%zmm18,%zmm13
vpord		%zmm13,%zmm17,%zmm13

vpsrlw		$4,%zmm18,%zmm14
vpsllw		$8,%zmm19,%zmm15
vpord		%zmm14,%zmm15,%zmm14

vpsrlw		$8,%zmm19,%zmm15
vpsllw		$4,%zmm20,%zmm16
vpord		%zmm15,%zmm16,%zmm15

#store
vmovdqu32		%zmm4,(\pos2)*2(%rdi)
vmovdqu32		%zmm5,(\pos2 +  32)*2(%rdi)
vmovdqu32		%zmm6,(\pos2 +  64)*2(%rdi)
vmovdqu32		%zmm7,(\pos2 +  96)*2(%rdi)
vmovdqu32		%zmm8,(\pos2 + 128)*2(%rdi)
vmovdqu32		%zmm9,(\pos2 + 160)*2(%rdi)
vmovdqu32		%zmm10,(\pos2 + 192)*2(%rdi)
vmovdqu32		%zmm11,(\pos2 + 224)*2(%rdi)
vmovdqu32		%zmm12,(\pos2 + 256)*2(%rdi)
vmovdqu32		%zmm13,(\pos2 + 288)*2(%rdi)
vmovdqu32		%zmm14,(\pos2 + 320)*2(%rdi)
vmovdqu32		%zmm15,(\pos2 + 352)*2(%rdi)

.endm


#
#   @brief: concatenate signicifant bits so that every coefficients is 12-bit-significant, then don't change the sequence of coefficients(which is different from original) \
#           as the NTT in one-way AVX2 Kyber is incomplete(this word is not accurate here,不完整是说其顺序与7层NTT后的顺序不同，Kyber中NTT都是不完整的) so they need shuffle to adapt other versions like C or so on.
#   @details: 16-way doesn't use shuffle1,2,4,8 as our NTT is complete. So after pack_sk(), indcpa_keypair() can still adapt to other version Kyber.
#             
#
.global cdecl(ntttobytes_avx512_32)
cdecl(ntttobytes_avx512_32):
#consts
vmovdqa32		_32XQ*2(%rdx),%zmm0

tobyte256b 0, 0
tobyte256b 512, 384
tobyte256b 1024, 768
tobyte256b 1536, 1152
tobyte256b 2048, 1536
tobyte256b 2560, 1920
tobyte256b 3072, 2304
tobyte256b 3584, 2688
tobyte256b 4096, 3072
tobyte256b 4608, 3456

tobyte256b 5120, 3840
tobyte256b 5632, 4224
tobyte256b 6144, 4608
tobyte256b 6656, 4992
tobyte256b 7168, 5376
tobyte256b 7680, 5760

ret


.macro frombyte256b pos1, pos2
#load
vmovdqu32		(\pos1)*2(%rsi),%zmm10  //前6个寄存器这样写的顺序是为了下面的bitunpack直接套用原来的代码
vmovdqu32		(\pos1 +  32)*2(%rsi),%zmm7
vmovdqu32		(\pos1 +  64)*2(%rsi),%zmm4
vmovdqu32		(\pos1 +  96)*2(%rsi),%zmm8
vmovdqu32		(\pos1 + 128)*2(%rsi),%zmm5
vmovdqu32		(\pos1 + 160)*2(%rsi),%zmm9
vmovdqu32		(\pos1 + 192)*2(%rsi),%zmm16
vmovdqu32		(\pos1 + 224)*2(%rsi),%zmm17
vmovdqu32		(\pos1 + 256)*2(%rsi),%zmm18
vmovdqu32		(\pos1 + 288)*2(%rsi),%zmm19
vmovdqu32		(\pos1 + 320)*2(%rsi),%zmm20
vmovdqu32		(\pos1 + 352)*2(%rsi),%zmm21

#bitunpack
vpsrlw		$12,%zmm10,%zmm11
vpsllw		$4,%zmm7,%zmm12
vpord		%zmm11,%zmm12,%zmm11
vpandd		%zmm0,%zmm10,%zmm10
vpandd		%zmm0,%zmm11,%zmm11

vpsrlw		$8,%zmm7,%zmm12
vpsllw		$8,%zmm4,%zmm13
vpord		%zmm12,%zmm13,%zmm12
vpandd		%zmm0,%zmm12,%zmm12

vpsrlw		$4,%zmm4,%zmm13
vpandd		%zmm0,%zmm13,%zmm13

vpsrlw		$12,%zmm8,%zmm14
vpsllw		$4,%zmm5,%zmm15
vpord		%zmm14,%zmm15,%zmm14
vpandd		%zmm0,%zmm8,%zmm8
vpandd		%zmm0,%zmm14,%zmm14

vpsrlw		$8,%zmm5,%zmm15
vpsllw		$8,%zmm9,%zmm1
vpord		%zmm15,%zmm1,%zmm15
vpandd		%zmm0,%zmm15,%zmm15

vpsrlw		$4,%zmm9,%zmm1
vpandd		%zmm0,%zmm1,%zmm1

///////////////////////////////////

vpsrlw		$12,%zmm16,%zmm22
vpsllw		$4,%zmm17,%zmm23
vpord		%zmm22,%zmm23,%zmm22
vpandd		%zmm0,%zmm16,%zmm16
vpandd		%zmm0,%zmm22,%zmm22

vpsrlw		$8,%zmm17,%zmm23
vpsllw		$8,%zmm18,%zmm24
vpord		%zmm23,%zmm24,%zmm23
vpandd		%zmm0,%zmm23,%zmm23

vpsrlw		$4,%zmm18,%zmm24
vpandd		%zmm0,%zmm24,%zmm24

vpsrlw		$12,%zmm19,%zmm25
vpsllw		$4,%zmm20,%zmm26
vpord		%zmm25,%zmm26,%zmm25
vpandd		%zmm0,%zmm19,%zmm19
vpandd		%zmm0,%zmm25,%zmm25

vpsrlw		$8,%zmm20,%zmm26
vpsllw		$8,%zmm21,%zmm27
vpord		%zmm26,%zmm27,%zmm26
vpandd		%zmm0,%zmm26,%zmm26

vpsrlw		$4,%zmm21,%zmm21
vpandd		%zmm0,%zmm21,%zmm21


#store
vmovdqa32		%zmm10,(\pos2)*2(%rdi)
vmovdqa32		%zmm11,(\pos2 +  32)*2(%rdi)
vmovdqa32		%zmm12,(\pos2 +  64)*2(%rdi)
vmovdqa32		%zmm13,(\pos2 +  96)*2(%rdi)
vmovdqa32		%zmm8,(\pos2 + 128)*2(%rdi)
vmovdqa32		%zmm14,(\pos2 + 160)*2(%rdi)
vmovdqa32		%zmm15,(\pos2 + 192)*2(%rdi)
vmovdqa32		%zmm1,(\pos2 + 224)*2(%rdi)

vmovdqa32		%zmm16,(\pos2 + 256)*2(%rdi)
vmovdqa32		%zmm22,(\pos2 + 288)*2(%rdi)
vmovdqa32		%zmm23,(\pos2 + 320)*2(%rdi)
vmovdqa32		%zmm24,(\pos2 + 352)*2(%rdi)
vmovdqa32		%zmm19,(\pos2 + 384)*2(%rdi)
vmovdqa32		%zmm25,(\pos2 + 416)*2(%rdi)
vmovdqa32		%zmm26,(\pos2 + 448)*2(%rdi)
vmovdqa32		%zmm21,(\pos2 + 480)*2(%rdi)

.endm

#
#   @brief: seperate the concatenation between coefficients, and this function is the inverse of ntttobytes_avx_32
#
#
.global cdecl(nttfrombytes_avx512_32)
cdecl(nttfrombytes_avx512_32):
#consts
vmovdqa32		_32XMASK*2(%rdx),%zmm0

frombyte256b 0, 0
frombyte256b 384, 512
frombyte256b 768, 1024
frombyte256b 1152, 1536
frombyte256b 1536, 2048
frombyte256b 1920, 2560
frombyte256b 2304, 3072
frombyte256b 2688, 3584
frombyte256b 3072, 4096
frombyte256b 3456, 4608

frombyte256b 3840, 5120
frombyte256b 4224, 5632
frombyte256b 4608, 6144
frombyte256b 4992, 6656
frombyte256b 5376, 7168
frombyte256b 5760, 7680

ret