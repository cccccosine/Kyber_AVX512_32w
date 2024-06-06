#include "consts_32.h"

.macro schoolbook off
vmovdqa32		_32XQINV*2(%rcx),%zmm0
vmovdqa32		(256*\off+  0)*2(%rsi),%zmm1		// a0   
vmovdqa32       (256*\off+ 32)*2(%rsi),%zmm2		// b0
vmovdqa32		(256*\off+ 64)*2(%rsi),%zmm3		// a1
vmovdqa32		(256*\off+ 96)*2(%rsi),%zmm4		// b1
vmovdqa32		(256*\off+128)*2(%rsi),%zmm5		// a2   
vmovdqa32       (256*\off+160)*2(%rsi),%zmm6		// b2
vmovdqa32		(256*\off+192)*2(%rsi),%zmm7		// a3
vmovdqa32		(256*\off+224)*2(%rsi),%zmm8		// b3

vpmullw		%zmm0,%zmm1,%zmm17			// a0.lo
vpmullw		%zmm0,%zmm2,%zmm18			// b0.lo
vpmullw		%zmm0,%zmm3,%zmm19			// a1.lo
vpmullw		%zmm0,%zmm4,%zmm20			// b1.lo
vpmullw		%zmm0,%zmm5,%zmm21			// a2.lo
vpmullw		%zmm0,%zmm6,%zmm22			// b2.lo
vpmullw		%zmm0,%zmm7,%zmm23			// a3.lo
vpmullw		%zmm0,%zmm8,%zmm24			// b3.lo

vmovdqa32		(256*\off+  0)*2(%rdx),%zmm9		// c0
vmovdqa32		(256*\off+ 32)*2(%rdx),%zmm10		// d0
vmovdqa32		(256*\off+ 64)*2(%rdx),%zmm11		// c1
vmovdqa32		(256*\off+ 96)*2(%rdx),%zmm12		// d1

vpmulhw		%zmm9,%zmm1,%zmm25			// a0c0.hi
vpmulhw		%zmm10,%zmm1,%zmm1			// a0d0.hi
vpmulhw		%zmm9,%zmm2,%zmm26			// b0c0.hi
vpmulhw		%zmm10,%zmm2,%zmm2			// b0d0.hi
vpmulhw		%zmm11,%zmm3,%zmm27			// a1c1.hi
vpmulhw		%zmm12,%zmm3,%zmm3			// a1d1.hi
vpmulhw		%zmm11,%zmm4,%zmm28			// b1c1.hi
vpmulhw		%zmm12,%zmm4,%zmm4			// b1d1.hi

vmovdqa32		(256*\off+128)*2(%rdx),%zmm13		// c2
vmovdqa32		(256*\off+160)*2(%rdx),%zmm14		// d2
vmovdqa32		(256*\off+192)*2(%rdx),%zmm15		// c3
vmovdqa32		(256*\off+224)*2(%rdx),%zmm16		// d3

vpmulhw		%zmm13,%zmm5,%zmm29			// a2c2.hi
vpmulhw		%zmm14,%zmm5,%zmm5			// a2d2.hi
vpmulhw		%zmm13,%zmm6,%zmm30			// b2c2.hi
vpmulhw		%zmm14,%zmm6,%zmm6			// b2d2.hi
vpmulhw		%zmm15,%zmm7,%zmm31			// a3c3.hi
vpmulhw		%zmm16,%zmm7,%zmm7			// a3d3.hi
vpmulhw		%zmm15,%zmm8,%zmm0			// b3c3.hi
vpmulhw		%zmm16,%zmm8,%zmm8			// b3d3.hi

vmovdqa32		%zmm25,(%rsp)
/**/
vpmullw		%zmm9,%zmm17,%zmm25			// a0c0.lo
vpmullw		%zmm10,%zmm17,%zmm17		// a0d0.lo
vpmullw		%zmm9,%zmm18,%zmm9			// b0c0.lo
vpmullw		%zmm10,%zmm18,%zmm18		// b0d0.lo
vpmullw		%zmm11,%zmm19,%zmm10		// a1c1.lo
vpmullw		%zmm12,%zmm19,%zmm19		// a1d1.lo
vpmullw		%zmm11,%zmm20,%zmm11		// b1c1.lo
vpmullw		%zmm12,%zmm20,%zmm20		// b1d1.lo

vpmullw		%zmm13,%zmm21,%zmm12		// a2c2.lo
vpmullw		%zmm14,%zmm21,%zmm21		// a2d2.lo
vpmullw		%zmm13,%zmm22,%zmm13		// b2c2.lo
vpmullw		%zmm14,%zmm22,%zmm22		// b2d2.lo
vpmullw		%zmm15,%zmm23,%zmm14		// a3c3.lo
vpmullw		%zmm16,%zmm23,%zmm23		// a3d3.lo
vpmullw		%zmm15,%zmm24,%zmm15		// b3c3.lo
vpmullw		%zmm16,%zmm24,%zmm24		// b3d3.lo

vmovdqa32		_32XQ*2(%rcx),%zmm16
vpmulhw		%zmm16,%zmm25,%zmm25
vpmulhw		%zmm16,%zmm17,%zmm17
vpmulhw		%zmm16,%zmm9,%zmm9
vpmulhw		%zmm16,%zmm18,%zmm18
vpmulhw		%zmm16,%zmm10,%zmm10
vpmulhw		%zmm16,%zmm19,%zmm19
vpmulhw		%zmm16,%zmm11,%zmm11
vpmulhw		%zmm16,%zmm20,%zmm20
vpmulhw		%zmm16,%zmm12,%zmm12
vpmulhw		%zmm16,%zmm21,%zmm21
vpmulhw		%zmm16,%zmm13,%zmm13
vpmulhw		%zmm16,%zmm22,%zmm22
vpmulhw		%zmm16,%zmm14,%zmm14
vpmulhw		%zmm16,%zmm23,%zmm23
vpmulhw		%zmm16,%zmm15,%zmm15
vpmulhw		%zmm16,%zmm24,%zmm24

vpsubw		(%rsp),%zmm25,%zmm25		// -a0c0
vpsubw		%zmm17,%zmm1,%zmm17			// a0d0
vpsubw		%zmm9,%zmm26,%zmm9			// b0c0
vpsubw		%zmm18,%zmm2,%zmm18			// b0d0
vpsubw		%zmm10,%zmm27,%zmm10		// a1c1
vpsubw		%zmm19,%zmm3,%zmm19			// a1d1
vpsubw		%zmm11,%zmm28,%zmm11		// b1c1
vpsubw		%zmm20,%zmm4,%zmm20			// b1d1

vpsubw		%zmm29,%zmm12,%zmm12		// -a2c2
vpsubw		%zmm21,%zmm5,%zmm21			// a2d2
vpsubw		%zmm13,%zmm30,%zmm13		// b2c2
vpsubw		%zmm22,%zmm6,%zmm22			// b2d2
vpsubw		%zmm14,%zmm31,%zmm14		// a3c3
vpsubw		%zmm23,%zmm7,%zmm23			// a3d3
vpsubw		%zmm15,%zmm0,%zmm15		    // b3c3
vpsubw		%zmm24,%zmm8,%zmm24			// b3d3

vmovdqa32		(%r9),%zmm0
vmovdqa32		64(%r9),%zmm1
vpmullw		%zmm0,%zmm18,%zmm2
vpmullw		%zmm0,%zmm20,%zmm3
vpmulhw		%zmm1,%zmm18,%zmm18
vpmulhw		%zmm1,%zmm20,%zmm20
vpmulhw		%zmm16,%zmm2,%zmm2
vpmulhw		%zmm16,%zmm3,%zmm3
vpsubw		%zmm2,%zmm18,%zmm18			// rb0d0
vpsubw		%zmm3,%zmm20,%zmm20			// rb1d1

vmovdqa32		128(%r9),%zmm0
vmovdqa32		192(%r9),%zmm1
vpmullw		%zmm0,%zmm22,%zmm4
vpmullw		%zmm0,%zmm24,%zmm5
vpmulhw		%zmm1,%zmm22,%zmm22
vpmulhw		%zmm1,%zmm24,%zmm24
vpmulhw		%zmm16,%zmm4,%zmm4
vpmulhw		%zmm16,%zmm5,%zmm5
vpsubw		%zmm4,%zmm22,%zmm22			// rb2d2
vpsubw		%zmm5,%zmm24,%zmm24			// rb3d3

vpaddw		%zmm9,%zmm17,%zmm17         //mont(a0d0+b0c0)
vpaddw		%zmm11,%zmm19,%zmm19        //mont(a1d1+b1c1)
vpaddw		%zmm13,%zmm21,%zmm21        //mont(a2d2+b2c2)
vpaddw		%zmm15,%zmm23,%zmm23        //mont(a3d3+b3c3)

vpsubw		%zmm25,%zmm18,%zmm25        //mont(b0d0zeta+a0c0)
vpsubw		%zmm20,%zmm10,%zmm10        //mont(a1c1+b1d1zeta)
vpsubw		%zmm12,%zmm22,%zmm12        //mont(b2d2zeta+a2c2)
vpsubw		%zmm24,%zmm14,%zmm14        //mont(a3c3+b3d3zeta)

vmovdqa32		%zmm25,(256*\off+  0)*2(%rdi)
vmovdqa32		%zmm17,(256*\off+ 32)*2(%rdi)
vmovdqa32		%zmm10,(256*\off+ 64)*2(%rdi)
vmovdqa32		%zmm19,(256*\off+ 96)*2(%rdi)
vmovdqa32		%zmm12,(256*\off+128)*2(%rdi)
vmovdqa32		%zmm21,(256*\off+160)*2(%rdi)
vmovdqa32		%zmm14,(256*\off+192)*2(%rdi)
vmovdqa32		%zmm23,(256*\off+224)*2(%rdi)
.endm

.text
.global cdecl(basemul_avx_32)
cdecl(basemul_avx_32):
mov		%rsp,%r8  //从AVX2到AVX512所有汇编都采用了AT&T格式，所以这里是将rsp中的值mov到r8中
and		$-64,%rsp  //这里用rsp and -64
sub		$64,%rsp  //这里是rsp = rsp - 64，让出了512bit的位置用于schoolbook中存放某个寄存器的值

lea		(_ZETAS_BASEMUL+   0)*2(%rcx),%r9
schoolbook	0

lea		(_ZETAS_BASEMUL+ 128)*2(%rcx),%r9
schoolbook	1

lea		(_ZETAS_BASEMUL+ 256)*2(%rcx),%r9
schoolbook	2

lea		(_ZETAS_BASEMUL+ 384)*2(%rcx),%r9
schoolbook	3

lea		(_ZETAS_BASEMUL+ 512)*2(%rcx),%r9
schoolbook	4

lea		(_ZETAS_BASEMUL+ 640)*2(%rcx),%r9
schoolbook	5

lea		(_ZETAS_BASEMUL+ 768)*2(%rcx),%r9
schoolbook	6

lea		(_ZETAS_BASEMUL+ 896)*2(%rcx),%r9
schoolbook	7

lea		(_ZETAS_BASEMUL+1024)*2(%rcx),%r9
schoolbook	8

lea		(_ZETAS_BASEMUL+1152)*2(%rcx),%r9
schoolbook	9

lea		(_ZETAS_BASEMUL+1280)*2(%rcx),%r9
schoolbook	10

lea		(_ZETAS_BASEMUL+1408)*2(%rcx),%r9
schoolbook	11

lea		(_ZETAS_BASEMUL+1536)*2(%rcx),%r9
schoolbook	12

lea		(_ZETAS_BASEMUL+1664)*2(%rcx),%r9
schoolbook	13

lea		(_ZETAS_BASEMUL+1792)*2(%rcx),%r9
schoolbook	14

lea		(_ZETAS_BASEMUL+1920)*2(%rcx),%r9
schoolbook	15

lea		(_ZETAS_BASEMUL+2048)*2(%rcx),%r9
schoolbook	16

lea		(_ZETAS_BASEMUL+2176)*2(%rcx),%r9
schoolbook	17

lea		(_ZETAS_BASEMUL+2304)*2(%rcx),%r9
schoolbook	18

lea		(_ZETAS_BASEMUL+2432)*2(%rcx),%r9
schoolbook	19

lea		(_ZETAS_BASEMUL+2560)*2(%rcx),%r9
schoolbook	20

lea		(_ZETAS_BASEMUL+2688)*2(%rcx),%r9
schoolbook	21

lea		(_ZETAS_BASEMUL+2816)*2(%rcx),%r9
schoolbook	22

lea		(_ZETAS_BASEMUL+2944)*2(%rcx),%r9
schoolbook	23

lea		(_ZETAS_BASEMUL+3072)*2(%rcx),%r9
schoolbook	24

lea		(_ZETAS_BASEMUL+3200)*2(%rcx),%r9
schoolbook	25

lea		(_ZETAS_BASEMUL+3328)*2(%rcx),%r9
schoolbook	26

lea		(_ZETAS_BASEMUL+3456)*2(%rcx),%r9
schoolbook	27

lea		(_ZETAS_BASEMUL+3584)*2(%rcx),%r9
schoolbook	28

lea		(_ZETAS_BASEMUL+3712)*2(%rcx),%r9
schoolbook	29

lea		(_ZETAS_BASEMUL+3840)*2(%rcx),%r9
schoolbook	30

lea		(_ZETAS_BASEMUL+3968)*2(%rcx),%r9
schoolbook	31


mov		%r8,%rsp
ret























