#include "consts_32.h"
.include "fq.inc"

.macro butterfly4 rl0,rl1,rl2,rl3,rl4,rl5,rl6,rl7,rh0,rh1,rh2,rh3,rh4,rh5,rh6,rh7,zl0=2,zl1=2,zh0=3,zh1=3,off  //butterfly4中的4表示整个butterfly需要2^4 = 16个vmov来装载zeta系数到寄存器中（宏调用前4个vmov + 宏定义中12个vmov）
vpsubw		%zmm\rl0,%zmm\rh0,%zmm21
vpaddw		%zmm\rh0,%zmm\rl0,%zmm\rl0
vpsubw		%zmm\rl1,%zmm\rh1,%zmm22

vpmullw		%zmm\zl0,%zmm21,%zmm\rh0
vpaddw		%zmm\rh1,%zmm\rl1,%zmm\rl1
vpsubw		%zmm\rl2,%zmm\rh2,%zmm23

vpmullw		%zmm\zl1,%zmm22,%zmm\rh1
vpaddw		%zmm\rh2,%zmm\rl2,%zmm\rl2
vpsubw		%zmm\rl3,%zmm\rh3,%zmm24

vmovdqa32     (_32XQINV-(28*(\off/32)+11)*32)*2(%rsi), %zmm\zl1
vpmullw		%zmm\zl1,%zmm23,%zmm\rh2
vpaddw		%zmm\rh3,%zmm\rl3,%zmm\rl3
vpsubw		%zmm\rl4,%zmm\rh4,%zmm25

vmovdqa32     (_32XQINV-(28*(\off/32)+12)*32)*2(%rsi), %zmm\zl1
vpmullw		%zmm\zl1,%zmm24,%zmm\rh3
vpaddw		%zmm\rh4,%zmm\rl4,%zmm\rl4
vpsubw		%zmm\rl5,%zmm\rh5,%zmm26

vmovdqa32     (_32XQINV-(28*(\off/32)+13)*32)*2(%rsi), %zmm\zl1
vpmullw		%zmm\zl1,%zmm25,%zmm\rh4
vpaddw		%zmm\rh5,%zmm\rl5,%zmm\rl5
vpsubw		%zmm\rl6,%zmm\rh6,%zmm27

vmovdqa32     (_32XQINV-(28*(\off/32)+14)*32)*2(%rsi), %zmm\zl1
vpmullw		%zmm\zl1,%zmm26,%zmm\rh5
vpaddw		%zmm\rh6,%zmm\rl6,%zmm\rl6
vpsubw		%zmm\rl7,%zmm\rh7,%zmm28

vmovdqa32     (_32XQINV-(28*(\off/32)+15)*32)*2(%rsi), %zmm\zl1
vpmullw		%zmm\zl1,%zmm27,%zmm\rh6
vpaddw		%zmm\rh7,%zmm\rl7,%zmm\rl7
vmovdqa32     (_32XQINV-(28*(\off/32)+16)*32)*2(%rsi), %zmm\zl1
vpmullw		%zmm\zl1,%zmm28,%zmm\rh7

vpmulhw		%zmm\zh0,%zmm21,%zmm21
vpmulhw		%zmm\zh1,%zmm22,%zmm22

vmovdqa32     (_32XQINV-(28*(\off/32)+3)*32)*2(%rsi), %zmm\zh0
vmovdqa32     (_32XQINV-(28*(\off/32)+4)*32)*2(%rsi), %zmm\zh1
vpmulhw		%zmm\zh0,%zmm23,%zmm23
vpmulhw		%zmm\zh1,%zmm24,%zmm24

vmovdqa32     (_32XQINV-(28*(\off/32)+5)*32)*2(%rsi), %zmm\zh0
vmovdqa32     (_32XQINV-(28*(\off/32)+6)*32)*2(%rsi), %zmm\zh1
vpmulhw		%zmm\zh0,%zmm25,%zmm25
vpmulhw		%zmm\zh1,%zmm26,%zmm26

vmovdqa32     (_32XQINV-(28*(\off/32)+7)*32)*2(%rsi), %zmm\zh0
vmovdqa32     (_32XQINV-(28*(\off/32)+8)*32)*2(%rsi), %zmm\zh1
vpmulhw		%zmm\zh0,%zmm27,%zmm27
vpmulhw		%zmm\zh1,%zmm28,%zmm28

vpmulhw		%zmm0,%zmm\rh0,%zmm\rh0
vpmulhw		%zmm0,%zmm\rh1,%zmm\rh1
vpmulhw		%zmm0,%zmm\rh2,%zmm\rh2
vpmulhw		%zmm0,%zmm\rh3,%zmm\rh3

vpmulhw		%zmm0,%zmm\rh4,%zmm\rh4
vpmulhw		%zmm0,%zmm\rh5,%zmm\rh5
vpmulhw		%zmm0,%zmm\rh6,%zmm\rh6
vpmulhw		%zmm0,%zmm\rh7,%zmm\rh7

#

#

vpsubw		%zmm\rh0,%zmm21,%zmm\rh0
vpsubw		%zmm\rh1,%zmm22,%zmm\rh1
vpsubw		%zmm\rh2,%zmm23,%zmm\rh2
vpsubw		%zmm\rh3,%zmm24,%zmm\rh3

vpsubw		%zmm\rh4,%zmm25,%zmm\rh4
vpsubw		%zmm\rh5,%zmm26,%zmm\rh5
vpsubw		%zmm\rh6,%zmm27,%zmm\rh6
vpsubw		%zmm\rh7,%zmm28,%zmm\rh7

.endm


.macro butterfly4_2 rl0,rl1,rl2,rl3,rl4,rl5,rl6,rl7,rh0,rh1,rh2,rh3,rh4,rh5,rh6,rh7,zl0=2,zl1=2,zh0=3,zh1=3  //butterfly4的第二个版本
vpsubw		%zmm\rl0,%zmm\rh0,%zmm21
vpaddw		%zmm\rh0,%zmm\rl0,%zmm\rl0
vpsubw		%zmm\rl1,%zmm\rh1,%zmm22

vpmullw		%zmm\zl0,%zmm21,%zmm\rh0
vpaddw		%zmm\rh1,%zmm\rl1,%zmm\rl1
vpsubw		%zmm\rl2,%zmm\rh2,%zmm23

vpmullw		%zmm\zl1,%zmm22,%zmm\rh1
vpaddw		%zmm\rh2,%zmm\rl2,%zmm\rl2
vpsubw		%zmm\rl3,%zmm\rh3,%zmm24

vmovdqa32     (992-11*32)*2(%rsi), %zmm\zl1
vpmullw		%zmm\zl1,%zmm23,%zmm\rh2
vpaddw		%zmm\rh3,%zmm\rl3,%zmm\rl3
vpsubw		%zmm\rl4,%zmm\rh4,%zmm25

vmovdqa32     (992-12*32)*2(%rsi), %zmm\zl1
vpmullw		%zmm\zl1,%zmm24,%zmm\rh3
vpaddw		%zmm\rh4,%zmm\rl4,%zmm\rl4
vpsubw		%zmm\rl5,%zmm\rh5,%zmm26

vmovdqa32     (992-13*32)*2(%rsi), %zmm\zl1
vpmullw		%zmm\zl1,%zmm25,%zmm\rh4
vpaddw		%zmm\rh5,%zmm\rl5,%zmm\rl5
vpsubw		%zmm\rl6,%zmm\rh6,%zmm27

vmovdqa32     (992-14*32)*2(%rsi), %zmm\zl1
vpmullw		%zmm\zl1,%zmm26,%zmm\rh5
vpaddw		%zmm\rh6,%zmm\rl6,%zmm\rl6
vpsubw		%zmm\rl7,%zmm\rh7,%zmm28

vmovdqa32     (992-15*32)*2(%rsi), %zmm\zl1
vpmullw		%zmm\zl1,%zmm27,%zmm\rh6
vpaddw		%zmm\rh7,%zmm\rl7,%zmm\rl7
vmovdqa32     (992-16*32)*2(%rsi), %zmm\zl1
vpmullw		%zmm\zl1,%zmm28,%zmm\rh7

vpmulhw		%zmm\zh0,%zmm21,%zmm21
vpmulhw		%zmm\zh1,%zmm22,%zmm22

vmovdqa32     (992-3*32)*2(%rsi), %zmm\zh0
vmovdqa32     (992-4*32)*2(%rsi), %zmm\zh1
vpmulhw		%zmm\zh0,%zmm23,%zmm23
vpmulhw		%zmm\zh1,%zmm24,%zmm24

vmovdqa32     (992-5*32)*2(%rsi), %zmm\zh0
vmovdqa32     (992-6*32)*2(%rsi), %zmm\zh1
vpmulhw		%zmm\zh0,%zmm25,%zmm25
vpmulhw		%zmm\zh1,%zmm26,%zmm26

vmovdqa32     (992-7*32)*2(%rsi), %zmm\zh0
vmovdqa32     (992-8*32)*2(%rsi), %zmm\zh1
vpmulhw		%zmm\zh0,%zmm27,%zmm27
vpmulhw		%zmm\zh1,%zmm28,%zmm28

vpmulhw		%zmm0,%zmm\rh0,%zmm\rh0
vpmulhw		%zmm0,%zmm\rh1,%zmm\rh1
vpmulhw		%zmm0,%zmm\rh2,%zmm\rh2
vpmulhw		%zmm0,%zmm\rh3,%zmm\rh3

vpmulhw		%zmm0,%zmm\rh4,%zmm\rh4
vpmulhw		%zmm0,%zmm\rh5,%zmm\rh5
vpmulhw		%zmm0,%zmm\rh6,%zmm\rh6
vpmulhw		%zmm0,%zmm\rh7,%zmm\rh7

#

#

vpsubw		%zmm\rh0,%zmm21,%zmm\rh0
vpsubw		%zmm\rh1,%zmm22,%zmm\rh1
vpsubw		%zmm\rh2,%zmm23,%zmm\rh2
vpsubw		%zmm\rh3,%zmm24,%zmm\rh3

vpsubw		%zmm\rh4,%zmm25,%zmm\rh4
vpsubw		%zmm\rh5,%zmm26,%zmm\rh5
vpsubw		%zmm\rh6,%zmm27,%zmm\rh6
vpsubw		%zmm\rh7,%zmm28,%zmm\rh7

.endm


.macro butterfly3 rl0,rl1,rl2,rl3,rl4,rl5,rl6,rl7,rh0,rh1,rh2,rh3,rh4,rh5,rh6,rh7,zl0=2,zl1=2,zh0=3,zh1=3,off    //butterfly3中的3表示整个butterfly需要2^3 = 8个vmov来装载zeta系数到寄存器中（宏调用前4个vmov + 宏定义中4个vmov）
vpsubw		%zmm\rl0,%zmm\rh0,%zmm21
vpaddw		%zmm\rh0,%zmm\rl0,%zmm\rl0
vpsubw		%zmm\rl1,%zmm\rh1,%zmm22

vpmullw		%zmm\zl0,%zmm21,%zmm\rh0
vpaddw		%zmm\rh1,%zmm\rl1,%zmm\rl1
vpsubw		%zmm\rl2,%zmm\rh2,%zmm23

vpmullw		%zmm\zl0,%zmm22,%zmm\rh1
vpaddw		%zmm\rh2,%zmm\rl2,%zmm\rl2
vpsubw		%zmm\rl3,%zmm\rh3,%zmm24

vpmullw		%zmm\zl1,%zmm23,%zmm\rh2
vpaddw		%zmm\rh3,%zmm\rl3,%zmm\rl3
vpsubw		%zmm\rl4,%zmm\rh4,%zmm25

vpmullw		%zmm\zl1,%zmm24,%zmm\rh3
vpaddw		%zmm\rh4,%zmm\rl4,%zmm\rl4
vpsubw		%zmm\rl5,%zmm\rh5,%zmm26

vmovdqa32     (_32XQINV-(28*(\off/32)+23)*32)*2(%rsi), %zmm\zl1
vpmullw		%zmm\zl1,%zmm25,%zmm\rh4
vpaddw		%zmm\rh5,%zmm\rl5,%zmm\rl5
vpsubw		%zmm\rl6,%zmm\rh6,%zmm27

vpmullw		%zmm\zl1,%zmm26,%zmm\rh5
vpaddw		%zmm\rh6,%zmm\rl6,%zmm\rl6
vpsubw		%zmm\rl7,%zmm\rh7,%zmm28

vmovdqa32     (_32XQINV-(28*(\off/32)+24)*32)*2(%rsi), %zmm\zl1
vpmullw		%zmm\zl1,%zmm27,%zmm\rh6
vpaddw		%zmm\rh7,%zmm\rl7,%zmm\rl7
vpmullw		%zmm\zl1,%zmm28,%zmm\rh7

vpmulhw		%zmm\zh0,%zmm21,%zmm21
vpmulhw		%zmm\zh0,%zmm22,%zmm22
vpmulhw		%zmm\zh1,%zmm23,%zmm23
vpmulhw		%zmm\zh1,%zmm24,%zmm24

vmovdqa32     (_32XQINV-(28*(\off/32)+19)*32)*2(%rsi), %zmm\zh0
vmovdqa32     (_32XQINV-(28*(\off/32)+20)*32)*2(%rsi), %zmm\zh1
vpmulhw		%zmm\zh0,%zmm25,%zmm25
vpmulhw		%zmm\zh0,%zmm26,%zmm26
vpmulhw		%zmm\zh1,%zmm27,%zmm27
vpmulhw		%zmm\zh1,%zmm28,%zmm28

vpmulhw		%zmm0,%zmm\rh0,%zmm\rh0
vpmulhw		%zmm0,%zmm\rh1,%zmm\rh1
vpmulhw		%zmm0,%zmm\rh2,%zmm\rh2
vpmulhw		%zmm0,%zmm\rh3,%zmm\rh3

vpmulhw		%zmm0,%zmm\rh4,%zmm\rh4
vpmulhw		%zmm0,%zmm\rh5,%zmm\rh5
vpmulhw		%zmm0,%zmm\rh6,%zmm\rh6
vpmulhw		%zmm0,%zmm\rh7,%zmm\rh7

#

#

vpsubw		%zmm\rh0,%zmm21,%zmm\rh0
vpsubw		%zmm\rh1,%zmm22,%zmm\rh1
vpsubw		%zmm\rh2,%zmm23,%zmm\rh2
vpsubw		%zmm\rh3,%zmm24,%zmm\rh3

vpsubw		%zmm\rh4,%zmm25,%zmm\rh4
vpsubw		%zmm\rh5,%zmm26,%zmm\rh5
vpsubw		%zmm\rh6,%zmm27,%zmm\rh6
vpsubw		%zmm\rh7,%zmm28,%zmm\rh7

.endm


.macro butterfly3_2 rl0,rl1,rl2,rl3,rl4,rl5,rl6,rl7,rh0,rh1,rh2,rh3,rh4,rh5,rh6,rh7,zl0=2,zl1=2,zh0=3,zh1=3  //butterfly3的第二个版本
vpsubw		%zmm\rl0,%zmm\rh0,%zmm21
vpaddw		%zmm\rh0,%zmm\rl0,%zmm\rl0
vpsubw		%zmm\rl1,%zmm\rh1,%zmm22

vpmullw		%zmm\zl0,%zmm21,%zmm\rh0
vpaddw		%zmm\rh1,%zmm\rl1,%zmm\rl1
vpsubw		%zmm\rl2,%zmm\rh2,%zmm23

vpmullw		%zmm\zl0,%zmm22,%zmm\rh1
vpaddw		%zmm\rh2,%zmm\rl2,%zmm\rl2
vpsubw		%zmm\rl3,%zmm\rh3,%zmm24

vpmullw		%zmm\zl1,%zmm23,%zmm\rh2
vpaddw		%zmm\rh3,%zmm\rl3,%zmm\rl3
vpsubw		%zmm\rl4,%zmm\rh4,%zmm25

vpmullw		%zmm\zl1,%zmm24,%zmm\rh3
vpaddw		%zmm\rh4,%zmm\rl4,%zmm\rl4
vpsubw		%zmm\rl5,%zmm\rh5,%zmm26

vmovdqa32     (480-7*32)*2(%rsi), %zmm\zl1
vpmullw		%zmm\zl1,%zmm25,%zmm\rh4
vpaddw		%zmm\rh5,%zmm\rl5,%zmm\rl5
vpsubw		%zmm\rl6,%zmm\rh6,%zmm27

vpmullw		%zmm\zl1,%zmm26,%zmm\rh5
vpaddw		%zmm\rh6,%zmm\rl6,%zmm\rl6
vpsubw		%zmm\rl7,%zmm\rh7,%zmm28

vmovdqa32     (480-8*32)*2(%rsi), %zmm\zl1
vpmullw		%zmm\zl1,%zmm27,%zmm\rh6
vpaddw		%zmm\rh7,%zmm\rl7,%zmm\rl7
vpmullw		%zmm\zl1,%zmm28,%zmm\rh7

vpmulhw		%zmm\zh0,%zmm21,%zmm21
vpmulhw		%zmm\zh0,%zmm22,%zmm22
vpmulhw		%zmm\zh1,%zmm23,%zmm23
vpmulhw		%zmm\zh1,%zmm24,%zmm24

vmovdqa32     (480-3*32)*2(%rsi), %zmm\zh0
vmovdqa32     (480-4*32)*2(%rsi), %zmm\zh1
vpmulhw		%zmm\zh0,%zmm25,%zmm25
vpmulhw		%zmm\zh0,%zmm26,%zmm26
vpmulhw		%zmm\zh1,%zmm27,%zmm27
vpmulhw		%zmm\zh1,%zmm28,%zmm28

vpmulhw		%zmm0,%zmm\rh0,%zmm\rh0
vpmulhw		%zmm0,%zmm\rh1,%zmm\rh1
vpmulhw		%zmm0,%zmm\rh2,%zmm\rh2
vpmulhw		%zmm0,%zmm\rh3,%zmm\rh3

vpmulhw		%zmm0,%zmm\rh4,%zmm\rh4
vpmulhw		%zmm0,%zmm\rh5,%zmm\rh5
vpmulhw		%zmm0,%zmm\rh6,%zmm\rh6
vpmulhw		%zmm0,%zmm\rh7,%zmm\rh7

#

#

vpsubw		%zmm\rh0,%zmm21,%zmm\rh0
vpsubw		%zmm\rh1,%zmm22,%zmm\rh1
vpsubw		%zmm\rh2,%zmm23,%zmm\rh2
vpsubw		%zmm\rh3,%zmm24,%zmm\rh3

vpsubw		%zmm\rh4,%zmm25,%zmm\rh4
vpsubw		%zmm\rh5,%zmm26,%zmm\rh5
vpsubw		%zmm\rh6,%zmm27,%zmm\rh6
vpsubw		%zmm\rh7,%zmm28,%zmm\rh7

.endm


.macro butterfly2 rl0,rl1,rl2,rl3,rl4,rl5,rl6,rl7,rh0,rh1,rh2,rh3,rh4,rh5,rh6,rh7,zl0=2,zl1=2,zh0=3,zh1=3    //butterfly2中的2表示整个butterfly需要2^2 = 4个vmov来装载zeta系数到寄存器中（宏调用前的4个vmov）
vpsubw		%zmm\rl0,%zmm\rh0,%zmm21
vpaddw		%zmm\rh0,%zmm\rl0,%zmm\rl0
vpsubw		%zmm\rl1,%zmm\rh1,%zmm22

vpmullw		%zmm\zl0,%zmm21,%zmm\rh0
vpaddw		%zmm\rh1,%zmm\rl1,%zmm\rl1
vpsubw		%zmm\rl2,%zmm\rh2,%zmm23

vpmullw		%zmm\zl0,%zmm22,%zmm\rh1
vpaddw		%zmm\rh2,%zmm\rl2,%zmm\rl2
vpsubw		%zmm\rl3,%zmm\rh3,%zmm24

vpmullw		%zmm\zl0,%zmm23,%zmm\rh2
vpaddw		%zmm\rh3,%zmm\rl3,%zmm\rl3
vpsubw		%zmm\rl4,%zmm\rh4,%zmm25

vpmullw		%zmm\zl0,%zmm24,%zmm\rh3
vpaddw		%zmm\rh4,%zmm\rl4,%zmm\rl4
vpsubw		%zmm\rl5,%zmm\rh5,%zmm26

vpmullw		%zmm\zl1,%zmm25,%zmm\rh4
vpaddw		%zmm\rh5,%zmm\rl5,%zmm\rl5
vpsubw		%zmm\rl6,%zmm\rh6,%zmm27

vpmullw		%zmm\zl1,%zmm26,%zmm\rh5
vpaddw		%zmm\rh6,%zmm\rl6,%zmm\rl6
vpsubw		%zmm\rl7,%zmm\rh7,%zmm28

vpmullw		%zmm\zl1,%zmm27,%zmm\rh6
vpaddw		%zmm\rh7,%zmm\rl7,%zmm\rl7
vpmullw		%zmm\zl1,%zmm28,%zmm\rh7


vpmulhw		%zmm\zh0,%zmm21,%zmm21
vpmulhw		%zmm\zh0,%zmm22,%zmm22
vpmulhw		%zmm\zh0,%zmm23,%zmm23
vpmulhw		%zmm\zh0,%zmm24,%zmm24

vpmulhw		%zmm\zh1,%zmm25,%zmm25
vpmulhw		%zmm\zh1,%zmm26,%zmm26
vpmulhw		%zmm\zh1,%zmm27,%zmm27
vpmulhw		%zmm\zh1,%zmm28,%zmm28

vpmulhw		%zmm0,%zmm\rh0,%zmm\rh0
vpmulhw		%zmm0,%zmm\rh1,%zmm\rh1
vpmulhw		%zmm0,%zmm\rh2,%zmm\rh2
vpmulhw		%zmm0,%zmm\rh3,%zmm\rh3

vpmulhw		%zmm0,%zmm\rh4,%zmm\rh4
vpmulhw		%zmm0,%zmm\rh5,%zmm\rh5
vpmulhw		%zmm0,%zmm\rh6,%zmm\rh6
vpmulhw		%zmm0,%zmm\rh7,%zmm\rh7

#

#

vpsubw		%zmm\rh0,%zmm21,%zmm\rh0
vpsubw		%zmm\rh1,%zmm22,%zmm\rh1
vpsubw		%zmm\rh2,%zmm23,%zmm\rh2
vpsubw		%zmm\rh3,%zmm24,%zmm\rh3

vpsubw		%zmm\rh4,%zmm25,%zmm\rh4
vpsubw		%zmm\rh5,%zmm26,%zmm\rh5
vpsubw		%zmm\rh6,%zmm27,%zmm\rh6
vpsubw		%zmm\rh7,%zmm28,%zmm\rh7

.endm

.macro intt_levels0t2 off
/* level 0 */
vmovdqa32		_32XFLO*2(%rsi),%zmm2
vmovdqa32		_32XFHI*2(%rsi),%zmm3

vmovdqa32         (\off+  0)*32*2(%rdi),%zmm4
vmovdqa32         (\off+  4)*32*2(%rdi),%zmm5
vmovdqa32         (\off+  8)*32*2(%rdi),%zmm6
vmovdqa32         (\off+  12)*32*2(%rdi),%zmm7
vmovdqa32         (\off+  16)*32*2(%rdi),%zmm8
vmovdqa32         (\off+  20)*32*2(%rdi),%zmm9
vmovdqa32         (\off+  24)*32*2(%rdi),%zmm10
vmovdqa32         (\off+  28)*32*2(%rdi),%zmm11

fqmulprecomp	2,3,4
fqmulprecomp	2,3,5
fqmulprecomp	2,3,6
fqmulprecomp	2,3,7
fqmulprecomp	2,3,8
fqmulprecomp	2,3,9
fqmulprecomp	2,3,10
fqmulprecomp	2,3,11

vmovdqa32         (\off+  2)*32*2(%rdi),%zmm12
vmovdqa32         (\off+ 6)*32*2(%rdi),%zmm13
vmovdqa32         (\off+ 10)*32*2(%rdi),%zmm14
vmovdqa32         (\off+ 14)*32*2(%rdi),%zmm15
vmovdqa32         (\off+ 18)*32*2(%rdi),%zmm16
vmovdqa32         (\off+ 22)*32*2(%rdi),%zmm17
vmovdqa32         (\off+ 26)*32*2(%rdi),%zmm18
vmovdqa32         (\off+ 30)*32*2(%rdi),%zmm19

fqmulprecomp	2,3,12
fqmulprecomp	2,3,13
fqmulprecomp	2,3,14
fqmulprecomp	2,3,15
fqmulprecomp	2,3,16
fqmulprecomp	2,3,17
fqmulprecomp	2,3,18
fqmulprecomp	2,3,19

vmovdqa32     (_32XQINV-(28*(\off/32)+9)*32)*2(%rsi), %zmm20
vmovdqa32     (_32XQINV-(28*(\off/32)+10)*32)*2(%rsi), %zmm1
vmovdqa32     (_32XQINV-(28*(\off/32)+1)*32)*2(%rsi), %zmm2
vmovdqa32     (_32XQINV-(28*(\off/32)+2)*32)*2(%rsi), %zmm3

butterfly4	4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1,2,3,\off

/* level 1 */

vmovdqa32     (_32XQINV-(28*(\off/32)+21)*32)*2(%rsi), %zmm20
vmovdqa32     (_32XQINV-(28*(\off/32)+22)*32)*2(%rsi), %zmm1
vmovdqa32     (_32XQINV-(28*(\off/32)+17)*32)*2(%rsi), %zmm2
vmovdqa32     (_32XQINV-(28*(\off/32)+18)*32)*2(%rsi), %zmm3

butterfly3	4,12,6,14,8,16,10,18,5,13,7,15,9,17,11,19,20,1,2,3,\off

/* level 2 */

vmovdqa32     (_32XQINV-(28*(\off/32)+27)*32)*2(%rsi), %zmm20
vmovdqa32     (_32XQINV-(28*(\off/32)+28)*32)*2(%rsi), %zmm1
vmovdqa32     (_32XQINV-(28*(\off/32)+25)*32)*2(%rsi), %zmm2
vmovdqa32     (_32XQINV-(28*(\off/32)+26)*32)*2(%rsi), %zmm3

butterfly2    4,12,5,13,8,16,9,17,6,14,7,15,10,18,11,19,20,1,2,3

vmovdqa32		_32XV*2(%rsi),%zmm1
red32		4
red32       8

vmovdqa32         %zmm4,(\off+  0)*32*2(%rdi)
vmovdqa32         %zmm5,(\off+  4)*32*2(%rdi)
vmovdqa32         %zmm6,(\off+  8)*32*2(%rdi)
vmovdqa32         %zmm7,(\off+ 12)*32*2(%rdi)
vmovdqa32         %zmm8,(\off+ 16)*32*2(%rdi)
vmovdqa32         %zmm9,(\off+ 20)*32*2(%rdi)
vmovdqa32         %zmm10,(\off+ 24)*32*2(%rdi)
vmovdqa32         %zmm11,(\off+ 28)*32*2(%rdi)
vmovdqa32         %zmm12,(\off+  2)*32*2(%rdi)
vmovdqa32         %zmm13,(\off+  6)*32*2(%rdi)
vmovdqa32         %zmm14,(\off+ 10)*32*2(%rdi)
vmovdqa32         %zmm15,(\off+ 14)*32*2(%rdi)
vmovdqa32         %zmm16,(\off+ 18)*32*2(%rdi)
vmovdqa32         %zmm17,(\off+ 22)*32*2(%rdi)
vmovdqa32         %zmm18,(\off+ 26)*32*2(%rdi)
vmovdqa32         %zmm19,(\off+ 30)*32*2(%rdi)

.endm


.macro intt_levels3t6 off
/* level 3 */
vmovdqa32         (\off+  0)*32*2(%rdi),%zmm4
vmovdqa32         (\off+ 16)*32*2(%rdi),%zmm5
vmovdqa32         (\off+ 32)*32*2(%rdi),%zmm6
vmovdqa32         (\off+ 48)*32*2(%rdi),%zmm7
vmovdqa32         (\off+ 64)*32*2(%rdi),%zmm8
vmovdqa32         (\off+ 80)*32*2(%rdi),%zmm9
vmovdqa32         (\off+ 96)*32*2(%rdi),%zmm10
vmovdqa32         (\off+112)*32*2(%rdi),%zmm11

vmovdqa32         (\off+128)*32*2(%rdi),%zmm12
vmovdqa32         (\off+144)*32*2(%rdi),%zmm13
vmovdqa32         (\off+160)*32*2(%rdi),%zmm14
vmovdqa32         (\off+176)*32*2(%rdi),%zmm15
vmovdqa32         (\off+192)*32*2(%rdi),%zmm16
vmovdqa32         (\off+208)*32*2(%rdi),%zmm17
vmovdqa32         (\off+224)*32*2(%rdi),%zmm18
vmovdqa32         (\off+240)*32*2(%rdi),%zmm19

vmovdqa32     (992-9*32)*2(%rsi), %zmm20
vmovdqa32     (992-10*32)*2(%rsi), %zmm1
vmovdqa32     (992-1*32)*2(%rsi), %zmm2
vmovdqa32     (992-2*32)*2(%rsi), %zmm3

butterfly4_2	4,6,8,10,12,14,16,18,5,7,9,11,13,15,17,19,20,1,2,3

vmovdqa32		_32XV*2(%rsi),%zmm1
.if \off == 2
red32		4
red32		6
red32		8
red32		10
red32		12
red32		14
red32		16
red32		18
.endif

.if \off == 3
red32		4
red32		6
red32		8
red32		10
red32		12
red32		14
red32		16
red32		18
.endif

/* level 4 */

vmovdqa32     (480-5*32)*2(%rsi), %zmm20
vmovdqa32     (480-6*32)*2(%rsi), %zmm1
vmovdqa32     (480-1*32)*2(%rsi), %zmm2
vmovdqa32     (480-2*32)*2(%rsi), %zmm3

butterfly3_2	4,5,8,9,12,13,16,17,6,7,10,11,14,15,18,19,20,1,2,3

vmovdqa32		_32XV*2(%rsi),%zmm1
.if \off == 4
red32		4
red32		8
red32		12
red32		16
.endif

.if \off == 5
red32		4
red32		8
red32		12
red32		16
.endif

.if \off == 6
red32		4
red32		8
red32		12
red32		16
.endif

.if \off == 7
red32		4
red32		8
red32		12
red32		16
.endif

/* level 5 */

vmovdqa32     (224-3*32)*2(%rsi), %zmm20
vmovdqa32     (224-4*32)*2(%rsi), %zmm1
vmovdqa32     (224-1*32)*2(%rsi), %zmm2
vmovdqa32     (224-2*32)*2(%rsi), %zmm3

butterfly2	4,5,6,7,12,13,14,15,8,9,10,11,16,17,18,19,20,1,2,3

vmovdqa32		_32XV*2(%rsi),%zmm1
.if \off == 8
red32		4
red32		12
.endif

.if \off == 9
red32		4
red32		12
.endif

.if \off == 10
red32		4
red32		12
.endif

.if \off == 11
red32		4
red32		12
.endif

.if \off == 12
red32		4
red32		12
.endif

.if \off == 13
red32		4
red32		12
.endif

.if \off == 14
red32		4
red32		12
.endif

.if \off == 15
red32		4
red32		12
.endif

.if \off == 0
red32   4
.endif

.if \off == 1
red32   4
.endif

/* level 6 */

vmovdqa32     (_ZETAS_EXP_32 + 0)*2(%rsi), %zmm20
vmovdqa32     (_ZETAS_EXP_32 +32)*2(%rsi), %zmm1

butterfly2	4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,20,1,1

vmovdqa32		_32XV*2(%rsi),%zmm1
/*  //后续研究一下lazy reduction，invntt只有一个数不对可能是范围过大，进行Montgomery约简时超出范围
.if \off == 0
red32   4
.endif

.if \off == 1
red32   4
.endif
*/
red32   5

vmovdqa32         %zmm4,(\off+  0)*32*2(%rdi)
vmovdqa32         %zmm5,(\off+ 16)*32*2(%rdi)
vmovdqa32         %zmm6,(\off+ 32)*32*2(%rdi)
vmovdqa32         %zmm7,(\off+ 48)*32*2(%rdi)
vmovdqa32         %zmm8,(\off+ 64)*32*2(%rdi)
vmovdqa32         %zmm9,(\off+ 80)*32*2(%rdi)
vmovdqa32         %zmm10,(\off+ 96)*32*2(%rdi)
vmovdqa32         %zmm11,(\off+112)*32*2(%rdi)
vmovdqa32         %zmm12,(\off+128)*32*2(%rdi)
vmovdqa32         %zmm13,(\off+144)*32*2(%rdi)
vmovdqa32         %zmm14,(\off+160)*32*2(%rdi)
vmovdqa32         %zmm15,(\off+176)*32*2(%rdi)
vmovdqa32         %zmm16,(\off+192)*32*2(%rdi)
vmovdqa32         %zmm17,(\off+208)*32*2(%rdi)
vmovdqa32         %zmm18,(\off+224)*32*2(%rdi)
vmovdqa32         %zmm19,(\off+240)*32*2(%rdi)
.endm


.text
.global cdecl(invntt_avx_32)
cdecl(invntt_avx_32):
vmovdqa32         _32XQ*2(%rsi),%zmm0

intt_levels0t2	0
intt_levels0t2	1

intt_levels0t2	32
intt_levels0t2	33

intt_levels0t2	64
intt_levels0t2	65

intt_levels0t2	96
intt_levels0t2	97

intt_levels0t2	128
intt_levels0t2	129

intt_levels0t2	160
intt_levels0t2	161

intt_levels0t2	192
intt_levels0t2	193

intt_levels0t2	224
intt_levels0t2	225

/**/
intt_levels3t6  0
intt_levels3t6	1
intt_levels3t6	2
intt_levels3t6	3
intt_levels3t6	4
intt_levels3t6	5
intt_levels3t6	6
intt_levels3t6	7
intt_levels3t6	8
intt_levels3t6	9
intt_levels3t6	10
intt_levels3t6	11
intt_levels3t6	12
intt_levels3t6	13
intt_levels3t6	14
intt_levels3t6	15


ret


