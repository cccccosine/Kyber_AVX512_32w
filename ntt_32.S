#include "consts_32.h"

.macro mul1 rh0,rh1,rh2,rh3,rh4,rh5,rh6,rh7,qnum  # mul1中的1表示整个宏需要2^1 = 2个vmovdqa32来装载qdata中的系数

vmovdqa32   (\qnum    )*2(%rsi),%zmm28
vpmullw		%zmm28,%zmm\rh0,%zmm20
vpmullw		%zmm28,%zmm\rh1,%zmm21
vpmullw		%zmm28,%zmm\rh2,%zmm22
vpmullw		%zmm28,%zmm\rh3,%zmm23
vpmullw		%zmm28,%zmm\rh4,%zmm24
vpmullw		%zmm28,%zmm\rh5,%zmm25
vpmullw		%zmm28,%zmm\rh6,%zmm26
vpmullw		%zmm28,%zmm\rh7,%zmm27

vmovdqa32   (\qnum+ 32)*2(%rsi),%zmm2
vpmulhw		%zmm2,%zmm\rh0,%zmm\rh0
vpmulhw		%zmm2,%zmm\rh1,%zmm\rh1
vpmulhw		%zmm2,%zmm\rh2,%zmm\rh2
vpmulhw		%zmm2,%zmm\rh3,%zmm\rh3
vpmulhw		%zmm2,%zmm\rh4,%zmm\rh4
vpmulhw		%zmm2,%zmm\rh5,%zmm\rh5
vpmulhw		%zmm2,%zmm\rh6,%zmm\rh6
vpmulhw		%zmm2,%zmm\rh7,%zmm\rh7

.endm


.macro mul2 rh0,rh1,rh2,rh3,rh4,rh5,rh6,rh7,qnum  # mul2中的2表示整个宏需要2^2 = 4个vmovdqa32来装载qdata中的系数

vmovdqa32   (\qnum    )*2(%rsi),%zmm28
vmovdqa32   (\qnum+ 32)*2(%rsi),%zmm2
vpmullw		%zmm28,%zmm\rh0,%zmm20
vpmullw		%zmm2,%zmm\rh1,%zmm21
vpmullw		%zmm28,%zmm\rh2,%zmm22
vpmullw		%zmm2,%zmm\rh3,%zmm23
vpmullw     %zmm28,%zmm\rh4,%zmm24
vpmullw     %zmm2,%zmm\rh5,%zmm25
vpmullw     %zmm28,%zmm\rh6,%zmm26
vpmullw     %zmm2,%zmm\rh7,%zmm27

vmovdqa32   (\qnum+ 64)*2(%rsi),%zmm2
vmovdqa32   (\qnum+ 96)*2(%rsi),%zmm1
vpmulhw		%zmm2,%zmm\rh0,%zmm\rh0
vpmulhw		%zmm1,%zmm\rh1,%zmm\rh1
vpmulhw		%zmm2,%zmm\rh2,%zmm\rh2
vpmulhw		%zmm1,%zmm\rh3,%zmm\rh3
vpmulhw		%zmm2,%zmm\rh4,%zmm\rh4
vpmulhw		%zmm1,%zmm\rh5,%zmm\rh5
vpmulhw		%zmm2,%zmm\rh6,%zmm\rh6
vpmulhw		%zmm1,%zmm\rh7,%zmm\rh7

.endm


.macro mul2_2 rh0,rh1,rh2,rh3,rh4,rh5,rh6,rh7,qnum  # mul2_2中的第一个2表示整个宏需要2^2 = 4个vmovdqa32来装载qdata中的系数，第二个2表示是第二种mul2宏

vmovdqa32   (\qnum    )*2(%rsi),%zmm28
vmovdqa32   (\qnum+ 32)*2(%rsi),%zmm2
vpmullw		%zmm28,%zmm\rh0,%zmm20
vpmullw		%zmm28,%zmm\rh1,%zmm21
vpmullw		%zmm28,%zmm\rh2,%zmm22
vpmullw		%zmm28,%zmm\rh3,%zmm23
vpmullw     %zmm2,%zmm\rh4,%zmm24
vpmullw     %zmm2,%zmm\rh5,%zmm25
vpmullw     %zmm2,%zmm\rh6,%zmm26
vpmullw     %zmm2,%zmm\rh7,%zmm27

vmovdqa32   (\qnum+ 64)*2(%rsi),%zmm2
vmovdqa32   (\qnum+ 96)*2(%rsi),%zmm1
vpmulhw		%zmm2,%zmm\rh0,%zmm\rh0
vpmulhw		%zmm2,%zmm\rh1,%zmm\rh1
vpmulhw		%zmm2,%zmm\rh2,%zmm\rh2
vpmulhw		%zmm2,%zmm\rh3,%zmm\rh3
vpmulhw		%zmm1,%zmm\rh4,%zmm\rh4
vpmulhw		%zmm1,%zmm\rh5,%zmm\rh5
vpmulhw		%zmm1,%zmm\rh6,%zmm\rh6
vpmulhw		%zmm1,%zmm\rh7,%zmm\rh7

.endm


.macro mul3 rh0,rh1,rh2,rh3,rh4,rh5,rh6,rh7,qnum  # mul3中的3表示整个宏需要2^3 = 8个vmovdqa32来装载qdata中的系数

vmovdqa32   (\qnum    )*2(%rsi),%zmm28
vpmullw		%zmm28,%zmm\rh0,%zmm20
vpmullw		%zmm28,%zmm\rh4,%zmm24

vmovdqa32   (\qnum+ 32)*2(%rsi),%zmm28
vpmullw		%zmm28,%zmm\rh1,%zmm21
vpmullw		%zmm28,%zmm\rh5,%zmm25

vmovdqa32   (\qnum+ 64)*2(%rsi),%zmm28
vpmullw		%zmm28,%zmm\rh2,%zmm22
vpmullw		%zmm28,%zmm\rh6,%zmm26

vmovdqa32   (\qnum+ 96)*2(%rsi),%zmm28
vpmullw		%zmm28,%zmm\rh3,%zmm23
vpmullw		%zmm28,%zmm\rh7,%zmm27

////////////////////////////////////////////

vmovdqa32   (\qnum+ 128)*2(%rsi),%zmm2
vpmulhw		%zmm2,%zmm\rh0,%zmm\rh0
vpmulhw		%zmm2,%zmm\rh4,%zmm\rh4

vmovdqa32   (\qnum+ 160)*2(%rsi),%zmm2
vpmulhw		%zmm2,%zmm\rh1,%zmm\rh1
vpmulhw		%zmm2,%zmm\rh5,%zmm\rh5

vmovdqa32   (\qnum+ 192)*2(%rsi),%zmm2
vpmulhw		%zmm2,%zmm\rh2,%zmm\rh2
vpmulhw		%zmm2,%zmm\rh6,%zmm\rh6

vmovdqa32   (\qnum+224)*2(%rsi),%zmm2
vpmulhw		%zmm2,%zmm\rh3,%zmm\rh3
vpmulhw		%zmm2,%zmm\rh7,%zmm\rh7

.endm


.macro mul3_2 rh0,rh1,rh2,rh3,rh4,rh5,rh6,rh7,qnum  # mul3中的3表示整个宏需要2^3 = 8个vmovdqa32来装载qdata中的系数，第二个2表示是第二种mul3宏

vmovdqa32   (\qnum    )*2(%rsi),%zmm28
vpmullw		%zmm28,%zmm\rh0,%zmm20
vpmullw		%zmm28,%zmm\rh2,%zmm22

vmovdqa32   (\qnum+ 32)*2(%rsi),%zmm28
vpmullw		%zmm28,%zmm\rh1,%zmm21
vpmullw		%zmm28,%zmm\rh3,%zmm23

vmovdqa32   (\qnum+ 64)*2(%rsi),%zmm28
vpmullw		%zmm28,%zmm\rh4,%zmm24
vpmullw		%zmm28,%zmm\rh6,%zmm26

vmovdqa32   (\qnum+ 96)*2(%rsi),%zmm28
vpmullw		%zmm28,%zmm\rh5,%zmm25
vpmullw		%zmm28,%zmm\rh7,%zmm27

////////////////////////////////////////////

vmovdqa32   (\qnum+ 128)*2(%rsi),%zmm2
vpmulhw		%zmm2,%zmm\rh0,%zmm\rh0
vpmulhw		%zmm2,%zmm\rh2,%zmm\rh2

vmovdqa32   (\qnum+ 160)*2(%rsi),%zmm2
vpmulhw		%zmm2,%zmm\rh1,%zmm\rh1
vpmulhw		%zmm2,%zmm\rh3,%zmm\rh3

vmovdqa32   (\qnum+ 192)*2(%rsi),%zmm2
vpmulhw		%zmm2,%zmm\rh4,%zmm\rh4
vpmulhw		%zmm2,%zmm\rh6,%zmm\rh6

vmovdqa32   (\qnum+224)*2(%rsi),%zmm2
vpmulhw		%zmm2,%zmm\rh5,%zmm\rh5
vpmulhw		%zmm2,%zmm\rh7,%zmm\rh7

.endm


.macro mul4 rh0,rh1,rh2,rh3,rh4,rh5,rh6,rh7,qnum  # mul4中的4表示整个宏需要2^4 = 16个vmovdqa32来装载qdata中的系数

vmovdqa32		(\qnum    )*2(%rsi),%zmm28
vpmullw		%zmm28,%zmm\rh0,%zmm20
vmovdqa32		(\qnum+ 32)*2(%rsi),%zmm28
vpmullw		%zmm28,%zmm\rh1,%zmm21
vmovdqa32		(\qnum+ 64)*2(%rsi),%zmm28
vpmullw		%zmm28,%zmm\rh2,%zmm22
vmovdqa32		(\qnum+ 96)*2(%rsi),%zmm28
vpmullw		%zmm28,%zmm\rh3,%zmm23
vmovdqa32		(\qnum+128)*2(%rsi),%zmm28
vpmullw		%zmm28,%zmm\rh4,%zmm24
vmovdqa32		(\qnum+160)*2(%rsi),%zmm28
vpmullw		%zmm28,%zmm\rh5,%zmm25
vmovdqa32		(\qnum+192)*2(%rsi),%zmm28
vpmullw		%zmm28,%zmm\rh6,%zmm26
vmovdqa32		(\qnum+224)*2(%rsi),%zmm28
vpmullw		%zmm28,%zmm\rh7,%zmm27

////////////////////////////////////////////

vmovdqa32		(\qnum+256)*2(%rsi),%zmm2
vpmulhw		%zmm2,%zmm\rh0,%zmm\rh0
vmovdqa32		(\qnum+288)*2(%rsi),%zmm2
vpmulhw		%zmm2,%zmm\rh1,%zmm\rh1
vmovdqa32		(\qnum+320)*2(%rsi),%zmm2
vpmulhw		%zmm2,%zmm\rh2,%zmm\rh2
vmovdqa32		(\qnum+352)*2(%rsi),%zmm2
vpmulhw		%zmm2,%zmm\rh3,%zmm\rh3
vmovdqa32		(\qnum+384)*2(%rsi),%zmm2
vpmulhw		%zmm2,%zmm\rh4,%zmm\rh4
vmovdqa32		(\qnum+416)*2(%rsi),%zmm2
vpmulhw		%zmm2,%zmm\rh5,%zmm\rh5
vmovdqa32		(\qnum+448)*2(%rsi),%zmm2
vpmulhw		%zmm2,%zmm\rh6,%zmm\rh6
vmovdqa32		(\qnum+480)*2(%rsi),%zmm2
vpmulhw		%zmm2,%zmm\rh7,%zmm\rh7

.endm


.macro reduce
vpmulhw		%zmm0,%zmm20,%zmm20
vpmulhw		%zmm0,%zmm21,%zmm21

vpmulhw		%zmm0,%zmm22,%zmm22
vpmulhw		%zmm0,%zmm23,%zmm23

vpmulhw		%zmm0,%zmm24,%zmm24
vpmulhw		%zmm0,%zmm25,%zmm25

vpmulhw		%zmm0,%zmm26,%zmm26
vpmulhw		%zmm0,%zmm27,%zmm27

.endm

.macro update rln,rl0,rl1,rl2,rl3,rl4,rl5,rl6,rl7,rh0,rh1,rh2,rh3,rh4,rh5,rh6,rh7
vpaddw		%zmm\rh0,%zmm\rl0,%zmm\rln
vpsubw		%zmm\rh0,%zmm\rl0,%zmm\rh0
vpaddw		%zmm\rh1,%zmm\rl1,%zmm\rl0

vpsubw		%zmm\rh1,%zmm\rl1,%zmm\rh1
vpaddw		%zmm\rh2,%zmm\rl2,%zmm\rl1
vpsubw		%zmm\rh2,%zmm\rl2,%zmm\rh2

vpaddw		%zmm\rh3,%zmm\rl3,%zmm\rl2
vpsubw		%zmm\rh3,%zmm\rl3,%zmm\rh3

vpaddw		%zmm\rh4,%zmm\rl4,%zmm\rl3
vpsubw		%zmm\rh4,%zmm\rl4,%zmm\rh4
vpaddw		%zmm\rh5,%zmm\rl5,%zmm\rl4

vpsubw		%zmm\rh5,%zmm\rl5,%zmm\rh5
vpaddw		%zmm\rh6,%zmm\rl6,%zmm\rl5
vpsubw		%zmm\rh6,%zmm\rl6,%zmm\rh6

vpaddw		%zmm\rh7,%zmm\rl7,%zmm\rl6
vpsubw		%zmm\rh7,%zmm\rl7,%zmm\rh7

////////////////////////////////////////////

vpsubw		%zmm20,%zmm\rln,%zmm\rln
vpaddw		%zmm20,%zmm\rh0,%zmm\rh0
vpsubw		%zmm21,%zmm\rl0,%zmm\rl0

vpaddw		%zmm21,%zmm\rh1,%zmm\rh1
vpsubw		%zmm22,%zmm\rl1,%zmm\rl1
vpaddw		%zmm22,%zmm\rh2,%zmm\rh2

vpsubw		%zmm23,%zmm\rl2,%zmm\rl2
vpaddw		%zmm23,%zmm\rh3,%zmm\rh3

vpsubw		%zmm24,%zmm\rl3,%zmm\rl3
vpaddw		%zmm24,%zmm\rh4,%zmm\rh4
vpsubw		%zmm25,%zmm\rl4,%zmm\rl4

vpaddw		%zmm25,%zmm\rh5,%zmm\rh5
vpsubw		%zmm26,%zmm\rl5,%zmm\rl5
vpaddw		%zmm26,%zmm\rh6,%zmm\rh6

vpsubw		%zmm27,%zmm\rl6,%zmm\rl6
vpaddw		%zmm27,%zmm\rh7,%zmm\rh7
.endm


.macro levels0t3 off
/* level0 */
vmovdqa32		(\off+128)*32*2(%rdi),%zmm12  #每行有*2是因为地址加偏移量单位是字节，一个Kyber系数是16b=2字节
vmovdqa32		(\off+144)*32*2(%rdi),%zmm13
vmovdqa32		(\off+160)*32*2(%rdi),%zmm14
vmovdqa32		(\off+176)*32*2(%rdi),%zmm15
vmovdqa32		(\off+192)*32*2(%rdi),%zmm16
vmovdqa32		(\off+208)*32*2(%rdi),%zmm17
vmovdqa32		(\off+224)*32*2(%rdi),%zmm18
vmovdqa32		(\off+240)*32*2(%rdi),%zmm19

mul1        12,13,14,15,16,17,18,19,32

vmovdqa32		(\off+  0)*32*2(%rdi),%zmm4
vmovdqa32		(\off+ 16)*32*2(%rdi),%zmm5
vmovdqa32		(\off+ 32)*32*2(%rdi),%zmm6
vmovdqa32		(\off+ 48)*32*2(%rdi),%zmm7
vmovdqa32		(\off+ 64)*32*2(%rdi),%zmm8
vmovdqa32		(\off+ 80)*32*2(%rdi),%zmm9
vmovdqa32		(\off+ 96)*32*2(%rdi),%zmm10
vmovdqa32		(\off+112)*32*2(%rdi),%zmm11

reduce
update		3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19

/* level1 */
mul2        7,16,8,17,9,18,10,19,96

reduce
update		11,3,12,4,13,5,14,6,15,7,16,8,17,9,18,10,19

/* level2 */
mul3        13,9,5,18,14,10,6,19,224

reduce
update		15,11,7,3,16,12,8,4,17,13,9,5,18,14,10,6,19

/* level3 */
mul4        16,14,12,10,8,6,4,19,480

reduce
update      17,15,13,11,9,7,5,3,18,16,14,12,10,8,6,4,19

vmovdqa32		%zmm17,(\off+  0)*32*2(%rdi)
vmovdqa32		%zmm16,(\off+ 16)*32*2(%rdi)
vmovdqa32		%zmm15,(\off+ 32)*32*2(%rdi)
vmovdqa32		%zmm14,(\off+ 48)*32*2(%rdi)
vmovdqa32		%zmm13,(\off+ 64)*32*2(%rdi)
vmovdqa32		%zmm12,(\off+ 80)*32*2(%rdi)
vmovdqa32		%zmm11,(\off+ 96)*32*2(%rdi)
vmovdqa32		%zmm10,(\off+112)*32*2(%rdi)
vmovdqa32		%zmm9,(\off+128)*32*2(%rdi)
vmovdqa32		%zmm8,(\off+144)*32*2(%rdi)
vmovdqa32		%zmm7,(\off+160)*32*2(%rdi)
vmovdqa32		%zmm6,(\off+176)*32*2(%rdi)
vmovdqa32		%zmm5,(\off+192)*32*2(%rdi)
vmovdqa32		%zmm4,(\off+208)*32*2(%rdi)
vmovdqa32		%zmm3,(\off+224)*32*2(%rdi)
vmovdqa32		%zmm19,(\off+240)*32*2(%rdi)

.endm


.macro levels4t6 off
/* level4 */
vmovdqa32		(\off+  8)*32*2(%rdi),%zmm12
vmovdqa32		(\off+ 10)*32*2(%rdi),%zmm13
vmovdqa32		(\off+ 12)*32*2(%rdi),%zmm14
vmovdqa32		(\off+ 14)*32*2(%rdi),%zmm15
vmovdqa32		(\off+ 24)*32*2(%rdi),%zmm16
vmovdqa32		(\off+ 26)*32*2(%rdi),%zmm17
vmovdqa32		(\off+ 28)*32*2(%rdi),%zmm18
vmovdqa32		(\off+ 30)*32*2(%rdi),%zmm19

mul2_2		12,13,14,15,16,17,18,19,(992+(\off/32)*896)

vmovdqa32		(\off+  0)*32*2(%rdi),%zmm4
vmovdqa32		(\off+  2)*32*2(%rdi),%zmm5
vmovdqa32		(\off+  4)*32*2(%rdi),%zmm6
vmovdqa32		(\off+  6)*32*2(%rdi),%zmm7
vmovdqa32		(\off+ 16)*32*2(%rdi),%zmm8
vmovdqa32		(\off+ 18)*32*2(%rdi),%zmm9
vmovdqa32		(\off+ 20)*32*2(%rdi),%zmm10
vmovdqa32		(\off+ 22)*32*2(%rdi),%zmm11

reduce
update		3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19

/* level5 */
mul3_2      5,14,6,15,9,18,10,19,(1120+(\off/32)*896)

reduce
update		11,3,12,4,13,7,16,8,17,5,14,6,15,9,18,10,19

/* level6 */
mul4        12,6,4,15,16,10,8,19,(1376+(\off/32)*896)

reduce
update		17,11,5,3,14,13,9,7,18,12,6,4,15,16,10,8,19

vmovdqa32		%zmm17,(\off+  0)*32*2(%rdi)
vmovdqa32		%zmm12,(\off+  2)*32*2(%rdi)
vmovdqa32		%zmm11,(\off+  4)*32*2(%rdi)
vmovdqa32		%zmm6,(\off+ 6)*32*2(%rdi)
vmovdqa32		%zmm5,(\off+ 8)*32*2(%rdi)
vmovdqa32		%zmm4,(\off+ 10)*32*2(%rdi)
vmovdqa32		%zmm3,(\off+ 12)*32*2(%rdi)
vmovdqa32		%zmm15,(\off+ 14)*32*2(%rdi)
vmovdqa32		%zmm14,(\off+ 16)*32*2(%rdi)
vmovdqa32		%zmm16,(\off+ 18)*32*2(%rdi)
vmovdqa32		%zmm13,(\off+ 20)*32*2(%rdi)
vmovdqa32		%zmm10,(\off+ 22)*32*2(%rdi)
vmovdqa32		%zmm9,(\off+ 24)*32*2(%rdi)
vmovdqa32		%zmm8,(\off+ 26)*32*2(%rdi)
vmovdqa32		%zmm7,(\off+ 28)*32*2(%rdi)
vmovdqa32		%zmm19,(\off+ 30)*32*2(%rdi)

.endm



.text
.global cdecl(ntt_avx_32)
cdecl(ntt_avx_32):
vmovdqa32		_32XQ*2(%rsi),%zmm0

levels0t3  0
levels0t3  1
levels0t3  2
levels0t3  3
levels0t3  4
levels0t3  5
levels0t3  6
levels0t3  7
levels0t3  8
levels0t3  9
levels0t3  10
levels0t3  11
levels0t3  12
levels0t3  13
levels0t3  14
levels0t3  15


levels4t6  0
levels4t6  1

levels4t6  32
levels4t6  33

levels4t6  64
levels4t6  65

levels4t6  96
levels4t6  97

levels4t6  128
levels4t6  129

levels4t6  160
levels4t6  161

levels4t6  192
levels4t6  193

levels4t6  224
levels4t6  225


ret
