#include "consts_32.h"
.include "shuffle.inc"

//下一行的256b表示处理一个多项式中16*16b = 256b，有一点点歧义，这个宏只是处理了32-way的各16个数，但是最后并没有变成32-way格式
.macro formseqto32_AVX512_256b p q r # param1表示输入的位置，param2表示输出的位置，param3表示该多项式有多少系数

vmovdqu32		(\p      )*2(%rdi),%zmm4  #使用vmovdqa需要先将数据对齐，使用vmovdqu32就不用考虑对齐
vmovdqu32		(\p+   \r)*2(%rdi),%zmm5
vmovdqu32		(\p+ \r*2)*2(%rdi),%zmm6
vmovdqu32		(\p+ \r*3)*2(%rdi),%zmm7
vmovdqu32		(\p+ \r*4)*2(%rdi),%zmm8
vmovdqu32		(\p+ \r*5)*2(%rdi),%zmm9
vmovdqu32		(\p+ \r*6)*2(%rdi),%zmm10
vmovdqu32		(\p+ \r*7)*2(%rdi),%zmm11
vmovdqu32		(\p+ \r*8)*2(%rdi),%zmm12
vmovdqu32		(\p+ \r*9)*2(%rdi),%zmm13
vmovdqu32		(\p+\r*10)*2(%rdi),%zmm14
vmovdqu32		(\p+\r*11)*2(%rdi),%zmm15
vmovdqu32		(\p+\r*12)*2(%rdi),%zmm16
vmovdqu32		(\p+\r*13)*2(%rdi),%zmm17
vmovdqu32		(\p+\r*14)*2(%rdi),%zmm18
vmovdqu32		(\p+\r*15)*2(%rdi),%zmm19

shuffle16   4,12,3,12
shuffle16   5,13,4,13
shuffle16   6,14,5,14
shuffle16   7,15,6,15
shuffle16   8,16,7,16
shuffle16   9,17,8,17
shuffle16   10,18,9,18
shuffle16   11,19,10,19

shuffle8    3,7,11,7
shuffle8	4,8,3,8
shuffle8	5,9,4,9
shuffle8	6,10,5,10
shuffle8    12,16,6,16
shuffle8    13,17,12,17
shuffle8    14,18,13,18
shuffle8    15,19,14,19

shuffle4	11,4,15,4
shuffle4	3,5,11,5
shuffle4	7,9,3,9
shuffle4	8,10,7,10
shuffle4	6,13,8,13
shuffle4	12,14,6,14
shuffle4	16,18,12,18
shuffle4	17,19,16,19

shuffle2	15,11,17,11
shuffle2	4,5,15,5
shuffle2	3,7,4,7
shuffle2	9,10,3,10
shuffle2	8,6,9,6
shuffle2	13,14,8,14
shuffle2	12,16,13,16
shuffle2	18,19,12,19
/*
shuffle1	10,4,14,4
shuffle1	15,6,10,6
shuffle1	11,12,15,12
shuffle1	3,13,11,13
shuffle1	5,16,3,16
shuffle1	7,17,5,17
shuffle1	8,18,7,18
shuffle1	9,19,8,19
*/

vmovdqa32       SEQ_SHUFIDX_32*2(%rcx),%zmm0

vpshufb     %zmm0, %zmm17, %zmm17
vpshufb     %zmm0, %zmm11, %zmm11
vpshufb     %zmm0, %zmm15, %zmm15
vpshufb     %zmm0, %zmm5, %zmm5
vpshufb     %zmm0, %zmm4, %zmm4
vpshufb     %zmm0, %zmm7, %zmm7
vpshufb     %zmm0, %zmm3, %zmm3
vpshufb     %zmm0, %zmm10, %zmm10
vpshufb     %zmm0, %zmm9, %zmm9
vpshufb     %zmm0, %zmm6, %zmm6
vpshufb     %zmm0, %zmm8, %zmm8
vpshufb     %zmm0, %zmm14, %zmm14
vpshufb     %zmm0, %zmm13, %zmm13
vpshufb     %zmm0, %zmm16, %zmm16
vpshufb     %zmm0, %zmm12, %zmm12
vpshufb     %zmm0, %zmm19, %zmm19

vmovdqa32   SEQ_PERMIDX_32*2(%rcx),%zmm1

vpermq      %zmm17,%zmm1,%zmm17
vpermq      %zmm11,%zmm1,%zmm11
vpermq      %zmm15,%zmm1,%zmm15
vpermq      %zmm5,%zmm1,%zmm5
vpermq      %zmm4,%zmm1,%zmm4
vpermq      %zmm7,%zmm1,%zmm7
vpermq      %zmm3,%zmm1,%zmm3
vpermq      %zmm10,%zmm1,%zmm10
vpermq      %zmm9,%zmm1,%zmm9
vpermq      %zmm6,%zmm1,%zmm6
vpermq      %zmm8,%zmm1,%zmm8
vpermq      %zmm14,%zmm1,%zmm14
vpermq      %zmm13,%zmm1,%zmm13
vpermq      %zmm16,%zmm1,%zmm16
vpermq      %zmm12,%zmm1,%zmm12
vpermq      %zmm19,%zmm1,%zmm19
/*
shuffle16    10,5,6,5  //shuffle16未改，但是也不需要shuffle16了
shuffle16    9,4,10,4
shuffle16    8,3,9,3
shuffle16    7,11,8,11
*/
vmovdqu32		%zmm17,(\q     )*2(%rsi)
vmovdqu32		%zmm11,(\q+  64)*2(%rsi)
vmovdqu32		%zmm15,(\q+ 128)*2(%rsi)
vmovdqu32		%zmm5, (\q+ 192)*2(%rsi)
vmovdqu32		%zmm4, (\q+ 256)*2(%rsi)
vmovdqu32		%zmm7, (\q+ 320)*2(%rsi)
vmovdqu32		%zmm3, (\q+ 384)*2(%rsi)
vmovdqu32		%zmm10,(\q+ 448)*2(%rsi)
vmovdqu32		%zmm9, (\q+ 512)*2(%rsi)
vmovdqu32		%zmm6, (\q+ 576)*2(%rsi)
vmovdqu32		%zmm8, (\q+ 640)*2(%rsi)
vmovdqu32		%zmm14,(\q+ 704)*2(%rsi)
vmovdqu32		%zmm13,(\q+ 768)*2(%rsi)
vmovdqu32		%zmm16,(\q+ 832)*2(%rsi)
vmovdqu32		%zmm12,(\q+ 896)*2(%rsi)
vmovdqu32		%zmm19,(\q+ 960)*2(%rsi)

.endm

.macro formseq_supp p q #补充的位置调换函数，每次调换32*30个系数的位置，对于to16/from16通用，因为使用两次formseq_supp就会相互抵消

vmovdqu32		(\p     )*2(%rsi),%zmm1
vmovdqu32		(\p+  32)*2(%rsi),%zmm2
vmovdqu32		(\p+  64)*2(%rsi),%zmm3
vmovdqu32		(\p+  96)*2(%rsi),%zmm4
vmovdqu32		(\p+ 128)*2(%rsi),%zmm5
vmovdqu32		(\p+ 160)*2(%rsi),%zmm6
vmovdqu32		(\p+ 192)*2(%rsi),%zmm7
vmovdqu32		(\p+ 224)*2(%rsi),%zmm8
vmovdqu32		(\p+ 256)*2(%rsi),%zmm9
vmovdqu32		(\p+ 288)*2(%rsi),%zmm10
vmovdqu32		(\p+ 320)*2(%rsi),%zmm11
vmovdqu32		(\p+ 352)*2(%rsi),%zmm12
vmovdqu32		(\p+ 384)*2(%rsi),%zmm13
vmovdqu32		(\p+ 416)*2(%rsi),%zmm14
vmovdqu32		(\p+ 448)*2(%rsi),%zmm15

vmovdqu32		(\p+ 480)*2(%rsi),%zmm16
vmovdqu32		(\p+ 512)*2(%rsi),%zmm17
vmovdqu32		(\p+ 544)*2(%rsi),%zmm18
vmovdqu32		(\p+ 576)*2(%rsi),%zmm19
vmovdqu32		(\p+ 608)*2(%rsi),%zmm20
vmovdqu32		(\p+ 640)*2(%rsi),%zmm21
vmovdqu32		(\p+ 672)*2(%rsi),%zmm22
vmovdqu32		(\p+ 704)*2(%rsi),%zmm23
vmovdqu32		(\p+ 736)*2(%rsi),%zmm24
vmovdqu32		(\p+ 768)*2(%rsi),%zmm25
vmovdqu32		(\p+ 800)*2(%rsi),%zmm26
vmovdqu32		(\p+ 832)*2(%rsi),%zmm27
vmovdqu32		(\p+ 864)*2(%rsi),%zmm28
vmovdqu32		(\p+ 896)*2(%rsi),%zmm29
vmovdqu32		(\p+ 928)*2(%rsi),%zmm30

shuffle16   1,2,0,1
shuffle16   3,4,2,3
shuffle16   5,6,4,5
shuffle16   7,8,6,7
shuffle16   9,10,8,9
shuffle16   11,12,10,11
shuffle16   13,14,12,13
shuffle16   15,16,14,15
shuffle16   17,18,16,17
shuffle16   19,20,18,19
shuffle16   21,22,20,21
shuffle16   23,24,22,23
shuffle16   25,26,24,25
shuffle16   27,28,26,27
shuffle16   29,30,28,29

vmovdqu32		%zmm0,(\q     )*2(%rdx)
vmovdqu32		%zmm1,(\q+  32)*2(%rdx)
vmovdqu32		%zmm2,(\q+  64)*2(%rdx)
vmovdqu32		%zmm3,(\q+  96)*2(%rdx)
vmovdqu32		%zmm4,(\q+ 128)*2(%rdx)
vmovdqu32		%zmm5,(\q+ 160)*2(%rdx)
vmovdqu32		%zmm6,(\q+ 192)*2(%rdx)
vmovdqu32		%zmm7,(\q+ 224)*2(%rdx)
vmovdqu32		%zmm8,(\q+ 256)*2(%rdx)
vmovdqu32		%zmm9,(\q+ 288)*2(%rdx)
vmovdqu32		%zmm10,(\q+ 320)*2(%rdx)
vmovdqu32		%zmm11,(\q+ 352)*2(%rdx)
vmovdqu32		%zmm12,(\q+ 384)*2(%rdx)
vmovdqu32		%zmm13,(\q+ 416)*2(%rdx)
vmovdqu32		%zmm14,(\q+ 448)*2(%rdx)

vmovdqu32		%zmm15,(\q+ 480)*2(%rdx)
vmovdqu32		%zmm16,(\q+ 512)*2(%rdx)
vmovdqu32		%zmm17,(\q+ 544)*2(%rdx)
vmovdqu32		%zmm18,(\q+ 576)*2(%rdx)
vmovdqu32		%zmm19,(\q+ 608)*2(%rdx)
vmovdqu32		%zmm20,(\q+ 640)*2(%rdx)
vmovdqu32		%zmm21,(\q+ 672)*2(%rdx)
vmovdqu32		%zmm22,(\q+ 704)*2(%rdx)
vmovdqu32		%zmm23,(\q+ 736)*2(%rdx)
vmovdqu32		%zmm24,(\q+ 768)*2(%rdx)
vmovdqu32		%zmm25,(\q+ 800)*2(%rdx)
vmovdqu32		%zmm26,(\q+ 832)*2(%rdx)
vmovdqu32		%zmm27,(\q+ 864)*2(%rdx)
vmovdqu32		%zmm28,(\q+ 896)*2(%rdx)
vmovdqu32		%zmm29,(\q+ 928)*2(%rdx)

.endm

.macro formseqfrom32_AVX512_256b p q r

vmovdqu32		(\p     )*2(%rdx),%zmm17
vmovdqu32		(\p+  64)*2(%rdx),%zmm11
vmovdqu32		(\p+ 128)*2(%rdx),%zmm15
vmovdqu32		(\p+ 192)*2(%rdx),%zmm5
vmovdqu32		(\p+ 256)*2(%rdx),%zmm4
vmovdqu32		(\p+ 320)*2(%rdx),%zmm7
vmovdqu32		(\p+ 384)*2(%rdx),%zmm3
vmovdqu32		(\p+ 448)*2(%rdx),%zmm10
vmovdqu32		(\p+ 512)*2(%rdx),%zmm9
vmovdqu32		(\p+ 576)*2(%rdx),%zmm6
vmovdqu32		(\p+ 640)*2(%rdx),%zmm8
vmovdqu32		(\p+ 704)*2(%rdx),%zmm14
vmovdqu32		(\p+ 768)*2(%rdx),%zmm13
vmovdqu32		(\p+ 832)*2(%rdx),%zmm16
vmovdqu32		(\p+ 896)*2(%rdx),%zmm12
vmovdqu32		(\p+ 960)*2(%rdx),%zmm19

vmovdqa32   SEQ_INVSHUFIDX_32*2(%rcx),%zmm0
vmovdqa32   SEQ_INVPERMIDX_32*2(%rcx),%zmm1
/*
shuffle8    8,11,7,11
shuffle8    9,3,8,3
shuffle8    10,4,9,4
shuffle8    6,5,10,5
*/
vpermq      %zmm17,%zmm1,%zmm17
vpermq      %zmm11,%zmm1,%zmm11
vpermq      %zmm15,%zmm1,%zmm15
vpermq      %zmm5,%zmm1,%zmm5
vpermq      %zmm4,%zmm1,%zmm4
vpermq      %zmm7,%zmm1,%zmm7
vpermq      %zmm3,%zmm1,%zmm3
vpermq      %zmm10,%zmm1,%zmm10
vpermq      %zmm9,%zmm1,%zmm9
vpermq      %zmm6,%zmm1,%zmm6
vpermq      %zmm8,%zmm1,%zmm8
vpermq      %zmm14,%zmm1,%zmm14
vpermq      %zmm13,%zmm1,%zmm13
vpermq      %zmm16,%zmm1,%zmm16
vpermq      %zmm12,%zmm1,%zmm12
vpermq      %zmm19,%zmm1,%zmm19

vpshufb     %zmm0, %zmm17, %zmm17
vpshufb     %zmm0, %zmm11, %zmm11
vpshufb     %zmm0, %zmm15, %zmm15
vpshufb     %zmm0, %zmm5, %zmm5
vpshufb     %zmm0, %zmm4, %zmm4
vpshufb     %zmm0, %zmm7, %zmm7
vpshufb     %zmm0, %zmm3, %zmm3
vpshufb     %zmm0, %zmm10, %zmm10
vpshufb     %zmm0, %zmm9, %zmm9
vpshufb     %zmm0, %zmm6, %zmm6
vpshufb     %zmm0, %zmm8, %zmm8
vpshufb     %zmm0, %zmm14, %zmm14
vpshufb     %zmm0, %zmm13, %zmm13
vpshufb     %zmm0, %zmm16, %zmm16
vpshufb     %zmm0, %zmm12, %zmm12
vpshufb     %zmm0, %zmm19, %zmm19

shuffle2	12,19,18,19
shuffle2	13,16,12,16
shuffle2	8,14,13,14
shuffle2	9,6,8,6
shuffle2	3,10,9,10
shuffle2	4,7,3,7
shuffle2	15,5,4,5
shuffle2	17,11,15,11

shuffle4	16,19,17,19
shuffle4	12,18,16,18
shuffle4	6,14,12,14
shuffle4	8,13,6,13
shuffle4    7,10,8,10
shuffle4    3,9,7,9
shuffle4    11,5,3,5
shuffle4    15,4,11,4

shuffle8	14,19,15,19
shuffle8	13,18,14,18
shuffle8	12,17,13,17
shuffle8    6,16,12,16
shuffle8	5,10,6,10
shuffle8	4,9,5,9
shuffle8	3,8,4,8
shuffle8    11,7,3,7

shuffle16   10,19,11,19
shuffle16   9,18,10,18
shuffle16   8,17,9,17
shuffle16   7,16,8,16
shuffle16   6,15,7,15
shuffle16   5,14,6,14
shuffle16   4,13,5,13
shuffle16   3,12,4,12

vmovdqu32		%zmm4,(\q      )*2(%rdi)
vmovdqu32		%zmm5,(\q+   \r)*2(%rdi)
vmovdqu32		%zmm6,(\q+ \r*2)*2(%rdi)
vmovdqu32		%zmm7,(\q+ \r*3)*2(%rdi)
vmovdqu32		%zmm8,(\q+ \r*4)*2(%rdi)
vmovdqu32		%zmm9,(\q+ \r*5)*2(%rdi)
vmovdqu32		%zmm10,(\q+ \r*6)*2(%rdi)
vmovdqu32		%zmm11,(\q+ \r*7)*2(%rdi)
vmovdqu32		%zmm12,(\q+ \r*8)*2(%rdi)
vmovdqu32		%zmm13,(\q+ \r*9)*2(%rdi)
vmovdqu32		%zmm14,(\q+ \r*10)*2(%rdi)
vmovdqu32		%zmm15,(\q+ \r*11)*2(%rdi)
vmovdqu32		%zmm16,(\q+ \r*12)*2(%rdi)
vmovdqu32		%zmm17,(\q+ \r*13)*2(%rdi)
vmovdqu32		%zmm18,(\q+ \r*14)*2(%rdi)
vmovdqu32		%zmm19,(\q+ \r*15)*2(%rdi)

.endm


.text

.global cdecl(poly_formseqto32_AVX512)
cdecl(poly_formseqto32_AVX512):
kmovw       K1IDX*2(%rcx),%k1 
kmovw       K2IDX*2(%rcx),%k2

# the first part of to32 conversion
formseqto32_AVX512_256b 0 0 256
formseqto32_AVX512_256b 4096 32 256 # 因为AVX512寄存器不够，不能32个全部用来存储数据，所以每次只能完成一半的转换。因此，宏后面的数字要这样写，即两行为一组。此外，这个宏运行完，并没有完全变成32-way格式，而是a(0)0~a(15)0,a(0)1~a(15)1,a(16)0~a(31)0, a(16)1~a(31)1,...,a(0)254~a(15)254,a(0)255~a(15)255,a(16)254~a(31)254,a(16)255~a(31)255,a(0)表示第0-way，a(1)表示第1-way，以此类推。因此还需要额外的一个宏定义formseq_supp才可以完成整体转换。
formseqto32_AVX512_256b 32 1024 256
formseqto32_AVX512_256b 4128 1056 256
formseqto32_AVX512_256b 64 2048 256
formseqto32_AVX512_256b 4160 2080 256
formseqto32_AVX512_256b 96 3072 256
formseqto32_AVX512_256b 4192 3104 256

formseqto32_AVX512_256b 128 4096 256
formseqto32_AVX512_256b 4224 4128 256
formseqto32_AVX512_256b 160 5120 256
formseqto32_AVX512_256b 4256 5152 256
formseqto32_AVX512_256b 192 6144 256
formseqto32_AVX512_256b 4288 6176 256
formseqto32_AVX512_256b 224 7168 256
formseqto32_AVX512_256b 4320 7200 256

# the second part of to32 conversion
formseq_supp 0 0  #每个formseq_supp转换30*32个系数
formseq_supp 960 960
formseq_supp 1920 1920
formseq_supp 2880 2880
formseq_supp 3840 3840
formseq_supp 4800 4800
formseq_supp 5760 5760
formseq_supp 6720 6720

vmovdqu32		(7680     )*2(%rsi),%zmm1 # 剩下的16*32个系数不够再调用一次formseq_supp，所以直接单独写
vmovdqu32		(7680+  32)*2(%rsi),%zmm2
vmovdqu32		(7680+  64)*2(%rsi),%zmm3
vmovdqu32		(7680+  96)*2(%rsi),%zmm4
vmovdqu32		(7680+ 128)*2(%rsi),%zmm5
vmovdqu32		(7680+ 160)*2(%rsi),%zmm6
vmovdqu32		(7680+ 192)*2(%rsi),%zmm7
vmovdqu32		(7680+ 224)*2(%rsi),%zmm8
vmovdqu32		(7680+ 256)*2(%rsi),%zmm9
vmovdqu32		(7680+ 288)*2(%rsi),%zmm10
vmovdqu32		(7680+ 320)*2(%rsi),%zmm11
vmovdqu32		(7680+ 352)*2(%rsi),%zmm12
vmovdqu32		(7680+ 384)*2(%rsi),%zmm13
vmovdqu32		(7680+ 416)*2(%rsi),%zmm14
vmovdqu32		(7680+ 448)*2(%rsi),%zmm15
vmovdqu32		(7680+ 480)*2(%rsi),%zmm16

shuffle16   1,2,0,1
shuffle16   3,4,2,3
shuffle16   5,6,4,5
shuffle16   7,8,6,7
shuffle16   9,10,8,9
shuffle16   11,12,10,11
shuffle16   13,14,12,13
shuffle16   15,16,14,15

vmovdqu32		%zmm0,(7680     )*2(%rdx)
vmovdqu32		%zmm1,(7680+  32)*2(%rdx)
vmovdqu32		%zmm2,(7680+  64)*2(%rdx)
vmovdqu32		%zmm3,(7680+  96)*2(%rdx)
vmovdqu32		%zmm4,(7680+ 128)*2(%rdx)
vmovdqu32		%zmm5,(7680+ 160)*2(%rdx)
vmovdqu32		%zmm6,(7680+ 192)*2(%rdx)
vmovdqu32		%zmm7,(7680+ 224)*2(%rdx)
vmovdqu32		%zmm8,(7680+ 256)*2(%rdx)
vmovdqu32		%zmm9,(7680+ 288)*2(%rdx)
vmovdqu32		%zmm10,(7680+ 320)*2(%rdx)
vmovdqu32		%zmm11,(7680+ 352)*2(%rdx)
vmovdqu32		%zmm12,(7680+ 384)*2(%rdx)
vmovdqu32		%zmm13,(7680+ 416)*2(%rdx)
vmovdqu32		%zmm14,(7680+ 448)*2(%rdx)
vmovdqu32		%zmm15,(7680+ 480)*2(%rdx)

ret


.global cdecl(pk_formseqfrom32_AVX512)
cdecl(pk_formseqfrom32_AVX512):
kmovw       K1IDX*2(%rcx),%k1 
kmovw       K2IDX*2(%rcx),%k2

# the first part of from32 conversion
formseq_supp 0 0  #每个formseq_supp转换30*32个系数
formseq_supp 960 960
formseq_supp 1920 1920
formseq_supp 2880 2880
formseq_supp 3840 3840
formseq_supp 4800 4800
formseq_supp 5760 5760
formseq_supp 6720 6720
formseq_supp 7680 7680
formseq_supp 8640 8640

formseq_supp 9600 9600
formseq_supp 10560 10560
formseq_supp 11520 11520
formseq_supp 12480 12480
formseq_supp 13440 13440
formseq_supp 14400 14400
formseq_supp 15360 15360
formseq_supp 16320 16320
formseq_supp 17280 17280

vmovdqu32		(18240     )*2(%rsi),%zmm1 # 剩下的6*32个系数不够再调用一次formseq_supp，所以直接单独写
vmovdqu32		(18240+  32)*2(%rsi),%zmm2
vmovdqu32		(18240+  64)*2(%rsi),%zmm3
vmovdqu32		(18240+  96)*2(%rsi),%zmm4
vmovdqu32		(18240+ 128)*2(%rsi),%zmm5
vmovdqu32		(18240+ 160)*2(%rsi),%zmm6

shuffle16    1,2,0,1
shuffle16    3,4,2,3
shuffle16    5,6,4,5

vmovdqu32		%zmm0,(18240     )*2(%rdx)
vmovdqu32		%zmm1,(18240+  32)*2(%rdx)
vmovdqu32		%zmm2,(18240+  64)*2(%rdx)
vmovdqu32		%zmm3,(18240+  96)*2(%rdx)
vmovdqu32		%zmm4,(18240+ 128)*2(%rdx)
vmovdqu32		%zmm5,(18240+ 160)*2(%rdx)

# the second part of from32 conversion
formseqfrom32_AVX512_256b 0 0 592
formseqfrom32_AVX512_256b 32 9472 592  //9472 = (3*384/2+32/2)*16 = 592*16
formseqfrom32_AVX512_256b 1024 32 592
formseqfrom32_AVX512_256b 1056 9504 592
formseqfrom32_AVX512_256b 2048 64 592
formseqfrom32_AVX512_256b 2080 9536 592
formseqfrom32_AVX512_256b 3072 96 592
formseqfrom32_AVX512_256b 3104 9568 592

formseqfrom32_AVX512_256b 4096 128 592
formseqfrom32_AVX512_256b 4128 9600 592
formseqfrom32_AVX512_256b 5120 160 592
formseqfrom32_AVX512_256b 5152 9632 592
formseqfrom32_AVX512_256b 6144 192 592
formseqfrom32_AVX512_256b 6176 9664 592
formseqfrom32_AVX512_256b 7168 224 592
formseqfrom32_AVX512_256b 7200 9696 592

formseqfrom32_AVX512_256b 8192 256 592
formseqfrom32_AVX512_256b 8224 9728 592
formseqfrom32_AVX512_256b 9216 288 592
formseqfrom32_AVX512_256b 9248 9760 592
formseqfrom32_AVX512_256b 10240 320 592
formseqfrom32_AVX512_256b 10272 9792 592
formseqfrom32_AVX512_256b 11264 352 592
formseqfrom32_AVX512_256b 11296 9824 592

formseqfrom32_AVX512_256b 12288 384 592
formseqfrom32_AVX512_256b 12320 9856 592
formseqfrom32_AVX512_256b 13312 416 592
formseqfrom32_AVX512_256b 13344 9888 592
formseqfrom32_AVX512_256b 14336 448 592
formseqfrom32_AVX512_256b 14368 9920 592
formseqfrom32_AVX512_256b 15360 480 592
formseqfrom32_AVX512_256b 15392 9952 592

formseqfrom32_AVX512_256b 16384 512 592
formseqfrom32_AVX512_256b 16416 9984 592
formseqfrom32_AVX512_256b 17408 544 592
formseqfrom32_AVX512_256b 17440 10016 592

ret


.global cdecl(pk_formseqto32_AVX512)
cdecl(pk_formseqto32_AVX512):
kmovw       K1IDX*2(%rcx),%k1 
kmovw       K2IDX*2(%rcx),%k2

# the first part of to32 conversion
formseqto32_AVX512_256b 0 0 592
formseqto32_AVX512_256b 9472 32 592  //9472 = (3*384/2+32/2)*16 = 592*16
formseqto32_AVX512_256b 32 1024 592
formseqto32_AVX512_256b 9504 1056 592
formseqto32_AVX512_256b 64 2048 592
formseqto32_AVX512_256b 9536 2080 592
formseqto32_AVX512_256b 96 3072 592
formseqto32_AVX512_256b 9568 3104 592

formseqto32_AVX512_256b 128 4096 592
formseqto32_AVX512_256b 9600 4128 592
formseqto32_AVX512_256b 160 5120 592
formseqto32_AVX512_256b 9632 5152 592
formseqto32_AVX512_256b 192 6144 592
formseqto32_AVX512_256b 9664 6176 592
formseqto32_AVX512_256b 224 7168 592
formseqto32_AVX512_256b 9696 7200 592

formseqto32_AVX512_256b 256 8192 592
formseqto32_AVX512_256b 9728 8224 592
formseqto32_AVX512_256b 288 9216 592
formseqto32_AVX512_256b 9760 9248 592
formseqto32_AVX512_256b 320 10240 592
formseqto32_AVX512_256b 9792 10272 592
formseqto32_AVX512_256b 352 11264 592
formseqto32_AVX512_256b 9824 11296 592

formseqto32_AVX512_256b 384 12288 592
formseqto32_AVX512_256b 9856 12320 592
formseqto32_AVX512_256b 416 13312 592
formseqto32_AVX512_256b 9888 13344 592
formseqto32_AVX512_256b 448 14336 592
formseqto32_AVX512_256b 9920 14368 592
formseqto32_AVX512_256b 480 15360 592
formseqto32_AVX512_256b 9952 15392 592

formseqto32_AVX512_256b 512 16384 592
formseqto32_AVX512_256b 9984 16416 592
formseqto32_AVX512_256b 544 17408 592
formseqto32_AVX512_256b 10016 17440 592


# the second part of to32 conversion
formseq_supp 0 0  #每个formseq_supp转换30*32个系数
formseq_supp 960 960
formseq_supp 1920 1920
formseq_supp 2880 2880
formseq_supp 3840 3840
formseq_supp 4800 4800
formseq_supp 5760 5760
formseq_supp 6720 6720
formseq_supp 7680 7680
formseq_supp 8640 8640

formseq_supp 9600 9600
formseq_supp 10560 10560
formseq_supp 11520 11520
formseq_supp 12480 12480
formseq_supp 13440 13440
formseq_supp 14400 14400
formseq_supp 15360 15360
formseq_supp 16320 16320
formseq_supp 17280 17280

vmovdqu32		(18240     )*2(%rsi),%zmm1 # 剩下的6*32个系数不够再调用一次formseq_supp，所以直接单独写
vmovdqu32		(18240+  32)*2(%rsi),%zmm2
vmovdqu32		(18240+  64)*2(%rsi),%zmm3
vmovdqu32		(18240+  96)*2(%rsi),%zmm4
vmovdqu32		(18240+ 128)*2(%rsi),%zmm5
vmovdqu32		(18240+ 160)*2(%rsi),%zmm6

shuffle16    1,2,0,1
shuffle16    3,4,2,3
shuffle16    5,6,4,5

vmovdqu32		%zmm0,(18240     )*2(%rdx)
vmovdqu32		%zmm1,(18240+  32)*2(%rdx)
vmovdqu32		%zmm2,(18240+  64)*2(%rdx)
vmovdqu32		%zmm3,(18240+  96)*2(%rdx)
vmovdqu32		%zmm4,(18240+ 128)*2(%rdx)
vmovdqu32		%zmm5,(18240+ 160)*2(%rdx)

ret


.global cdecl(sk_formseqfrom32_AVX512)
cdecl(sk_formseqfrom32_AVX512):
kmovw       K1IDX*2(%rcx),%k1 
kmovw       K2IDX*2(%rcx),%k2

# the first part of from32 conversion
formseq_supp 0 0  #每个formseq_supp转换30*32个系数
formseq_supp 960 960
formseq_supp 1920 1920
formseq_supp 2880 2880
formseq_supp 3840 3840
formseq_supp 4800 4800
formseq_supp 5760 5760
formseq_supp 6720 6720
formseq_supp 7680 7680
formseq_supp 8640 8640

formseq_supp 9600 9600
formseq_supp 10560 10560
formseq_supp 11520 11520
formseq_supp 12480 12480
formseq_supp 13440 13440
formseq_supp 14400 14400
formseq_supp 15360 15360
formseq_supp 16320 16320
formseq_supp 17280 17280

vmovdqu32		(18240     )*2(%rsi),%zmm1 # 剩下的6*32个系数不够再调用一次formseq_supp，所以直接单独写
vmovdqu32		(18240+  32)*2(%rsi),%zmm2
vmovdqu32		(18240+  64)*2(%rsi),%zmm3
vmovdqu32		(18240+  96)*2(%rsi),%zmm4
vmovdqu32		(18240+ 128)*2(%rsi),%zmm5
vmovdqu32		(18240+ 160)*2(%rsi),%zmm6

shuffle16    1,2,0,1
shuffle16    3,4,2,3
shuffle16    5,6,4,5

vmovdqu32		%zmm0,(18240     )*2(%rdx)
vmovdqu32		%zmm1,(18240+  32)*2(%rdx)
vmovdqu32		%zmm2,(18240+  64)*2(%rdx)
vmovdqu32		%zmm3,(18240+  96)*2(%rdx)
vmovdqu32		%zmm4,(18240+ 128)*2(%rdx)
vmovdqu32		%zmm5,(18240+ 160)*2(%rdx)

# the second part of from32 conversion
formseqfrom32_AVX512_256b 0 0 576
formseqfrom32_AVX512_256b 32 9216 576  // 9216 = 576*16
formseqfrom32_AVX512_256b 1024 32 576
formseqfrom32_AVX512_256b 1056 9248 576
formseqfrom32_AVX512_256b 2048 64 576
formseqfrom32_AVX512_256b 2080 9280 576
formseqfrom32_AVX512_256b 3072 96 576
formseqfrom32_AVX512_256b 3104 9312 576

formseqfrom32_AVX512_256b 4096 128 576
formseqfrom32_AVX512_256b 4128 9344 576
formseqfrom32_AVX512_256b 5120 160 576
formseqfrom32_AVX512_256b 5152 9376 576
formseqfrom32_AVX512_256b 6144 192 576
formseqfrom32_AVX512_256b 6176 9408 576
formseqfrom32_AVX512_256b 7168 224 576
formseqfrom32_AVX512_256b 7200 9440 576

formseqfrom32_AVX512_256b 8192 256 576
formseqfrom32_AVX512_256b 8224 9472 576
formseqfrom32_AVX512_256b 9216 288 576
formseqfrom32_AVX512_256b 9248 9504 576
formseqfrom32_AVX512_256b 10240 320 576
formseqfrom32_AVX512_256b 10272 9536 576
formseqfrom32_AVX512_256b 11264 352 576
formseqfrom32_AVX512_256b 11296 9568 576

formseqfrom32_AVX512_256b 12288 384 576
formseqfrom32_AVX512_256b 12320 9600 576
formseqfrom32_AVX512_256b 13312 416 576
formseqfrom32_AVX512_256b 13344 9632 576
formseqfrom32_AVX512_256b 14336 448 576
formseqfrom32_AVX512_256b 14368 9664 576
formseqfrom32_AVX512_256b 15360 480 576
formseqfrom32_AVX512_256b 15392 9696 576

formseqfrom32_AVX512_256b 16384 512 576
formseqfrom32_AVX512_256b 16416 9728 576
formseqfrom32_AVX512_256b 17408 544 576
formseqfrom32_AVX512_256b 17440 9760 576

ret


.global cdecl(sk_formseqto32_AVX512)
cdecl(sk_formseqto32_AVX512):
kmovw       K1IDX*2(%rcx),%k1 
kmovw       K2IDX*2(%rcx),%k2

# the first part of to32 conversion
formseqto32_AVX512_256b 0 0 576
formseqto32_AVX512_256b 9216 32 576  // 9216 = 576*16
formseqto32_AVX512_256b 32 1024 576
formseqto32_AVX512_256b 9248 1056 576
formseqto32_AVX512_256b 64 2048 576
formseqto32_AVX512_256b 9280 2080 576
formseqto32_AVX512_256b 96 3072 576
formseqto32_AVX512_256b 9312 3104 576

formseqto32_AVX512_256b 128 4096 576
formseqto32_AVX512_256b 9344 4128 576
formseqto32_AVX512_256b 160 5120 576
formseqto32_AVX512_256b 9376 5152 576
formseqto32_AVX512_256b 192 6144 576
formseqto32_AVX512_256b 9408 6176 576
formseqto32_AVX512_256b 224 7168 576
formseqto32_AVX512_256b 9440 7200 576

formseqto32_AVX512_256b 256 8192 576
formseqto32_AVX512_256b 9472 8224 576
formseqto32_AVX512_256b 288 9216 576
formseqto32_AVX512_256b 9504 9248 576
formseqto32_AVX512_256b 320 10240 576
formseqto32_AVX512_256b 9536 10272 576
formseqto32_AVX512_256b 352 11264 576
formseqto32_AVX512_256b 9568 11296 576

formseqto32_AVX512_256b 384 12288 576
formseqto32_AVX512_256b 9600 12320 576
formseqto32_AVX512_256b 416 13312 576
formseqto32_AVX512_256b 9632 13344 576
formseqto32_AVX512_256b 448 14336 576
formseqto32_AVX512_256b 9664 14368 576
formseqto32_AVX512_256b 480 15360 576
formseqto32_AVX512_256b 9696 15392 576

formseqto32_AVX512_256b 512 16384 576
formseqto32_AVX512_256b 9728 16416 576
formseqto32_AVX512_256b 544 17408 576
formseqto32_AVX512_256b 9760 17440 576


# the second part of to32 conversion
formseq_supp 0 0  #每个formseq_supp转换30*32个系数
formseq_supp 960 960
formseq_supp 1920 1920
formseq_supp 2880 2880
formseq_supp 3840 3840
formseq_supp 4800 4800
formseq_supp 5760 5760
formseq_supp 6720 6720
formseq_supp 7680 7680
formseq_supp 8640 8640

formseq_supp 9600 9600
formseq_supp 10560 10560
formseq_supp 11520 11520
formseq_supp 12480 12480
formseq_supp 13440 13440
formseq_supp 14400 14400
formseq_supp 15360 15360
formseq_supp 16320 16320
formseq_supp 17280 17280

vmovdqu32		(18240     )*2(%rsi),%zmm1 # 剩下的6*32个系数不够再调用一次formseq_supp，所以直接单独写
vmovdqu32		(18240+  32)*2(%rsi),%zmm2
vmovdqu32		(18240+  64)*2(%rsi),%zmm3
vmovdqu32		(18240+  96)*2(%rsi),%zmm4
vmovdqu32		(18240+ 128)*2(%rsi),%zmm5
vmovdqu32		(18240+ 160)*2(%rsi),%zmm6

shuffle16    1,2,0,1
shuffle16    3,4,2,3
shuffle16    5,6,4,5

vmovdqu32		%zmm0,(18240     )*2(%rdx)
vmovdqu32		%zmm1,(18240+  32)*2(%rdx)
vmovdqu32		%zmm2,(18240+  64)*2(%rdx)
vmovdqu32		%zmm3,(18240+  96)*2(%rdx)
vmovdqu32		%zmm4,(18240+ 128)*2(%rdx)
vmovdqu32		%zmm5,(18240+ 160)*2(%rdx)

ret


.global cdecl(cipher_formseqto32_AVX512)
cdecl(cipher_formseqto32_AVX512):
kmovw       K1IDX*2(%rcx),%k1 
kmovw       K2IDX*2(%rcx),%k2

# the first part of to32 conversion
formseqto32_AVX512_256b 0 0 544  # 544 = 3*320/2 + 128/2
formseqto32_AVX512_256b 8704 32 544  # 8704 = 544*16
formseqto32_AVX512_256b 32 1024 544
formseqto32_AVX512_256b 8736 1056 544
formseqto32_AVX512_256b 64 2048 544
formseqto32_AVX512_256b 8768 2080 544
formseqto32_AVX512_256b 96 3072 544
formseqto32_AVX512_256b 8800 3104 544
formseqto32_AVX512_256b 128 4096 544
formseqto32_AVX512_256b 8832 4128 544

formseqto32_AVX512_256b 160 5120 544
formseqto32_AVX512_256b 8864 5152 544
formseqto32_AVX512_256b 192 6144 544
formseqto32_AVX512_256b 8896 6176 544
formseqto32_AVX512_256b 224 7168 544
formseqto32_AVX512_256b 8928 7200 544
formseqto32_AVX512_256b 256 8192 544
formseqto32_AVX512_256b 8960 8224 544
formseqto32_AVX512_256b 288 9216 544
formseqto32_AVX512_256b 8992 9248 544

formseqto32_AVX512_256b 320 10240 544
formseqto32_AVX512_256b 9024 10272 544
formseqto32_AVX512_256b 352 11264 544
formseqto32_AVX512_256b 9056 11296 544
formseqto32_AVX512_256b 384 12288 544
formseqto32_AVX512_256b 9088 12320 544
formseqto32_AVX512_256b 416 13312 544
formseqto32_AVX512_256b 9120 13344 544
formseqto32_AVX512_256b 448 14336 544
formseqto32_AVX512_256b 9152 14368 544  # 3*160*32*2*8B结束

########################################

formseqto32_AVX512_256b 480 15360 544
formseqto32_AVX512_256b 9184 15392 544
formseqto32_AVX512_256b 512 16384 544
formseqto32_AVX512_256b 9216 16416 544

# the second part of to32 conversion
formseq_supp 0 0  #每个formseq_supp转换30*32个系数
formseq_supp 960 960
formseq_supp 1920 1920
formseq_supp 2880 2880
formseq_supp 3840 3840
formseq_supp 4800 4800
formseq_supp 5760 5760
formseq_supp 6720 6720
formseq_supp 7680 7680
formseq_supp 8640 8640

formseq_supp 9600 9600
formseq_supp 10560 10560
formseq_supp 11520 11520
formseq_supp 12480 12480
formseq_supp 13440 13440
formseq_supp 14400 14400
formseq_supp 15360 15360
formseq_supp 16320 16320

vmovdqu32		(17280     )*2(%rsi),%zmm1
vmovdqu32		(17280+  32)*2(%rsi),%zmm2
vmovdqu32		(17280+  64)*2(%rsi),%zmm3
vmovdqu32		(17280+  96)*2(%rsi),%zmm4

shuffle16    1,2,0,1
shuffle16    3,4,2,3

vmovdqu32		%zmm0,(17280     )*2(%rdx)
vmovdqu32		%zmm1,(17280+  32)*2(%rdx)
vmovdqu32		%zmm2,(17280+  64)*2(%rdx)
vmovdqu32		%zmm3,(17280+  96)*2(%rdx)

ret


.global cdecl(cipher_formseqfrom32_AVX512)
cdecl(cipher_formseqfrom32_AVX512):
kmovw       K1IDX*2(%rcx),%k1 
kmovw       K2IDX*2(%rcx),%k2

# the first part of from32 conversion
formseq_supp 0 0  #每个formseq_supp转换30*32个系数
formseq_supp 960 960
formseq_supp 1920 1920
formseq_supp 2880 2880
formseq_supp 3840 3840
formseq_supp 4800 4800
formseq_supp 5760 5760
formseq_supp 6720 6720
formseq_supp 7680 7680
formseq_supp 8640 8640

formseq_supp 9600 9600
formseq_supp 10560 10560
formseq_supp 11520 11520
formseq_supp 12480 12480
formseq_supp 13440 13440
formseq_supp 14400 14400
formseq_supp 15360 15360
formseq_supp 16320 16320

vmovdqu32		(17280     )*2(%rsi),%zmm1
vmovdqu32		(17280+  32)*2(%rsi),%zmm2
vmovdqu32		(17280+  64)*2(%rsi),%zmm3
vmovdqu32		(17280+  96)*2(%rsi),%zmm4

shuffle16    1,2,0,1
shuffle16    3,4,2,3

vmovdqu32		%zmm0,(17280     )*2(%rdx)
vmovdqu32		%zmm1,(17280+  32)*2(%rdx)
vmovdqu32		%zmm2,(17280+  64)*2(%rdx)
vmovdqu32		%zmm3,(17280+  96)*2(%rdx)

# the second part of from32 conversion
formseqfrom32_AVX512_256b 0 0 544
formseqfrom32_AVX512_256b 32 8704 544  // 8704 = 544*16
formseqfrom32_AVX512_256b 1024 32 544
formseqfrom32_AVX512_256b 1056 8736 544
formseqfrom32_AVX512_256b 2048 64 544
formseqfrom32_AVX512_256b 2080 8768 544
formseqfrom32_AVX512_256b 3072 96 544
formseqfrom32_AVX512_256b 3104 8800 544

formseqfrom32_AVX512_256b 4096 128 544
formseqfrom32_AVX512_256b 4128 8832 544
formseqfrom32_AVX512_256b 5120 160 544
formseqfrom32_AVX512_256b 5152 8864 544
formseqfrom32_AVX512_256b 6144 192 544
formseqfrom32_AVX512_256b 6176 8896 544
formseqfrom32_AVX512_256b 7168 224 544
formseqfrom32_AVX512_256b 7200 8928 544

formseqfrom32_AVX512_256b 8192 256 544
formseqfrom32_AVX512_256b 8224 8960 544
formseqfrom32_AVX512_256b 9216 288 544
formseqfrom32_AVX512_256b 9248 8992 544
formseqfrom32_AVX512_256b 10240 320 544
formseqfrom32_AVX512_256b 10272 9024 544
formseqfrom32_AVX512_256b 11264 352 544
formseqfrom32_AVX512_256b 11296 9056 544

formseqfrom32_AVX512_256b 12288 384 544
formseqfrom32_AVX512_256b 12320 9088 544
formseqfrom32_AVX512_256b 13312 416 544
formseqfrom32_AVX512_256b 13344 9120 544
formseqfrom32_AVX512_256b 14336 448 544
formseqfrom32_AVX512_256b 14368 9152 544  // 3*160*32*2*8B结束，cipher from32之前是(3*320*32+128*32)个uint8_t，那么到这里就需要3*320*8/ 256 = 30个

########################################

formseqfrom32_AVX512_256b 15360 480 544
formseqfrom32_AVX512_256b 15392 9184 544
formseqfrom32_AVX512_256b 16384 512 544
formseqfrom32_AVX512_256b 16416 9216 544


ret