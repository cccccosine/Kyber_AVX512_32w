#include "consts_32.h"
.include "fq.inc"


.text
reduce256_avx:  //reduce256中的256表示一次reduce 32路其中1路的 256bit <==> reduce 32路其中1路的16个系数
#load
vmovdqa32		(%rdi),%zmm2
vmovdqa32		32*2(%rdi),%zmm3  //这里的32*2 = 32*2*8bit，多项式每个系数是16bit，2个字节，所以这里其实是vmovdqa了32个系数进去
vmovdqa32		64*2(%rdi),%zmm4
vmovdqa32		96*2(%rdi),%zmm5
vmovdqa32		128*2(%rdi),%zmm6
vmovdqa32		160*2(%rdi),%zmm7
vmovdqa32		192*2(%rdi),%zmm8
vmovdqa32		224*2(%rdi),%zmm9
vmovdqa32		256*2(%rdi),%zmm10
vmovdqa32		288*2(%rdi),%zmm11
vmovdqa32		320*2(%rdi),%zmm12
vmovdqa32		352*2(%rdi),%zmm13
vmovdqa32		384*2(%rdi),%zmm14
vmovdqa32		416*2(%rdi),%zmm15
vmovdqa32		448*2(%rdi),%zmm16
vmovdqa32		480*2(%rdi),%zmm17

red32		2
red32		3
red32		4
red32		5
red32		6
red32		7
red32		8
red32		9
red32		10
red32		11
red32		12
red32		13
red32		14
red32		15
red32		16
red32		17

#store
vmovdqa32		%zmm2,(%rdi)
vmovdqa32		%zmm3,32*2(%rdi)
vmovdqa32		%zmm4,64*2(%rdi)
vmovdqa32		%zmm5,96*2(%rdi)
vmovdqa32		%zmm6,128*2(%rdi)
vmovdqa32		%zmm7,160*2(%rdi)
vmovdqa32		%zmm8,192*2(%rdi)
vmovdqa32		%zmm9,224*2(%rdi)
vmovdqa32		%zmm10,256*2(%rdi)
vmovdqa32		%zmm11,288*2(%rdi)
vmovdqa32		%zmm12,320*2(%rdi)
vmovdqa32		%zmm13,352*2(%rdi)
vmovdqa32		%zmm14,384*2(%rdi)
vmovdqa32		%zmm15,416*2(%rdi)
vmovdqa32		%zmm16,448*2(%rdi)
vmovdqa32		%zmm17,480*2(%rdi)

ret

.global cdecl(reduce_avx_32)
cdecl(reduce_avx_32):
#consts
vmovdqa32		_32XQ*2(%rsi),%zmm0
vmovdqa32		_32XV*2(%rsi),%zmm1

call		reduce256_avx
add		$1024,%rdi
call		reduce256_avx
add		$1024,%rdi
call		reduce256_avx
add		$1024,%rdi
call		reduce256_avx
add		$1024,%rdi
call		reduce256_avx
add		$1024,%rdi
call		reduce256_avx
add		$1024,%rdi
call		reduce256_avx
add		$1024,%rdi
call		reduce256_avx
add		$1024,%rdi
call		reduce256_avx
add		$1024,%rdi
call		reduce256_avx
add		$1024,%rdi
call		reduce256_avx
add		$1024,%rdi
call		reduce256_avx
add		$1024,%rdi
call		reduce256_avx
add		$1024,%rdi
call		reduce256_avx
add		$1024,%rdi
call		reduce256_avx
add		$1024,%rdi
call		reduce256_avx



ret


tomont256_avx:
#load
vmovdqa32		(%rdi),%zmm3
vmovdqa32		32*2(%rdi),%zmm4
vmovdqa32		64*2(%rdi),%zmm5
vmovdqa32		96*2(%rdi),%zmm6
vmovdqa32		128*2(%rdi),%zmm7
vmovdqa32		160*2(%rdi),%zmm8
vmovdqa32		192*2(%rdi),%zmm9
vmovdqa32		224*2(%rdi),%zmm10
vmovdqa32		256*2(%rdi),%zmm11
vmovdqa32		288*2(%rdi),%zmm12
vmovdqa32		320*2(%rdi),%zmm13
vmovdqa32		352*2(%rdi),%zmm14
vmovdqa32		384*2(%rdi),%zmm15
vmovdqa32		416*2(%rdi),%zmm16
vmovdqa32		448*2(%rdi),%zmm17
vmovdqa32		480*2(%rdi),%zmm18

fqmulprecomp	1,2,3,19
fqmulprecomp	1,2,4,20
fqmulprecomp	1,2,5,21
fqmulprecomp	1,2,6,22
fqmulprecomp	1,2,7,23
fqmulprecomp	1,2,8,24
fqmulprecomp	1,2,9,25
fqmulprecomp	1,2,10,26
fqmulprecomp	1,2,11,27
fqmulprecomp	1,2,12,28
fqmulprecomp	1,2,13,29
fqmulprecomp	1,2,14,30
fqmulprecomp	1,2,15,31
fqmulprecomp	1,2,16,19
fqmulprecomp	1,2,17,20
fqmulprecomp	1,2,18,21


#store
vmovdqa32		%zmm3,(%rdi)
vmovdqa32		%zmm4,32*2(%rdi)
vmovdqa32		%zmm5,64*2(%rdi)
vmovdqa32		%zmm6,96*2(%rdi)
vmovdqa32		%zmm7,128*2(%rdi)
vmovdqa32		%zmm8,160*2(%rdi)
vmovdqa32		%zmm9,192*2(%rdi)
vmovdqa32		%zmm10,224*2(%rdi)
vmovdqa32		%zmm11,256*2(%rdi)
vmovdqa32		%zmm12,288*2(%rdi)
vmovdqa32		%zmm13,320*2(%rdi)
vmovdqa32		%zmm14,352*2(%rdi)
vmovdqa32		%zmm15,384*2(%rdi)
vmovdqa32		%zmm16,416*2(%rdi)
vmovdqa32		%zmm17,448*2(%rdi)
vmovdqa32		%zmm18,480*2(%rdi)

ret

.global cdecl(tomont_avx_32)
cdecl(tomont_avx_32):
#consts
vmovdqa32		_32XQ*2(%rsi),%zmm0
vmovdqa32		_32XMONTSQLO*2(%rsi),%zmm1
vmovdqa32		_32XMONTSQHI*2(%rsi),%zmm2
call		tomont256_avx
add		$1024,%rdi
call		tomont256_avx
add		$1024,%rdi
call		tomont256_avx
add		$1024,%rdi
call		tomont256_avx
add		$1024,%rdi
call		tomont256_avx
add		$1024,%rdi
call		tomont256_avx
add		$1024,%rdi
call		tomont256_avx
add		$1024,%rdi
call		tomont256_avx
add		$1024,%rdi
call		tomont256_avx
add		$1024,%rdi
call		tomont256_avx
add		$1024,%rdi
call		tomont256_avx
add		$1024,%rdi
call		tomont256_avx
add		$1024,%rdi
call		tomont256_avx
add		$1024,%rdi
call		tomont256_avx
add		$1024,%rdi
call		tomont256_avx
add		$1024,%rdi
call		tomont256_avx


ret



