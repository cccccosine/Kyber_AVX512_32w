.macro red32 r,rs=0,x=20
vpmulhw         %zmm1,%zmm\r,%zmm\x
.if \rs
vpmulhrsw	%zmm\rs,%zmm\x,%zmm\x
.else
vpsraw          $10,%zmm\x,%zmm\x
.endif
vpmullw         %zmm0,%zmm\x,%zmm\x
vpsubw          %zmm\x,%zmm\r,%zmm\r
.endm

.macro csubq r,x=20
vpsubw		%zmm0,%zmm\r,%zmm\r
vpsraw		$15,%zmm\r,%zmm\x
vpandd		%zmm0,%zmm\x,%zmm\x
vpaddw		%zmm\x,%zmm\r,%zmm\r
.endm

.macro caddq r,x=12
vpsraw		$15,%zmm\r,%zmm\x
vpandd		%zmm0,%zmm\x,%zmm\x
vpaddw		%zmm\x,%zmm\r,%zmm\r
.endm

.macro fqmulprecomp al,ah,b,x=20
vpmullw		%zmm\al,%zmm\b,%zmm\x
vpmulhw		%zmm\ah,%zmm\b,%zmm\b
vpmulhw		%zmm0,%zmm\x,%zmm\x
vpsubw		%zmm\x,%zmm\b,%zmm\b
.endm
