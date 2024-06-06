.macro shuffle16 r0,r1,r2,r3
#vpermd %zmm\r1,%zmm\r512,%zmm\r1
vshufi64x2 $0x44,%zmm\r1,%zmm\r0,%zmm\r2
vshufi64x2 $0xEE,%zmm\r1,%zmm\r0,%zmm\r3
#vpermi2d  %zmm\r0,%zmm\r1,%zmm\r2      
#vpermi2d  %zmm\r0,%zmm\r1,%zmm\r3      
.endm

.macro shuffle8 r0,r1,r2,r3
vpermq    $0x44,%zmm\r1,%zmm\r2
vpblendmd %zmm\r2,%zmm\r0,%zmm\r2   {%k2}
vpermq    $0xEE,%zmm\r0,%zmm\r0
vpblendmd %zmm\r1,%zmm\r0,%zmm\r3   {%k2}
.endm

.macro shuffle4 r0,r1,r2,r3
vshufpd   $0x00,%zmm\r1,%zmm\r0,%zmm\r2
vshufpd   $0xFF,%zmm\r1,%zmm\r0,%zmm\r3      
.endm

.macro shuffle2 r0,r1,r2,r3
vmovsldup	%zmm\r1,%zmm\r2
vpblendmd	%zmm\r2,%zmm\r0,%zmm\r2     {%k1}
vmovshdup	%zmm\r0,%zmm\r0
vpblendmd	%zmm\r1,%zmm\r0,%zmm\r3     {%k1}
.endm

//.macro shuffle1 r0,r1,r2,r3  //目前未用到
//vpslld		$16,%zmm\r1,%zmm\r2
//vpblendmw   %zmm\r2,%zmm\r0,%zmm\r2     {%k3}   //k3 = 0b(0101)*8个
//vpsrld		$16,%zmm\r0,%zmm\r0
//vpblendmw   %zmm\r1,%zmm\r0,%zmm\r3     {%k3}
//.endm
