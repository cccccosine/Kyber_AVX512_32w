CC ?= /usr/bin/cc
CFLAGS += -g -Wall -Wextra -Wpedantic -Wmissing-prototypes -Wredundant-decls -Wshadow -Wpointer-arith -mavx2 -mavx512f -mavx512bw -mbmi2 -mpopcnt -maes -march=native -mtune=native -O3 -fomit-frame-pointer
# LDFLAGS=-lcrypto

SOURCES= consts_32.c ntt_32.S basemul_32.S cpucycles.c speed_print.c invntt_32.S \
		 cbd.c fips202.c fips202x8.c fq_32.S indcpa_32.c poly_32.c polyvec_32.c \
		 rejsample.c keccak8x/KeccakP-1600-times8-SIMD512.c symmetric-shake.c shuffle_32.S \
		 clocks.c kem_32.c verify_32.c formseq_32.S
		 

HEADERS= consts_32.h ntt_32.h params.h align.h cpucycles.h speed_print.h fq.inc shuffle.inc \
		 reduce.h cbd.h fips202.h fips202x8.h indcpa_32.h poly_32.h polyvec_32.h \
		 randombytes.h rejsample.h symmetric.h clocks.h kem_32.h verify_32.h
		 
all: $(HEADERS) $(SOURCES) main.c randombytes.c 
	$(CC) $(CFLAGS) $(SOURCES) main.c randombytes.c -o main

test_vectors: $(HEADERS) $(SOURCES) test_vectors.c
	$(CC) $(CFLAGS) $(SOURCES) test_vectors.c -o test_vectors

.PHONY: clean

clean:
	-rm all





# ntt : ntt.o
# 	gcc -o ntt ntt.o

# ntt.o : ntt.c
# 	gcc -c ntt.c

# clean : 
# 	rm *.o ntt 



# ntt : ntt.o
# 	gcc -o ntt ntt.o

# ntt.o : ntt.c
# 	gcc -c ntt.c

# clean : 
# 	rm *.o ntt 