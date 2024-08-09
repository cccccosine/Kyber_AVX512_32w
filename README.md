# Kyber_AVX512_32W

The work is the 32-way implementation of cryptographic scheme Kyber on AVX512.

## How to run?

Files:
* main.c : clock cycles and throughput test file
* test_vectors.c : correctness test file

```
#compile and run
make main; ./main  // compile the speed test files and run
make test_vectors; ./test_vectors  //  compile the correctness test files and run
```
