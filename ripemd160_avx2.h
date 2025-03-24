#ifndef RIPEMD160_AVX2_H
#define RIPEMD160_AVX2_H

#include <immintrin.h>
#include <cstdint>

namespace ripemd160avx2 {

// Initialize Ripemd160
void Initialize(__m256i *state);

// Transform AVX2
void Transform(__m256i *state, uint8_t *blocks[8]);

// Hashing function
void ripemd160avx2_32(
    const unsigned char* i0, const unsigned char* i1,
    const unsigned char* i2, const unsigned char* i3,
    const unsigned char* i4, const unsigned char* i5,
    const unsigned char* i6, const unsigned char* i7,
    unsigned char* d0, unsigned char* d1,
    unsigned char* d2, unsigned char* d3,
    unsigned char* d4, unsigned char* d5,
    unsigned char* d6, unsigned char* d7);

} // namespace ripemd160avx2

#endif  // RIPEMD160_AVX2_H