#include "ripemd160_avx2.h"
#include <immintrin.h>
#include <cstring>
#include <cstdint>

namespace ripemd160avx2 {

// Aligned initial state constants
#ifdef _MSC_VER
__declspec(align(32)) static const uint32_t _init[40] = {
#else
static const uint32_t _init[40] __attribute__((aligned(32))) = {
#endif
    // 8 copies of each initial state variable (A-E)
    0x67452301, 0x67452301, 0x67452301, 0x67452301, 0x67452301, 0x67452301, 0x67452301, 0x67452301, // A
    0xEFCDAB89, 0xEFCDAB89, 0xEFCDAB89, 0xEFCDAB89, 0xEFCDAB89, 0xEFCDAB89, 0xEFCDAB89, 0xEFCDAB89, // B
    0x98BADCFE, 0x98BADCFE, 0x98BADCFE, 0x98BADCFE, 0x98BADCFE, 0x98BADCFE, 0x98BADCFE, 0x98BADCFE, // C
    0x10325476, 0x10325476, 0x10325476, 0x10325476, 0x10325476, 0x10325476, 0x10325476, 0x10325476, // D
    0xC3D2E1F0, 0xC3D2E1F0, 0xC3D2E1F0, 0xC3D2E1F0, 0xC3D2E1F0, 0xC3D2E1F0, 0xC3D2E1F0, 0xC3D2E1F0  // E
};

// Helper macros
#define _mm256_not_si256(x) _mm256_xor_si256((x), _mm256_set1_epi32(-1))
#define ROL(x, n) _mm256_or_si256(_mm256_slli_epi32((x), (n)), _mm256_srli_epi32((x), 32-(n)))
#define add3(x, y, z) _mm256_add_epi32(_mm256_add_epi32((x), (y)), (z))
#define add4(x, y, z, w) _mm256_add_epi32(_mm256_add_epi32((x), (y)), _mm256_add_epi32((z), (w)))

// RIPEMD-160 functions
#define f1(x, y, z) _mm256_xor_si256((x), _mm256_xor_si256((y), (z)))
#define f2(x, y, z) _mm256_or_si256(_mm256_and_si256((x), (y)), _mm256_andnot_si256((x), (z)))
#define f3(x, y, z) _mm256_xor_si256(_mm256_or_si256((x), _mm256_not_si256((y))), (z))
#define f4(x, y, z) _mm256_or_si256(_mm256_and_si256((x), (z)), _mm256_andnot_si256((z), (y)))
#define f5(x, y, z) _mm256_xor_si256((x), _mm256_or_si256((y), _mm256_not_si256((z))))

// Round macros with forced inlining
#define Round(a, b, c, d, e, f, x, k, r) do { \
    __m256i u = add4((a), (f), (x), _mm256_set1_epi32(k)); \
    (a) = _mm256_add_epi32(ROL(u, (r)), (e)); \
    (c) = ROL((c), 10); \
} while(0)

// Round function macros
#define R11(a, b, c, d, e, x, r) Round(a, b, c, d, e, f1(b, c, d), x, 0, r)
#define R21(a, b, c, d, e, x, r) Round(a, b, c, d, e, f2(b, c, d), x, 0x5A827999, r)
#define R31(a, b, c, d, e, x, r) Round(a, b, c, d, e, f3(b, c, d), x, 0x6ED9EBA1, r)
#define R41(a, b, c, d, e, x, r) Round(a, b, c, d, e, f4(b, c, d), x, 0x8F1BBCDC, r)
#define R51(a, b, c, d, e, x, r) Round(a, b, c, d, e, f5(b, c, d), x, 0xA953FD4E, r)
#define R12(a, b, c, d, e, x, r) Round(a, b, c, d, e, f5(b, c, d), x, 0x50A28BE6, r)
#define R22(a, b, c, d, e, x, r) Round(a, b, c, d, e, f4(b, c, d), x, 0x5C4DD124, r)
#define R32(a, b, c, d, e, x, r) Round(a, b, c, d, e, f3(b, c, d), x, 0x6D703EF3, r)
#define R42(a, b, c, d, e, x, r) Round(a, b, c, d, e, f2(b, c, d), x, 0x7A6D76E9, r)
#define R52(a, b, c, d, e, x, r) Round(a, b, c, d, e, f1(b, c, d), x, 0, r)

// Load message words efficiently
#define LOADW(i) _mm256_setr_epi32( \
    *reinterpret_cast<const uint32_t*>(blk[0] + (i)*4), \
    *reinterpret_cast<const uint32_t*>(blk[1] + (i)*4), \
    *reinterpret_cast<const uint32_t*>(blk[2] + (i)*4), \
    *reinterpret_cast<const uint32_t*>(blk[3] + (i)*4), \
    *reinterpret_cast<const uint32_t*>(blk[4] + (i)*4), \
    *reinterpret_cast<const uint32_t*>(blk[5] + (i)*4), \
    *reinterpret_cast<const uint32_t*>(blk[6] + (i)*4), \
    *reinterpret_cast<const uint32_t*>(blk[7] + (i)*4))

// Initialize state from precomputed constants
inline void Initialize(__m256i* s) {
    _mm256_store_si256(s + 0, _mm256_load_si256(reinterpret_cast<const __m256i*>(_init + 0)));
    _mm256_store_si256(s + 1, _mm256_load_si256(reinterpret_cast<const __m256i*>(_init + 8)));
    _mm256_store_si256(s + 2, _mm256_load_si256(reinterpret_cast<const __m256i*>(_init + 16)));
    _mm256_store_si256(s + 3, _mm256_load_si256(reinterpret_cast<const __m256i*>(_init + 24)));
    _mm256_store_si256(s + 4, _mm256_load_si256(reinterpret_cast<const __m256i*>(_init + 32)));
}

// Process one block for each of 8 messages
inline void Transform(__m256i* s, const uint8_t* blk[8]) {
    __m256i a1 = _mm256_load_si256(s + 0);
    __m256i b1 = _mm256_load_si256(s + 1);
    __m256i c1 = _mm256_load_si256(s + 2);
    __m256i d1 = _mm256_load_si256(s + 3);
    __m256i e1 = _mm256_load_si256(s + 4);
    
    __m256i a2 = a1, b2 = b1, c2 = c1, d2 = d1, e2 = e1;
    
    // Load all message words upfront
    __m256i w[16];
    #ifdef __clang__
    #pragma unroll(16)
    #elif defined(__GNUC__)
    #pragma GCC unroll 16
    #endif
    for (int i = 0; i < 16; ++i) {
        w[i] = LOADW(i);
    }

    // Round 1-16
    R11(a1, b1, c1, d1, e1, w[ 0], 11); R12(a2, b2, c2, d2, e2, w[ 5],  8);
    R11(e1, a1, b1, c1, d1, w[ 1], 14); R12(e2, a2, b2, c2, d2, w[14],  9);
    R11(d1, e1, a1, b1, c1, w[ 2], 15); R12(d2, e2, a2, b2, c2, w[ 7],  9);
    R11(c1, d1, e1, a1, b1, w[ 3], 12); R12(c2, d2, e2, a2, b2, w[ 0], 11);
    R11(b1, c1, d1, e1, a1, w[ 4],  5); R12(b2, c2, d2, e2, a2, w[ 9], 13);
    R11(a1, b1, c1, d1, e1, w[ 5],  8); R12(a2, b2, c2, d2, e2, w[ 2], 15);
    R11(e1, a1, b1, c1, d1, w[ 6],  7); R12(e2, a2, b2, c2, d2, w[11], 15);
    R11(d1, e1, a1, b1, c1, w[ 7],  9); R12(d2, e2, a2, b2, c2, w[ 4],  5);
    R11(c1, d1, e1, a1, b1, w[ 8], 11); R12(c2, d2, e2, a2, b2, w[13],  7);
    R11(b1, c1, d1, e1, a1, w[ 9], 13); R12(b2, c2, d2, e2, a2, w[ 6],  7);
    R11(a1, b1, c1, d1, e1, w[10], 14); R12(a2, b2, c2, d2, e2, w[15],  8);
    R11(e1, a1, b1, c1, d1, w[11], 15); R12(e2, a2, b2, c2, d2, w[ 8], 11);
    R11(d1, e1, a1, b1, c1, w[12],  6); R12(d2, e2, a2, b2, c2, w[ 1], 14);
    R11(c1, d1, e1, a1, b1, w[13],  7); R12(c2, d2, e2, a2, b2, w[10], 14);
    R11(b1, c1, d1, e1, a1, w[14],  9); R12(b2, c2, d2, e2, a2, w[ 3], 12);
    R11(a1, b1, c1, d1, e1, w[15],  8); R12(a2, b2, c2, d2, e2, w[12],  6);

    // Round 17-32
    R21(e1, a1, b1, c1, d1, w[ 7],  7); R22(e2, a2, b2, c2, d2, w[ 6],  9);
    R21(d1, e1, a1, b1, c1, w[ 4],  6); R22(d2, e2, a2, b2, c2, w[11], 13);
    R21(c1, d1, e1, a1, b1, w[13],  8); R22(c2, d2, e2, a2, b2, w[ 3], 15);
    R21(b1, c1, d1, e1, a1, w[ 1], 13); R22(b2, c2, d2, e2, a2, w[ 7],  7);
    R21(a1, b1, c1, d1, e1, w[10], 11); R22(a2, b2, c2, d2, e2, w[ 0], 12);
    R21(e1, a1, b1, c1, d1, w[ 6],  9); R22(e2, a2, b2, c2, d2, w[13],  8);
    R21(d1, e1, a1, b1, c1, w[15],  7); R22(d2, e2, a2, b2, c2, w[ 5],  9);
    R21(c1, d1, e1, a1, b1, w[ 3], 15); R22(c2, d2, e2, a2, b2, w[10], 11);
    R21(b1, c1, d1, e1, a1, w[12],  7); R22(b2, c2, d2, e2, a2, w[14],  7);
    R21(a1, b1, c1, d1, e1, w[ 0], 12); R22(a2, b2, c2, d2, e2, w[15],  7);
    R21(e1, a1, b1, c1, d1, w[ 9], 15); R22(e2, a2, b2, c2, d2, w[ 8], 12);
    R21(d1, e1, a1, b1, c1, w[ 5],  9); R22(d2, e2, a2, b2, c2, w[12],  7);
    R21(c1, d1, e1, a1, b1, w[ 2], 11); R22(c2, d2, e2, a2, b2, w[ 4],  6);
    R21(b1, c1, d1, e1, a1, w[14],  7); R22(b2, c2, d2, e2, a2, w[ 9], 15);
    R21(a1, b1, c1, d1, e1, w[11], 13); R22(a2, b2, c2, d2, e2, w[ 1], 13);
    R21(e1, a1, b1, c1, d1, w[ 8], 12); R22(e2, a2, b2, c2, d2, w[ 2], 11);

    // Round 33-48
    R31(d1, e1, a1, b1, c1, w[ 3], 11); R32(d2, e2, a2, b2, c2, w[15],  9);
    R31(c1, d1, e1, a1, b1, w[10], 13); R32(c2, d2, e2, a2, b2, w[ 5],  7);
    R31(b1, c1, d1, e1, a1, w[14],  6); R32(b2, c2, d2, e2, a2, w[ 1], 15);
    R31(a1, b1, c1, d1, e1, w[ 4],  7); R32(a2, b2, c2, d2, e2, w[ 3], 11);
    R31(e1, a1, b1, c1, d1, w[ 9], 14); R32(e2, a2, b2, c2, d2, w[ 7],  8);
    R31(d1, e1, a1, b1, c1, w[15],  9); R32(d2, e2, a2, b2, c2, w[14],  6);
    R31(c1, d1, e1, a1, b1, w[ 8], 13); R32(c2, d2, e2, a2, b2, w[ 6],  6);
    R31(b1, c1, d1, e1, a1, w[ 1], 15); R32(b2, c2, d2, e2, a2, w[ 9], 14);
    R31(a1, b1, c1, d1, e1, w[ 2], 14); R32(a2, b2, c2, d2, e2, w[11], 12);
    R31(e1, a1, b1, c1, d1, w[ 7],  8); R32(e2, a2, b2, c2, d2, w[ 8], 13);
    R31(d1, e1, a1, b1, c1, w[ 0], 13); R32(d2, e2, a2, b2, c2, w[12],  5);
    R31(c1, d1, e1, a1, b1, w[ 6],  6); R32(c2, d2, e2, a2, b2, w[ 2], 14);
    R31(b1, c1, d1, e1, a1, w[13],  5); R32(b2, c2, d2, e2, a2, w[10], 13);
    R31(a1, b1, c1, d1, e1, w[11], 12); R32(a2, b2, c2, d2, e2, w[ 0], 13);
    R31(e1, a1, b1, c1, d1, w[ 5],  7); R32(e2, a2, b2, c2, d2, w[ 4],  7);
    R31(d1, e1, a1, b1, c1, w[12],  5); R32(d2, e2, a2, b2, c2, w[13],  5);

    // Round 49-64
    R41(c1, d1, e1, a1, b1, w[ 1], 11); R42(c2, d2, e2, a2, b2, w[ 8], 15);
    R41(b1, c1, d1, e1, a1, w[ 9], 12); R42(b2, c2, d2, e2, a2, w[ 6],  5);
    R41(a1, b1, c1, d1, e1, w[11], 14); R42(a2, b2, c2, d2, e2, w[ 4],  8);
    R41(e1, a1, b1, c1, d1, w[10], 15); R42(e2, a2, b2, c2, d2, w[ 1], 11);
    R41(d1, e1, a1, b1, c1, w[ 0], 14); R42(d2, e2, a2, b2, c2, w[ 3], 14);
    R41(c1, d1, e1, a1, b1, w[ 8], 15); R42(c2, d2, e2, a2, b2, w[11], 14);
    R41(b1, c1, d1, e1, a1, w[12],  9); R42(b2, c2, d2, e2, a2, w[15],  6);
    R41(a1, b1, c1, d1, e1, w[ 4],  8); R42(a2, b2, c2, d2, e2, w[ 0], 14);
    R41(e1, a1, b1, c1, d1, w[13],  9); R42(e2, a2, b2, c2, d2, w[ 5],  6);
    R41(d1, e1, a1, b1, c1, w[ 3], 14); R42(d2, e2, a2, b2, c2, w[12],  9);
    R41(c1, d1, e1, a1, b1, w[ 7],  5); R42(c2, d2, e2, a2, b2, w[ 2], 12);
    R41(b1, c1, d1, e1, a1, w[15],  6); R42(b2, c2, d2, e2, a2, w[13],  9);
    R41(a1, b1, c1, d1, e1, w[14],  8); R42(a2, b2, c2, d2, e2, w[ 9], 12);
    R41(e1, a1, b1, c1, d1, w[ 5],  6); R42(e2, a2, b2, c2, d2, w[ 7],  5);
    R41(d1, e1, a1, b1, c1, w[ 6],  5); R42(d2, e2, a2, b2, c2, w[10], 15);
    R41(c1, d1, e1, a1, b1, w[ 2], 12); R42(c2, d2, e2, a2, b2, w[14],  8);

    // Round 65-80
    R51(b1, c1, d1, e1, a1, w[ 4],  9); R52(b2, c2, d2, e2, a2, w[12],  8);
    R51(a1, b1, c1, d1, e1, w[ 0], 15); R52(a2, b2, c2, d2, e2, w[15],  5);
    R51(e1, a1, b1, c1, d1, w[ 5],  5); R52(e2, a2, b2, c2, d2, w[10], 12);
    R51(d1, e1, a1, b1, c1, w[ 9], 11); R52(d2, e2, a2, b2, c2, w[ 4],  9);
    R51(c1, d1, e1, a1, b1, w[ 7],  6); R52(c2, d2, e2, a2, b2, w[ 1], 12);
    R51(b1, c1, d1, e1, a1, w[12],  8); R52(b2, c2, d2, e2, a2, w[ 5],  5);
    R51(a1, b1, c1, d1, e1, w[ 2], 13); R52(a2, b2, c2, d2, e2, w[ 8], 14);
    R51(e1, a1, b1, c1, d1, w[10], 12); R52(e2, a2, b2, c2, d2, w[ 7],  6);
    R51(d1, e1, a1, b1, c1, w[14],  5); R52(d2, e2, a2, b2, c2, w[ 6],  8);
    R51(c1, d1, e1, a1, b1, w[ 1], 12); R52(c2, d2, e2, a2, b2, w[ 2], 13);
    R51(b1, c1, d1, e1, a1, w[ 3], 13); R52(b2, c2, d2, e2, a2, w[13],  6);
    R51(a1, b1, c1, d1, e1, w[ 8], 14); R52(a2, b2, c2, d2, e2, w[14],  5);
    R51(e1, a1, b1, c1, d1, w[11], 11); R52(e2, a2, b2, c2, d2, w[ 0], 15);
    R51(d1, e1, a1, b1, c1, w[ 6],  8); R52(d2, e2, a2, b2, c2, w[ 3], 13);
    R51(c1, d1, e1, a1, b1, w[15],  5); R52(c2, d2, e2, a2, b2, w[ 9], 11);
    R51(b1, c1, d1, e1, a1, w[13],  6); R52(b2, c2, d2, e2, a2, w[11], 11);

    // Combine results
    __m256i t = _mm256_load_si256(s + 0);
    _mm256_store_si256(s + 0, add3(_mm256_load_si256(s + 1), c1, d2));
    _mm256_store_si256(s + 1, add3(_mm256_load_si256(s + 2), d1, e2));
    _mm256_store_si256(s + 2, add3(_mm256_load_si256(s + 3), e1, a2));
    _mm256_store_si256(s + 3, add3(_mm256_load_si256(s + 4), a1, b2));
    _mm256_store_si256(s + 4, add3(t, b1, c2));
}

// Extract results to output buffers
#ifdef _MSC_VER
#define DEPACK(d, i) do { \
    ((uint32_t*)(d))[0] = _mm256_extract_epi32(s[0], (i)); \
    ((uint32_t*)(d))[1] = _mm256_extract_epi32(s[1], (i)); \
    ((uint32_t*)(d))[2] = _mm256_extract_epi32(s[2], (i)); \
    ((uint32_t*)(d))[3] = _mm256_extract_epi32(s[3], (i)); \
    ((uint32_t*)(d))[4] = _mm256_extract_epi32(s[4], (i)); \
} while(0)
#else
#define DEPACK(d, i) do { \
    uint32_t* s0 = (uint32_t*)&s[0]; \
    uint32_t* s1 = (uint32_t*)&s[1]; \
    uint32_t* s2 = (uint32_t*)&s[2]; \
    uint32_t* s3 = (uint32_t*)&s[3]; \
    uint32_t* s4 = (uint32_t*)&s[4]; \
    ((uint32_t*)(d))[0] = s0[(i)]; \
    ((uint32_t*)(d))[1] = s1[(i)]; \
    ((uint32_t*)(d))[2] = s2[(i)]; \
    ((uint32_t*)(d))[3] = s3[(i)]; \
    ((uint32_t*)(d))[4] = s4[(i)]; \
} while(0)
#endif

// Main hash function for 8 32-byte messages
void ripemd160avx2_32(
    const unsigned char* i0, const unsigned char* i1,
    const unsigned char* i2, const unsigned char* i3,
    const unsigned char* i4, const unsigned char* i5,
    const unsigned char* i6, const unsigned char* i7,
    unsigned char* d0, unsigned char* d1,
    unsigned char* d2, unsigned char* d3,
    unsigned char* d4, unsigned char* d5,
    unsigned char* d6, unsigned char* d7)
{
    __m256i s[5];
    const uint8_t* bs[] = {i0, i1, i2, i3, i4, i5, i6, i7};
    
    // Prepare padded blocks (64 bytes each)
    uint8_t blocks[8][64];
    const uint64_t sizedesc = 32 << 3;
    
    #ifdef __clang__
    #pragma unroll(16)
    #elif defined(__GNUC__)
    #pragma GCC unroll 16
    #endif
    for (int i = 0; i < 8; ++i) {
        memcpy(blocks[i], bs[i], 32);
        memset(blocks[i] + 32, 0x80, 1); // Padding
        memset(blocks[i] + 33, 0, 23);
        memcpy(blocks[i] + 56, &sizedesc, 8);
        bs[i] = blocks[i];
    }

    // Process blocks
    Initialize(s);
    Transform(s, bs);

    // Store results
    DEPACK(d0, 7);
    DEPACK(d1, 6);
    DEPACK(d2, 5);
    DEPACK(d3, 4);
    DEPACK(d4, 3);
    DEPACK(d5, 2);
    DEPACK(d6, 1);
    DEPACK(d7, 0);
}

} // namespace ripemd160avx2