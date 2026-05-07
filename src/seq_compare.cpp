#include "seq_compare.hpp"
#include "utils.hpp"
#include <cstring>
// DNA 碱基到 4bit 编码的映射表
std::vector<uint64_t> encode4bit(const std::string& s) {
    size_t L = s.size();
    size_t nword = (L + 15) / 16;

    std::vector<uint64_t> out(nword, 0);

    for (size_t i = 0; i < L; ++i) {
        uint8_t code = DNA4BIT[(unsigned char)s[i]];
        size_t w = i >> 4;          // i / 16 word pos
        size_t shift = (i % 16) * 4; // 每个 nt 占 4 bit
        out[w] |= (uint64_t(code) << shift);
    }
    return out;
}


// -------------------------------
// 标量版本
// -------------------------------
bool hamming_k_scalar(const uint64_t* A, const uint64_t* B,
    size_t nword, size_t max_diff)
{
    size_t mism = 0;
    for (size_t i = 0; i < nword; ++i) {
        uint64_t diff = A[i] ^ B[i];
        if (!diff) continue;
        // the lowest bit in 4 bit mismatch：
        // bit0 → base0 mismatch?
        // bit4 → base1 mismatch?
        // ...
        // bit60 → base15 mismatch?
        uint64_t t = diff | (diff >> 1) | (diff >> 2) | (diff >> 3);
        uint64_t nibbleMask = t & 0x1111111111111111ULL;
        mism += __builtin_popcountll(nibbleMask);
		if (mism > max_diff) return false;
    }
    return mism <= max_diff;
}


// -------------------------------
// AVX2 version
// -------------------------------
#if METAEDNA_X86 && (defined(__GNUC__) || defined(__clang__))

__attribute__((target("avx2")))
bool hamming_k_avx(const uint64_t* A, const uint64_t* B,
    size_t nword, size_t max_diff)
{
    size_t mism = 0;
    size_t i = 0;

    const __m256i nibmask = _mm256_set1_epi64x(0x1111111111111111ULL);

    for (; i + 4 <= nword; i += 4)
    {
        // load 4×uint64
        __m256i va = _mm256_loadu_si256((const __m256i*)(A + i));
        __m256i vb = _mm256_loadu_si256((const __m256i*)(B + i));
        __m256i diff = _mm256_xor_si256(va, vb);

        // t = diff | (diff>>1) | (diff>>2) | (diff>>3)
        __m256i t1 = _mm256_or_si256(diff, _mm256_srli_epi64(diff, 1));
        __m256i t2 = _mm256_or_si256(t1, _mm256_srli_epi64(diff, 2));
        __m256i t = _mm256_or_si256(t2, _mm256_srli_epi64(diff, 3));

        // nibbleMask = t & 0x1111111111111111ULL
        __m256i nib = _mm256_and_si256(t, nibmask);

        // --- popcount 全在寄存器里，不写回内存 ---

        // 将 nibbleMask 中的 4 个 64bit 部分独立 popcount
        // AVX2 没有 64-bit popcount，拆成两个 128-bit lane
        __m128i lo128 = _mm256_castsi256_si128(nib);
        __m128i hi128 = _mm256_extracti128_si256(nib, 1);

        // 128-bit → 2×64-bit
        uint64_t x0 = (uint64_t)_mm_cvtsi128_si64(lo128);
        uint64_t x1 = (uint64_t)_mm_cvtsi128_si64(_mm_srli_si128(lo128, 8));
        uint64_t x2 = (uint64_t)_mm_cvtsi128_si64(hi128);
        uint64_t x3 = (uint64_t)_mm_cvtsi128_si64(_mm_srli_si128(hi128, 8));

        mism += __builtin_popcountll(x0);
        if (mism > max_diff) return false;
        mism += __builtin_popcountll(x1);
        if (mism > max_diff) return false;
        mism += __builtin_popcountll(x2);
        if (mism > max_diff) return false;
        mism += __builtin_popcountll(x3);
        if (mism > max_diff) return false;
    }

    // scalar 部分保持一致性
    for (; i < nword; ++i) {
        uint64_t diff = A[i] ^ B[i];
        uint64_t t =
            diff |
            (diff >> 1) |
            (diff >> 2) |
            (diff >> 3);

        mism += __builtin_popcountll(t & 0x1111111111111111ULL);
        if (mism > max_diff) return false;
    }

    return mism <= max_diff;
}
#endif

// -------------------------------
// auto fallback
// -------------------------------
bool hamming_k(const uint64_t* A, const uint64_t* B,
    size_t nword, size_t max_diff, bool enable_avx2)
{

    if (CPU_SUPPORTS_AVX2() && enable_avx2) {
        return hamming_k_avx(A, B, nword, max_diff);
    }
    else {
        return hamming_k_scalar(A, B, nword, max_diff);
    }
}
