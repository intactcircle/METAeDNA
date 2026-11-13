
#include <vector>
#include <string>
#include <cstdint>
#include <array>

// 4-bit：A/C/G/T/N 以及简并碱基
static const std::array<uint8_t, 256> DNA4BIT = []() {
    std::array<uint8_t, 256> t{};
    t.fill(0xF); // 0xF unknown

    t['A'] = 0x1; t['C'] = 0x2; t['G'] = 0x3; t['T'] = 0x4;
    t['a'] = 0x1; t['c'] = 0x2; t['g'] = 0x3; t['t'] = 0x4;

    // 
    t['N'] = 0xF; t['n'] = 0xF;

    return t;
    }();


//  4bit-packed buffer（每 16 nt = 64 bit）
std::vector<uint64_t> encode4bit(const std::string& s);


// 4-bit packed buffer
bool hamming_k_scalar(const uint64_t* A, const uint64_t* B,
    size_t nword,  size_t max_diff);

// AVX2 
bool hamming_k_avx(const uint64_t* A, const uint64_t* B,
    size_t nword,  size_t max_diff);

// auto AVX2
bool hamming_k(const uint64_t* A, const uint64_t* B,
    size_t nword, size_t max_diff, bool enable_avx2);




