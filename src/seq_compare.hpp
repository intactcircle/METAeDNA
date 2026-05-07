
#include <vector>
#include <string>
#include <cstdint>
#include <array>

// 4-bit 编码表：A/C/G/T/N 以及简并碱基
static const std::array<uint8_t, 256> DNA4BIT = []() {
    std::array<uint8_t, 256> t{};
    t.fill(0xF); // 0xF 表示未知

    t['A'] = 0x1; t['C'] = 0x2; t['G'] = 0x3; t['T'] = 0x4;
    t['a'] = 0x1; t['c'] = 0x2; t['g'] = 0x3; t['t'] = 0x4;

    // 
    t['N'] = 0xF; t['n'] = 0xF;

    return t;
    }();


// 序列编码成 4bit-packed buffer（每 16 nt = 64 bit）
std::vector<uint64_t> encode4bit(const std::string& s);


// 标量版本：4-bit packed buffer
bool hamming_k_scalar(const uint64_t* A, const uint64_t* B,
    size_t nword,  size_t max_diff);

// AVX2 版本（如果支持）
bool hamming_k_avx(const uint64_t* A, const uint64_t* B,
    size_t nword,  size_t max_diff);

// 自动检测 AVX2
bool hamming_k(const uint64_t* A, const uint64_t* B,
    size_t nword, size_t max_diff, bool enable_avx2);




