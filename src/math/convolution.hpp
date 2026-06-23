#pragma once

#include <vector>
#include <algorithm>
#include <cassert>

namespace internal {

constexpr long long CHAR_MOD = 998244353;

// 高速冪乗算 (内部用)
constexpr long long pow_mod_const(long long x, long long n) {
    long long r = 1;
    long long base = x % CHAR_MOD;
    while (n > 0) {
        if (n & 1) r = (r * base) % CHAR_MOD;
        base = (base * base) % CHAR_MOD;
        n >>= 1;
    }
    return r;
}

// 逆元計算 (内部用)
constexpr long long inv_mod_const(long long x) {
    return pow_mod_const(x, CHAR_MOD - 2);
}

// インプレース数論変換 (DIF NTT)
void butterfly(std::vector<long long>& a) {
    int n = int(a.size());
    int h = 0;
    while ((1 << h) < n) h++;
    
    static bool first = true;
    static long long sum_e[30];
    if (first) {
        first = false;
        long long es[30], ies[30];
        long long e = 911660635; // g^((MOD-1)/2^23)
        long long ie = 86583718; // e^-1
        for (int i = 23; i >= 2; i--) {
            es[i - 2] = e;
            ies[i - 2] = ie;
            e = (e * e) % CHAR_MOD;
            ie = (ie * ie) % CHAR_MOD;
        }
        long long now = 1;
        for (int i = 0; i <= 21; i++) {
            sum_e[i] = (es[i] * now) % CHAR_MOD;
            now = (now * ies[i]) % CHAR_MOD;
        }
    }
    
    for (int ph = 1; ph <= h; ph++) {
        int w = 1 << (ph - 1), p = 1 << (h - ph);
        long long now = 1;
        for (int s = 0; s < w; s++) {
            int offset = s << (h - ph + 1);
            for (int i = 0; i < p; i++) {
                auto l = a[i + offset];
                auto r = (a[i + offset + p] * now) % CHAR_MOD;
                a[i + offset] = (l + r >= CHAR_MOD ? l + r - CHAR_MOD : l + r);
                a[i + offset + p] = (l - r < 0 ? l - r + CHAR_MOD : l - r);
            }
            now = (now * sum_e[__builtin_ctz(~(unsigned int)(s))]) % CHAR_MOD;
        }
    }
}

// インプレース逆数論変換 (DIT INTT)
void butterfly_inv(std::vector<long long>& a) {
    int n = int(a.size());
    int h = 0;
    while ((1 << h) < n) h++;
    
    static bool first = true;
    static long long sum_ie[30];
    if (first) {
        first = false;
        long long es[30], ies[30];
        long long e = 911660635;
        long long ie = 86583718;
        for (int i = 23; i >= 2; i--) {
            es[i - 2] = e;
            ies[i - 2] = ie;
            e = (e * e) % CHAR_MOD;
            ie = (ie * ie) % CHAR_MOD;
        }
        long long now = 1;
        for (int i = 0; i <= 21; i++) {
            sum_ie[i] = (ies[i] * now) % CHAR_MOD;
            now = (now * es[i]) % CHAR_MOD;
        }
    }
    
    for (int ph = h; ph >= 1; ph--) {
        int w = 1 << (ph - 1), p = 1 << (h - ph);
        long long now = 1;
        for (int s = 0; s < w; s++) {
            int offset = s << (h - ph + 1);
            for (int i = 0; i < p; i++) {
                auto l = a[i + offset];
                auto r = a[i + offset + p];
                a[i + offset] = (l + r >= CHAR_MOD ? l + r - CHAR_MOD : l + r);
                a[i + offset + p] = (l - r < 0 ? l - r + CHAR_MOD : l - r);
                a[i + offset + p] = (a[i + offset + p] * now) % CHAR_MOD;
            }
            now = (now * sum_ie[__builtin_ctz(~(unsigned int)(s))]) % CHAR_MOD;
        }
    }
}

} // namespace internal

// 畳み込み関数本体
std::vector<long long> convolution(std::vector<long long> a, std::vector<long long> b) {
    int n = int(a.size()), m = int(b.size());
    if (n == 0 || m == 0) return {};
    
    int z = 1;
    while (z < n + m - 1) z <<= 1;
    
    a.resize(z, 0);
    b.resize(z, 0);
    
    internal::butterfly(a);
    internal::butterfly(b);
    for (int i = 0; i < z; i++) {
        a[i] = (a[i] * b[i]) % internal::CHAR_MOD;
    }
    internal::butterfly_inv(a);
    
    a.resize(n + m - 1);
    long long iz = internal::inv_mod_const(z);
    for (int i = 0; i < n + m - 1; i++) {
        a[i] = (a[i] * iz) % internal::CHAR_MOD;
    }
    return a;
}