#pragma once

#include <vector>
#include <algorithm>
#include <cassert>

namespace convolution_internal {

// 汎用冪乗算
constexpr long long pow_mod_general(long long x, long long n, long long m) {
    long long r = 1;
    long long base = x % m;
    while (n > 0) {
        if (n & 1) r = (long long)((__int128)r * base % m);
        base = (long long)((__int128)base * base % m);
        n >>= 1;
    }
    return r;
}

// 汎用逆元計算
constexpr long long inv_mod_general(long long x, long long m) {
    return pow_mod_general(x, m - 2, m);
}

// 任意のMOD, 原始根Gに対応するインプレース数論変換 (DIF NTT)
template <long long MOD, long long G>
void butterfly_ext(std::vector<long long>& a) {
    int n = int(a.size());
    int h = 0;
    while ((1 << h) < n) h++;
    
    static bool first = true;
    static long long sum_e[30];
    if (first) {
        first = false;
        long long es[30], ies[30];
        int cnt = 23; 
        long long e = pow_mod_general(G, (MOD - 1) >> cnt, MOD);
        long long ie = inv_mod_general(e, MOD);
        for (int i = cnt; i >= 2; i--) {
            es[i - 2] = e;
            ies[i - 2] = ie;
            e = (long long)((__int128)e * e % MOD);
            ie = (long long)((__int128)ie * ie % MOD);
        }
        long long now = 1;
        for (int i = 0; i <= cnt - 2; i++) {
            sum_e[i] = (long long)((__int128)es[i] * now % MOD);
            now = (long long)((__int128)now * ies[i] % MOD);
        }
    }
    
    for (int ph = 1; ph <= h; ph++) {
        int w = 1 << (ph - 1), p = 1 << (h - ph);
        long long now = 1;
        for (int s = 0; s < w; s++) {
            int offset = s << (h - ph + 1);
            for (int i = 0; i < p; i++) {
                auto l = a[i + offset];
                auto r = (long long)((__int128)a[i + offset + p] * now % MOD);
                a[i + offset] = (l + r >= MOD ? l + r - MOD : l + r);
                a[i + offset + p] = (l - r < 0 ? l - r + MOD : l - r);
            }
            now = (long long)((__int128)now * sum_e[__builtin_ctz(~(unsigned int)(s))] % MOD);
        }
    }
}

// 任意のMOD, 原始根Gに対応するインプレース逆数論変換 (DIT INTT)
template <long long MOD, long long G>
void butterfly_inv_ext(std::vector<long long>& a) {
    int n = int(a.size());
    int h = 0;
    while ((1 << h) < n) h++;
    
    static bool first = true;
    static long long sum_ie[30];
    if (first) {
        first = false;
        long long es[30], ies[30];
        int cnt = 23; 
        long long e = pow_mod_general(G, (MOD - 1) >> cnt, MOD);
        long long ie = inv_mod_general(e, MOD);
        for (int i = cnt; i >= 2; i--) {
            es[i - 2] = e;
            ies[i - 2] = ie;
            e = (long long)((__int128)e * e % MOD);
            ie = (long long)((__int128)ie * ie % MOD);
        }
        long long now = 1;
        for (int i = 0; i <= cnt - 2; i++) {
            sum_ie[i] = (long long)((__int128)ies[i] * now % MOD);
            now = (long long)((__int128)now * es[i] % MOD);
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
                a[i + offset] = (l + r >= MOD ? l + r - MOD : l + r);
                a[i + offset + p] = (l - r < 0 ? l - r + MOD : l - r);
                a[i + offset + p] = (long long)((__int128)a[i + offset + p] * now % MOD);
            }
            now = (long long)((__int128)now * sum_ie[__builtin_ctz(~(unsigned int)(s))] % MOD);
        }
    }
}

// 任意MOD用内部畳み込み関数
template <long long MOD, long long G>
std::vector<long long> convolution_ext(std::vector<long long> a, std::vector<long long> b) {
    int n = int(a.size()), m = int(b.size());
    if (n == 0 || m == 0) return {};
    int z = 1;
    while (z < n + m - 1) z <<= 1;
    
    a.resize(z, 0);
    b.resize(z, 0);
    
    butterfly_ext<MOD, G>(a);
    butterfly_ext<MOD, G>(b);
    for (int i = 0; i < z; i++) {
        a[i] = (long long)((__int128)a[i] * b[i] % MOD);
    }
    butterfly_inv_ext<MOD, G>(a);
    
    a.resize(n + m - 1);
    long long iz = inv_mod_general(z, MOD);
    for (int i = 0; i < n + m - 1; i++) {
        a[i] = (long long)((__int128)a[i] * iz % MOD);
    }
    return a;
}

} // namespace convolution_internal

/// 任意 modint 用の畳み込み。Garner で復元する。
template <class mint>
std::vector<mint> any_mod_convolution(const std::vector<mint>& a, const std::vector<mint>& b) {
    int n = int(a.size()), m = int(b.size());
    if (n == 0 || m == 0) return {};
    
    std::vector<long long> a_ll(n), b_ll(m);
    for (int i = 0; i < n; i++) a_ll[i] = a[i].val();
    for (int i = 0; i < m; i++) b_ll[i] = b[i].val();
    
    constexpr long long m1 = 998244353;
    constexpr long long m2 = 754974721;
    constexpr long long m3 = 469762049;
    
    auto c1 = convolution_internal::convolution_ext<m1, 3>(a_ll, b_ll);
    auto c2 = convolution_internal::convolution_ext<m2, 11>(a_ll, b_ll);
    auto c3 = convolution_internal::convolution_ext<m3, 3>(a_ll, b_ll);
    
    int sz = n + m - 1;
    std::vector<mint> res(sz);
    
    constexpr long long m1_inv_m2 = convolution_internal::inv_mod_general(m1, m2);
    constexpr long long m1_inv_m3 = convolution_internal::inv_mod_general(m1, m3);
    constexpr long long m2_inv_m3 = convolution_internal::inv_mod_general(m2, m3);
    
    mint mint_m1 = mint(m1);
    mint mint_m1m2 = mint_m1 * mint(m2);
    
    for (int i = 0; i < sz; i++) {
        long long v1 = c1[i];
        long long v2 = (long long)((__int128)(c2[i] - v1 + m2) * m1_inv_m2 % m2);
        long long v3 = (long long)((__int128)((c3[i] - v1 + m3) * m1_inv_m3 % m3 - v2 + m3) * m2_inv_m3 % m3);
        
        res[i] = mint(v1) + mint(v2) * mint_m1 + mint(v3) * mint_m1m2;
    }
    return res;
}

/// modint 配列の畳み込みを返す。998244353 なら NTT、それ以外は Garner。
template <class mint>
std::vector<mint> convolution(const std::vector<mint>& a, const std::vector<mint>& b) {
    if constexpr (mint::mod() == 998244353) {
        int n = int(a.size()), m = int(b.size());
        if (n == 0 || m == 0) return {};
        std::vector<long long> a_ll(n), b_ll(m);
        for (int i = 0; i < n; i++) a_ll[i] = a[i].val();
        for (int i = 0; i < m; i++) b_ll[i] = b[i].val();
        auto c_ll = convolution_internal::convolution_ext<998244353, 3>(a_ll, b_ll);
        std::vector<mint> res(n + m - 1);
        for (int i = 0; i < n + m - 1; i++) res[i] = mint(c_ll[i]);
        return res;
    } else {
        // MOD が 998244353 以外なら自動で Garner へフォールバック
        return any_mod_convolution<mint>(a, b);
    }
}


/// 998244353 上で long long 配列の畳み込みを返す。
std::vector<long long> convolution(const std::vector<long long>& a, const std::vector<long long>& b) {
    return convolution_internal::convolution_ext<998244353, 3>(a, b);
}


/// long long の範囲で正確に畳み込みを返す。
std::vector<long long> convolution_ll(const std::vector<long long>& a, const std::vector<long long>& b) {
    int n = int(a.size()), m = int(b.size());
    if (n == 0 || m == 0) return {};
    
    constexpr long long m1 = 998244353;
    constexpr long long m2 = 754974721;
    constexpr long long m3 = 469762049;
    
    auto c1 = convolution_internal::convolution_ext<m1, 3>(a, b);
    auto c2 = convolution_internal::convolution_ext<m2, 11>(a, b);
    auto c3 = convolution_internal::convolution_ext<m3, 3>(a, b);
    
    int sz = n + m - 1;
    std::vector<long long> res(sz);
    
    constexpr long long m1_inv_m2 = convolution_internal::inv_mod_general(m1, m2);
    constexpr long long m1_inv_m3 = convolution_internal::inv_mod_general(m1, m3);
    constexpr long long m2_inv_m3 = convolution_internal::inv_mod_general(m2, m3);
    
    for (int i = 0; i < sz; i++) {
        long long v1 = c1[i];
        long long v2 = (long long)((__int128)(c2[i] - v1 + m2) * m1_inv_m2 % m2);
        long long v3 = (long long)((__int128)((c3[i] - v1 + m3) * m1_inv_m3 % m3 - v2 + m3) * m2_inv_m3 % m3);
        
        res[i] = (long long)((__int128)v1 + (__int128)v2 * m1 + (__int128)v3 * m1 * m2);
    }
    return res;
}
