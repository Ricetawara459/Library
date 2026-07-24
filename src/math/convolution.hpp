#ifndef LIBRARY_MATH_CONVOLUTION_HPP
#define LIBRARY_MATH_CONVOLUTION_HPP

#pragma once

#include <algorithm>
#include <cassert>
#include <limits>
#include <vector>

#include "modint.hpp"

namespace convolution_internal {

constexpr int naive_threshold = 60;

template <class mint>
std::vector<mint> convolution_naive(const std::vector<mint>& a,
                                    const std::vector<mint>& b) {
    if (a.empty() || b.empty()) return {};
    std::vector<mint> res(a.size() + b.size() - 1);
    if (a.size() < b.size()) {
        for (int i = 0; i < int(a.size()); i++) {
            for (int j = 0; j < int(b.size()); j++) res[i + j] += a[i] * b[j];
        }
    } else {
        for (int j = 0; j < int(b.size()); j++) {
            for (int i = 0; i < int(a.size()); i++) res[i + j] += a[i] * b[j];
        }
    }
    return res;
}

template <class mint, int primitive_root>
void butterfly(std::vector<mint>& a) {
    int n = int(a.size());
    int h = 0;
    while ((1 << h) < n) h++;

    static const std::vector<mint> sum_e = [] {
        int cnt = __builtin_ctz(mint::mod() - 1);
        std::vector<mint> es(cnt - 1), ies(cnt - 1), sum(cnt - 1);
        mint e = mint(primitive_root).pow((mint::mod() - 1) >> cnt);
        mint ie = e.inv();
        for (int i = cnt; i >= 2; i--) {
            es[i - 2] = e;
            ies[i - 2] = ie;
            e *= e;
            ie *= ie;
        }
        mint now = 1;
        for (int i = 0; i <= cnt - 2; i++) {
            sum[i] = es[i] * now;
            now *= ies[i];
        }
        return sum;
    }();

    for (int ph = 1; ph <= h; ph++) {
        int w = 1 << (ph - 1);
        int p = 1 << (h - ph);
        mint now = 1;
        for (int s = 0; s < w; s++) {
            int offset = s << (h - ph + 1);
            for (int i = 0; i < p; i++) {
                mint l = a[i + offset];
                mint r = a[i + offset + p] * now;
                a[i + offset] = l + r;
                a[i + offset + p] = l - r;
            }
            if (s + 1 < w) now *= sum_e[__builtin_ctz(~unsigned(s))];
        }
    }
}

template <class mint, int primitive_root>
void butterfly_inv(std::vector<mint>& a) {
    int n = int(a.size());
    int h = 0;
    while ((1 << h) < n) h++;

    static const std::vector<mint> sum_ie = [] {
        int cnt = __builtin_ctz(mint::mod() - 1);
        std::vector<mint> es(cnt - 1), ies(cnt - 1), sum(cnt - 1);
        mint e = mint(primitive_root).pow((mint::mod() - 1) >> cnt);
        mint ie = e.inv();
        for (int i = cnt; i >= 2; i--) {
            es[i - 2] = e;
            ies[i - 2] = ie;
            e *= e;
            ie *= ie;
        }
        mint now = 1;
        for (int i = 0; i <= cnt - 2; i++) {
            sum[i] = ies[i] * now;
            now *= es[i];
        }
        return sum;
    }();

    for (int ph = h; ph >= 1; ph--) {
        int w = 1 << (ph - 1);
        int p = 1 << (h - ph);
        mint now = 1;
        for (int s = 0; s < w; s++) {
            int offset = s << (h - ph + 1);
            for (int i = 0; i < p; i++) {
                mint l = a[i + offset];
                mint r = a[i + offset + p];
                a[i + offset] = l + r;
                unsigned long long diff =
                    (unsigned long long)(mint::mod()) + l.val() - r.val();
                a[i + offset + p] =
                    mint::raw(int(diff * now.val() % mint::mod()));
            }
            if (s + 1 < w) now *= sum_ie[__builtin_ctz(~unsigned(s))];
        }
    }
}

template <class mint, int primitive_root>
std::vector<mint> convolution_ntt(std::vector<mint> a, std::vector<mint> b) {
    int n = int(a.size());
    int m = int(b.size());
    if (n == 0 || m == 0) return {};

    int z = 1;
    while (z < n + m - 1) z <<= 1;
    assert((mint::mod() - 1) % z == 0);
    a.resize(z);
    b.resize(z);

    butterfly<mint, primitive_root>(a);
    butterfly<mint, primitive_root>(b);
    for (int i = 0; i < z; i++) a[i] *= b[i];
    butterfly_inv<mint, primitive_root>(a);

    mint iz = mint(z).inv();
    a.resize(n + m - 1);
    for (auto& x : a) x *= iz;
    return a;
}

constexpr long long pow_mod_general(long long x, long long n, long long mod) {
    long long res = 1;
    x %= mod;
    while (n > 0) {
        if (n & 1) res = res * x % mod;
        x = x * x % mod;
        n >>= 1;
    }
    return res;
}

constexpr long long inv_mod_general(long long x, long long mod) {
    return pow_mod_general(x, mod - 2, mod);
}

inline long long mod_norm(long long x, long long mod) {
    x %= mod;
    if (x < 0) x += mod;
    return x;
}

template <int mod, int primitive_root>
std::vector<long long> convolution_prime(const std::vector<long long>& a,
                                         const std::vector<long long>& b) {
    using mint = static_modint<mod>;
    std::vector<mint> am(a.size()), bm(b.size());
    for (int i = 0; i < int(a.size()); i++) am[i] = mint(a[i]);
    for (int i = 0; i < int(b.size()); i++) bm[i] = mint(b[i]);

    std::vector<mint> cm;
    if (std::min(am.size(), bm.size()) <= naive_threshold) {
        cm = convolution_naive(am, bm);
    } else {
        cm = convolution_ntt<mint, primitive_root>(std::move(am), std::move(bm));
    }

    std::vector<long long> res(cm.size());
    for (int i = 0; i < int(cm.size()); i++) res[i] = cm[i].val();
    return res;
}

}  // namespace convolution_internal

/// 任意 modint 用の畳み込み。短い場合は愚直法、それ以外は3つのNTTとGarnerで復元する。
template <class mint>
std::vector<mint> any_mod_convolution(const std::vector<mint>& a,
                                      const std::vector<mint>& b) {
    int n = int(a.size());
    int m = int(b.size());
    if (n == 0 || m == 0) return {};
    if (std::min(n, m) <= convolution_internal::naive_threshold) {
        return convolution_internal::convolution_naive(a, b);
    }

    std::vector<long long> al(n), bl(m);
    for (int i = 0; i < n; i++) al[i] = a[i].val();
    for (int i = 0; i < m; i++) bl[i] = b[i].val();

    constexpr long long m1 = 998244353;
    constexpr long long m2 = 754974721;
    constexpr long long m3 = 469762049;
    auto c1 = convolution_internal::convolution_prime<m1, 3>(al, bl);
    auto c2 = convolution_internal::convolution_prime<m2, 11>(al, bl);
    auto c3 = convolution_internal::convolution_prime<m3, 3>(al, bl);

    constexpr long long m1_inv_m2 = convolution_internal::inv_mod_general(m1, m2);
    constexpr long long m1_inv_m3 = convolution_internal::inv_mod_general(m1, m3);
    constexpr long long m2_inv_m3 = convolution_internal::inv_mod_general(m2, m3);
    const mint mint_m1 = mint(m1);
    const mint mint_m1m2 = mint_m1 * mint(m2);

    std::vector<mint> res(n + m - 1);
    for (int i = 0; i < n + m - 1; i++) {
        long long v1 = c1[i];
        long long v2 = convolution_internal::mod_norm(
            (c2[i] - v1) * m1_inv_m2, m2);
        long long t3 = convolution_internal::mod_norm(
            (c3[i] - v1) * m1_inv_m3, m3);
        t3 = convolution_internal::mod_norm(t3 - v2, m3);
        long long v3 = t3 * m2_inv_m3 % m3;
        res[i] = mint(v1) + mint(v2) * mint_m1 + mint(v3) * mint_m1m2;
    }
    return res;
}

/// modint 配列の畳み込みを返す。短い場合は愚直法を使う。
template <class mint>
std::vector<mint> convolution(const std::vector<mint>& a,
                              const std::vector<mint>& b) {
    if (a.empty() || b.empty()) return {};
    if (std::min(a.size(), b.size()) <= convolution_internal::naive_threshold) {
        return convolution_internal::convolution_naive(a, b);
    }
    if constexpr (mint::mod() == 998244353) {
        return convolution_internal::convolution_ntt<mint, 3>(a, b);
    } else {
        return any_mod_convolution(a, b);
    }
}

/// 所有権を渡せる配列の畳み込み。入力領域を再利用してコピーを減らす。
template <class mint>
std::vector<mint> convolution(std::vector<mint>&& a,
                              std::vector<mint>&& b) {
    if (a.empty() || b.empty()) return {};
    if (std::min(a.size(), b.size()) <= convolution_internal::naive_threshold) {
        return convolution_internal::convolution_naive(a, b);
    }
    if constexpr (mint::mod() == 998244353) {
        return convolution_internal::convolution_ntt<mint, 3>(
            std::move(a), std::move(b));
    } else {
        return any_mod_convolution(a, b);
    }
}

/// 998244353 上で long long 配列の畳み込みを返す。
inline std::vector<long long> convolution(const std::vector<long long>& a,
                                          const std::vector<long long>& b) {
    return convolution_internal::convolution_prime<998244353, 3>(a, b);
}

/// long long の範囲で正確に畳み込みを返す。
inline std::vector<long long> convolution_ll(const std::vector<long long>& a,
                                             const std::vector<long long>& b) {
    int n = int(a.size());
    int m = int(b.size());
    if (n == 0 || m == 0) return {};
    if (std::min(n, m) <= convolution_internal::naive_threshold) {
        std::vector<__int128> wide(n + m - 1);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) wide[i + j] += (__int128)a[i] * b[j];
        }
        std::vector<long long> res(n + m - 1);
        for (int i = 0; i < n + m - 1; i++) {
            assert(wide[i] >= std::numeric_limits<long long>::min());
            assert(wide[i] <= std::numeric_limits<long long>::max());
            res[i] = (long long)wide[i];
        }
        return res;
    }

    constexpr long long m1 = 998244353;
    constexpr long long m2 = 754974721;
    constexpr long long m3 = 469762049;
    auto c1 = convolution_internal::convolution_prime<m1, 3>(a, b);
    auto c2 = convolution_internal::convolution_prime<m2, 11>(a, b);
    auto c3 = convolution_internal::convolution_prime<m3, 3>(a, b);

    constexpr long long m1_inv_m2 = convolution_internal::inv_mod_general(m1, m2);
    constexpr long long m1_inv_m3 = convolution_internal::inv_mod_general(m1, m3);
    constexpr long long m2_inv_m3 = convolution_internal::inv_mod_general(m2, m3);

    std::vector<long long> res(n + m - 1);
    for (int i = 0; i < n + m - 1; i++) {
        long long v1 = c1[i];
        long long v2 = convolution_internal::mod_norm(
            (c2[i] - v1) * m1_inv_m2, m2);
        long long t3 = convolution_internal::mod_norm(
            (c3[i] - v1) * m1_inv_m3, m3);
        t3 = convolution_internal::mod_norm(t3 - v2, m3);
        long long v3 = t3 * m2_inv_m3 % m3;
        __int128 value = (__int128)v1 + (__int128)v2 * m1 +
                         (__int128)v3 * m1 * m2;
        constexpr __int128 mod_product = (__int128)m1 * m2 * m3;
        if (value > mod_product / 2) value -= mod_product;
        assert(value >= std::numeric_limits<long long>::min());
        assert(value <= std::numeric_limits<long long>::max());
        res[i] = (long long)value;
    }
    return res;
}

#endif  // LIBRARY_MATH_CONVOLUTION_HPP
