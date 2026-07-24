#ifndef LIBRARY_FPS_FPS_HPP
#define LIBRARY_FPS_FPS_HPP

#line 1 "src/fps/fps.hpp"


#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>

#ifndef LIBRARY_MATH_CONVOLUTION_HPP
#define LIBRARY_MATH_CONVOLUTION_HPP

#line 1 "src/math/convolution.hpp"


#include <algorithm>
#include <cassert>
#include <limits>
#include <vector>

#ifndef LIBRARY_MATH_MODINT_HPP
#define LIBRARY_MATH_MODINT_HPP

#line 1 "src/math/modint.hpp"


#include <iostream>

template <int MOD>
struct static_modint {
    using mint = static_modint;

  private:
    int _v;

  public:
    static_modint() : _v(0) {}

    template <class T>
    static_modint(T v) {
        long long x = (long long)(v % MOD);
        if (x < 0) x += MOD;
        _v = int(x);
    }

    /// 0 以上 mod 未満の整数値を返す。
    int val() const { return _v; }
    /// 法 MOD を返す。
    static constexpr int mod() { return MOD; }

    /// n 乗を返す。n >= 0。
    mint pow(long long n) const {
        mint res(1), mul(*this);
        while (n > 0) {
            if (n & 1) res *= mul;
            mul *= mul;
            n >>= 1;
        }
        return res;
    }

    /// 逆元を返す。MOD は素数で、値は 0 でないこと。
    mint inv() const {
        return pow(MOD - 2);
    }

    mint& operator+=(const mint& a) {
        _v += a._v;
        if (_v >= MOD) _v -= MOD;
        return *this;
    }

    mint& operator-=(const mint& a) {
        _v -= a._v;
        if (_v < 0) _v += MOD;
        return *this;
    }

    mint& operator*=(const mint& a) {
        _v = int((long long)_v * a._v % MOD);
        return *this;
    }

    mint& operator/=(const mint& a) {
        return *this *= a.inv();
    }

    mint operator+() const { return *this; }
    mint operator-() const { return mint(0) - *this; }

    friend mint operator+(const mint& a, const mint& b) { return mint(a) += b; }
    friend mint operator-(const mint& a, const mint& b) { return mint(a) -= b; }
    friend mint operator*(const mint& a, const mint& b) { return mint(a) *= b; }
    friend mint operator/(const mint& a, const mint& b) { return mint(a) /= b; }
    friend bool operator==(const mint& a, const mint& b) { return a._v == b._v; }
    friend bool operator!=(const mint& a, const mint& b) { return a._v != b._v; }
    friend std::ostream& operator<<(std::ostream& os, const mint& a) { return os << a._v; }
    friend std::istream& operator>>(std::istream& is, mint& a) {
        long long v;
        is >> v;
        a = mint(v);
        return is;
    }
};

using modint998244353 = static_modint<998244353>;
using modint1000000007 = static_modint<1000000007>;

#endif  // LIBRARY_MATH_MODINT_HPP


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
                a[i + offset + p] = (l - r) * now;
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


namespace fps_internal {

template <class mint>
mint mod_sqrt(mint a) {
    if (a == mint(0)) return 0;
    if (mint::mod() == 2) return a;
    if (a.pow((mint::mod() - 1) / 2) != mint(1)) return 0;
    if (mint::mod() % 4 == 3) return a.pow((mint::mod() + 1) / 4);

    int s = __builtin_ctz(mint::mod() - 1);
    long long q = (mint::mod() - 1) >> s;
    mint z = 2;
    while (z.pow((mint::mod() - 1) / 2) != mint(-1)) z += 1;

    mint c = z.pow(q);
    mint x = a.pow((q + 1) / 2);
    mint t = a.pow(q);
    int m = s;
    while (t != mint(1)) {
        int i = 1;
        mint tt = t * t;
        while (tt != mint(1)) {
            tt *= tt;
            i++;
        }
        mint b = c.pow(1LL << (m - i - 1));
        x *= b;
        c = b * b;
        t *= c;
        m = i;
    }
    return x;
}

}  // namespace fps_internal

template <class mint>
struct formal_power_series : std::vector<mint> {
    using std::vector<mint>::vector;
    using fps = formal_power_series;

    formal_power_series(const std::vector<mint>& v) : std::vector<mint>(v) {}
    formal_power_series(std::vector<mint>&& v) : std::vector<mint>(std::move(v)) {}

    /// n 次未満、つまり先頭 n 項だけを返す。
    fps pre(int n) const {
        return fps(this->begin(), this->begin() + std::min(n, int(this->size())));
    }

    /// 末尾の 0 係数を削る。
    void shrink() {
        while (!this->empty() && this->back() == mint(0)) this->pop_back();
    }

    /// 多項式を x に代入した値を返す。
    mint eval(mint x) const {
        mint res = 0;
        for (int i = int(this->size()) - 1; i >= 0; i--) {
            res = res * x + (*this)[i];
        }
        return res;
    }

    fps& operator+=(const fps& rhs) {
        if (this->size() < rhs.size()) this->resize(rhs.size());
        for (int i = 0; i < int(rhs.size()); i++) (*this)[i] += rhs[i];
        return *this;
    }

    fps& operator-=(const fps& rhs) {
        if (this->size() < rhs.size()) this->resize(rhs.size());
        for (int i = 0; i < int(rhs.size()); i++) (*this)[i] -= rhs[i];
        return *this;
    }

    fps& operator*=(const fps& rhs) {
        auto res = convolution<mint>(*this, rhs);
        this->swap(res);
        return *this;
    }

    fps& operator*=(const mint& rhs) {
        for (auto& x : *this) x *= rhs;
        return *this;
    }

    fps& operator/=(const mint& rhs) {
        mint inv = rhs.inv();
        for (auto& x : *this) x *= inv;
        return *this;
    }

    friend fps operator+(fps lhs, const fps& rhs) { return lhs += rhs; }
    friend fps operator-(fps lhs, const fps& rhs) { return lhs -= rhs; }
    friend fps operator*(fps lhs, const fps& rhs) { return lhs *= rhs; }
    friend fps operator*(fps lhs, const mint& rhs) { return lhs *= rhs; }
    friend fps operator*(const mint& lhs, fps rhs) { return rhs *= lhs; }
    friend fps operator/(fps lhs, const mint& rhs) { return lhs /= rhs; }

    /// rhs と掛けて deg 項までに切る。
    fps mul_pre(const fps& rhs, int deg) const {
        if (deg <= 0) return {};
        fps lhs = pre(deg);
        fps rhs_pre = rhs.pre(deg);
        auto values = convolution<mint>(lhs, rhs_pre);
        fps res(std::move(values));
        if (int(res.size()) > deg) res.resize(deg);
        return res;
    }

    /// 形式微分を返す。
    fps diff() const {
        int n = int(this->size());
        if (n <= 1) return {};
        fps res(n - 1);
        for (int i = 1; i < n; i++) res[i - 1] = (*this)[i] * i;
        return res;
    }

    /// 定数項 0 の形式積分を返す。
    fps integral() const {
        int n = int(this->size());
        fps res(n + 1);
        assert(n < mint::mod());
        static std::vector<mint> invs{0, 1};
        if (int(invs.size()) <= n) {
            int old = int(invs.size());
            invs.resize(n + 1);
            for (int i = old; i <= n; i++) {
                invs[i] = -mint(mint::mod() / i) * invs[mint::mod() % i];
            }
        }
        for (int i = 0; i < n; i++) res[i + 1] = (*this)[i] * invs[i + 1];
        return res;
    }

    /// 逆数を deg 項まで返す。定数項が 0 でないこと。
    fps inv(int deg = -1) const {
        assert(!this->empty() && (*this)[0] != mint(0));
        if (deg == -1) deg = int(this->size());
        fps res{(*this)[0].inv()};
        for (int m = 1; m < deg; m <<= 1) {
            int nxt_deg = std::min(2 * m, deg);
            fps f = pre(nxt_deg);
            fps g = f.mul_pre(res, nxt_deg);
            g.resize(nxt_deg);
            g[0] = mint(2) - g[0];
            for (int i = 1; i < nxt_deg; i++) g[i] = -g[i];
            res = res.mul_pre(g, nxt_deg);
        }
        res.resize(deg);
        return res;
    }

    /// log(f) を deg 項まで返す。定数項が 1 であること。
    fps log(int deg = -1) const {
        assert(!this->empty() && (*this)[0] == mint(1));
        if (deg == -1) deg = int(this->size());
        if (deg == 0) return {};
        fps res = diff().mul_pre(inv(deg), deg - 1).integral();
        res.resize(deg);
        return res;
    }

    /// exp(f) を deg 項まで返す。定数項が 0 であること。
    fps exp(int deg = -1) const {
        assert(this->empty() || (*this)[0] == mint(0));
        if (deg == -1) deg = int(this->size());
        fps res{1};
        for (int m = 1; m < deg; m <<= 1) {
            int nxt_deg = std::min(2 * m, deg);
            fps f = pre(nxt_deg);
            fps g = res.log(nxt_deg);
            f.resize(nxt_deg);
            for (int i = 0; i < nxt_deg; i++) f[i] -= g[i];
            f[0] += 1;
            res = res.mul_pre(f, nxt_deg);
        }
        res.resize(deg);
        return res;
    }

    /// f^k を deg 項まで返す。k >= 0。
    fps pow(long long k, int deg = -1) const {
        if (deg == -1) deg = int(this->size());
        if (k == 0) {
            fps res(deg);
            if (deg > 0) res[0] = 1;
            return res;
        }

        int n = int(this->size());
        int lz = 0;
        while (lz < n && (*this)[lz] == mint(0)) lz++;
        if (lz == n || lz > (deg - 1) / k) return fps(deg);

        mint c = (*this)[lz];
        fps f(this->begin() + lz, this->end());
        f /= c;
        f = (f.log(deg - lz * k) * mint(k)).exp(deg - lz * k);
        f *= c.pow(k);

        fps res(deg);
        for (int i = 0; i < int(f.size()) && i + lz * k < deg; i++) {
            res[i + lz * k] = f[i];
        }
        return res;
    }

    /// 平方根を deg 項まで返す。存在しない場合は空を返す。
    fps sqrt(int deg = -1) const {
        if (deg == -1) deg = int(this->size());
        if (this->empty()) return fps(deg);
        if ((*this)[0] == mint(0)) {
            for (int i = 1; i < int(this->size()); i++) {
                if ((*this)[i] != mint(0)) {
                    if (i & 1) return {};
                    if (i / 2 >= deg) return fps(deg);
                    fps f(this->begin() + i, this->end());
                    fps res = f.sqrt(deg - i / 2);
                    if (res.empty() && deg - i / 2 > 0) return {};
                    fps shifted(deg);
                    for (int j = 0; j < int(res.size()) && j + i / 2 < deg; j++) {
                        shifted[j + i / 2] = res[j];
                    }
                    return shifted;
                }
            }
            return fps(deg);
        }

        mint s = fps_internal::mod_sqrt((*this)[0]);
        if (s * s != (*this)[0]) return {};

        fps res{s};
        mint inv2 = mint(2).inv();
        for (int m = 1; m < deg; m <<= 1) {
            int nxt_deg = std::min(2 * m, deg);
            fps f = pre(nxt_deg);
            res = (res + f.mul_pre(res.inv(nxt_deg), nxt_deg)) * inv2;
            res.resize(nxt_deg);
        }
        res.resize(deg);
        return res;
    }
};

template <class mint>
using fps = formal_power_series<mint>;

/// Bostan-Mori 法。P(x) / Q(x) の x^n 係数を O(d log d log n) で返す。Q[0] != 0。
template <class mint>
mint bostan_mori(formal_power_series<mint> p, formal_power_series<mint> q, long long n) {
    assert(!q.empty() && q[0] != mint(0));
    p.shrink();
    q.shrink();
    while (n > 0) {
        formal_power_series<mint> q_neg = q;
        for (int i = 1; i < int(q_neg.size()); i += 2) q_neg[i] = -q_neg[i];

        formal_power_series<mint> s = p * q_neg;
        formal_power_series<mint> t = q * q_neg;

        int parity = int(n & 1);
        int np_size = int(s.size()) <= parity ? 0 : (int(s.size()) - parity + 1) / 2;
        formal_power_series<mint> np(np_size);
        for (int i = 0; i < np_size; i++) np[i] = s[2 * i + parity];

        int nq_size = (int(t.size()) + 1) / 2;
        formal_power_series<mint> nq(nq_size);
        for (int i = 0; i < nq_size; i++) nq[i] = t[2 * i];

        p = std::move(np);
        q = std::move(nq);
        n >>= 1;
    }
    return p.empty() ? mint(0) : p[0] / q[0];
}

/// Berlekamp-Massey。s[n] = sum_{i=0}^{d-1} c[i] * s[n-i-1] となる最小漸化式 c を返す。
template <class mint>
std::vector<mint> berlekamp_massey(const std::vector<mint>& s) {
    std::vector<mint> c{1}, b{1};
    c.reserve(s.size() + 1);
    b.reserve(s.size() + 1);
    int l = 0, m = 1;
    mint bb = 1;

    for (int n = 0; n < int(s.size()); n++) {
        mint d = s[n];
        for (int i = 1; i <= l; i++) d += c[i] * s[n - i];
        if (d == mint(0)) {
            m++;
            continue;
        }

        std::vector<mint> t = c;
        mint coef = d / bb;
        if (int(c.size()) < int(b.size()) + m) c.resize(b.size() + m);
        for (int i = 0; i < int(b.size()); i++) c[i + m] -= coef * b[i];

        if (2 * l <= n) {
            l = n + 1 - l;
            b = t;
            bb = d;
            m = 1;
        } else {
            m++;
        }
    }

    std::vector<mint> res(l);
    for (int i = 0; i < l; i++) res[i] = -c[i + 1];
    return res;
}

/// 線形漸化式の k 項目を返す。coef は s[n] = sum coef[i] * s[n-i-1] の形。
template <class mint>
mint linear_recurrence_kth(const std::vector<mint>& init, const std::vector<mint>& coef, long long k) {
    int d = int(coef.size());
    assert(d >= 0);
    if (k < int(init.size())) return init[k];
    if (d == 0) return mint(0);
    assert(int(init.size()) >= d);

    formal_power_series<mint> a(init.begin(), init.begin() + d);
    formal_power_series<mint> q(d + 1);
    q[0] = 1;
    for (int i = 0; i < d; i++) q[i + 1] = -coef[i];

    formal_power_series<mint> p = (a * q).pre(d);
    return bostan_mori(p, q, k);
}

/// 先頭項列 s から BM で漸化式を推定し、Bostan-Mori で k 項目を返す。
template <class mint>
mint bmbm(const std::vector<mint>& s, long long k) {
    auto coef = berlekamp_massey(s);
    return linear_recurrence_kth(s, coef, k);
}

#endif  // LIBRARY_FPS_FPS_HPP
