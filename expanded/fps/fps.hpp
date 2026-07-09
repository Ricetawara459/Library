#line 2 "src/fps/fps.hpp"

#include <algorithm>
#include <cassert>
#include <vector>

#line 2 "src/math/convolution.hpp"

#line 6 "src/math/convolution.hpp"

namespace internal {

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

} // namespace internal

// 1. 任意の static_modint<MOD> と明示的に連携する畳み込み (Garner)
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
    
    auto c1 = internal::convolution_ext<m1, 3>(a_ll, b_ll);
    auto c2 = internal::convolution_ext<m2, 11>(a_ll, b_ll);
    auto c3 = internal::convolution_ext<m3, 3>(a_ll, b_ll);
    
    int sz = n + m - 1;
    std::vector<mint> res(sz);
    
    constexpr long long m1_inv_m2 = internal::inv_mod_general(m1, m2);
    constexpr long long m1_inv_m3 = internal::inv_mod_general(m1, m3);
    constexpr long long m2_inv_m3 = internal::inv_mod_general(m2, m3);
    
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

// 2. 標準の mint 版 (★ if constexpr により MOD に応じて自動で最適な処理に分岐)
template <class mint>
std::vector<mint> convolution(const std::vector<mint>& a, const std::vector<mint>& b) {
    if constexpr (mint::mod() == 998244353) {
        int n = int(a.size()), m = int(b.size());
        if (n == 0 || m == 0) return {};
        std::vector<long long> a_ll(n), b_ll(m);
        for (int i = 0; i < n; i++) a_ll[i] = a[i].val();
        for (int i = 0; i < m; i++) b_ll[i] = b[i].val();
        auto c_ll = internal::convolution_ext<998244353, 3>(a_ll, b_ll);
        std::vector<mint> res(n + m - 1);
        for (int i = 0; i < n + m - 1; i++) res[i] = mint(c_ll[i]);
        return res;
    } else {
        // MOD が 998244353 以外なら自動で Garner へフォールバック
        return any_mod_convolution<mint>(a, b);
    }
}

// 3. 標準の long long 版 (Mod 998244353 固定)
std::vector<long long> convolution(const std::vector<long long>& a, const std::vector<long long>& b) {
    return internal::convolution_ext<998244353, 3>(a, b);
}

// 4. 拡張：long long の範囲で正確に計算する畳み込み
std::vector<long long> convolution_ll(const std::vector<long long>& a, const std::vector<long long>& b) {
    int n = int(a.size()), m = int(b.size());
    if (n == 0 || m == 0) return {};
    
    constexpr long long m1 = 998244353;
    constexpr long long m2 = 754974721;
    constexpr long long m3 = 469762049;
    
    auto c1 = internal::convolution_ext<m1, 3>(a, b);
    auto c2 = internal::convolution_ext<m2, 11>(a, b);
    auto c3 = internal::convolution_ext<m3, 3>(a, b);
    
    int sz = n + m - 1;
    std::vector<long long> res(sz);
    
    constexpr long long m1_inv_m2 = internal::inv_mod_general(m1, m2);
    constexpr long long m1_inv_m3 = internal::inv_mod_general(m1, m3);
    constexpr long long m2_inv_m3 = internal::inv_mod_general(m2, m3);
    
    for (int i = 0; i < sz; i++) {
        long long v1 = c1[i];
        long long v2 = (long long)((__int128)(c2[i] - v1 + m2) * m1_inv_m2 % m2);
        long long v3 = (long long)((__int128)((c3[i] - v1 + m3) * m1_inv_m3 % m3 - v2 + m3) * m2_inv_m3 % m3);
        
        res[i] = (long long)((__int128)v1 + (__int128)v2 * m1 + (__int128)v3 * m1 * m2);
    }
    return res;
}
#line 8 "src/fps/fps.hpp"

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

    fps pre(int n) const {
        return fps(this->begin(), this->begin() + std::min(n, int(this->size())));
    }

    void shrink() {
        while (!this->empty() && this->back() == mint(0)) this->pop_back();
    }

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
        this->assign(res.begin(), res.end());
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

    fps diff() const {
        int n = int(this->size());
        if (n <= 1) return {};
        fps res(n - 1);
        for (int i = 1; i < n; i++) res[i - 1] = (*this)[i] * i;
        return res;
    }

    fps integral() const {
        int n = int(this->size());
        fps res(n + 1);
        static std::vector<mint> invs{0, 1};
        if (int(invs.size()) <= n) {
            int old = int(invs.size());
            invs.resize(n + 1);
            for (int i = old; i <= n; i++) invs[i] = mint(i).inv();
        }
        for (int i = 0; i < n; i++) res[i + 1] = (*this)[i] * invs[i + 1];
        return res;
    }

    fps inv(int deg = -1) const {
        assert(!this->empty() && (*this)[0] != mint(0));
        if (deg == -1) deg = int(this->size());
        fps res{(*this)[0].inv()};
        for (int m = 1; m < deg; m <<= 1) {
            fps f = pre(2 * m);
            fps nxt = res * res * f;
            res *= mint(2);
            res -= nxt;
            res.resize(std::min(2 * m, deg));
        }
        return res.pre(deg);
    }

    fps log(int deg = -1) const {
        assert(!this->empty() && (*this)[0] == mint(1));
        if (deg == -1) deg = int(this->size());
        if (deg == 0) return {};
        fps res = (diff() * inv(deg)).pre(deg - 1).integral();
        res.resize(deg);
        return res;
    }

    fps exp(int deg = -1) const {
        assert(this->empty() || (*this)[0] == mint(0));
        if (deg == -1) deg = int(this->size());
        fps res{1};
        for (int m = 1; m < deg; m <<= 1) {
            fps f = pre(2 * m);
            fps g = res.log(2 * m);
            f.resize(2 * m);
            for (int i = 0; i < 2 * m; i++) f[i] -= g[i];
            f[0] += 1;
            res = (res * f).pre(2 * m);
        }
        return res.pre(deg);
    }

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
            fps f = pre(2 * m);
            res = ((res + f * res.inv(2 * m)) * inv2).pre(2 * m);
        }
        return res.pre(deg);
    }
};

template <class mint>
using fps = formal_power_series<mint>;
