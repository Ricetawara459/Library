#ifndef LIBRARY_FPS_FPS_HPP
#define LIBRARY_FPS_FPS_HPP

#pragma once

#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>

#include "../math/convolution.hpp"

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
