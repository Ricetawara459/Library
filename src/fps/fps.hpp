#pragma once

#include <algorithm>
#include <cassert>
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
