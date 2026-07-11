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
