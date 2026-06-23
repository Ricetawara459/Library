#line 2 "src/math/modint.hpp"

#include <iostream>

template <int MOD>
struct static_modint {
    using mint = static_modint;
  private:
    long long _v; // 変数名を _v に変更し、private に隠蔽
  public:
    static_modint() : _v(0) {}
    static_modint(long long v) {
        _v = v % MOD;
        if (_v < 0) _v += MOD;
    }
    
    // ★ ACL互換のゲッター関数を追加
    int val() const { return static_cast<int>(_v); }
    
    static int mod() { return MOD; }
    mint pow(long long n) const {
        mint res(1), mul(_v);
        while (n > 0) {
            if (n & 1) res *= mul;
            mul *= mul;
            n >>= 1;
        }
        return res;
    }
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
        _v = (_v * a._v) % MOD;
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
