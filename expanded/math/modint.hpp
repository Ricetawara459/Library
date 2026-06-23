#line 2 "src/math/modint.hpp"

#include <iostream>

template <int MOD>
struct static_modint {
    using mint = static_modint;
  public:
    long long val;
    static_modint() : val(0) {}
    static_modint(long long v) {
        val = v % MOD;
        if (val < 0) val += MOD;
    }
    static int mod() { return MOD; }
    mint pow(long long n) const {
        mint res(1), mul(val);
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
        val += a.val;
        if (val >= MOD) val -= MOD;
        return *this;
    }
    mint& operator-=(const mint& a) {
        val -= a.val;
        if (val < 0) val += MOD;
        return *this;
    }
    mint& operator*=(const mint& a) {
        val = (val * a.val) % MOD;
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
    friend bool operator==(const mint& a, const mint& b) { return a.val == b.val; }
    friend bool operator!=(const mint& a, const mint& b) { return a.val != b.val; }
    friend std::ostream& operator<<(std::ostream& os, const mint& a) { return os << a.val; }
    friend std::istream& operator>>(std::istream& is, mint& a) {
        long long v;
        is >> v;
        a = mint(v);
        return is;
    }
};

using mint998244353 = static_modint<998244353>;
