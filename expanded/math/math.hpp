#line 2 "src/math/math.hpp"

#include <vector>
#include <algorithm>
#include <utility>
#include <cassert>

// 高速冪乗算: (x^n) % m
long long pow_mod(long long x, long long n, int m) {
    assert(0 <= n);
    assert(1 <= m);
    if (m == 1) return 0;
    unsigned int _m = (unsigned int)(m);
    unsigned long long r = 1;
    unsigned long long base = (x % m + m) % m;
    while (n > 0) {
        if (n & 1) r = (r * base) % _m;
        base = (base * base) % _m;
        n >>= 1;
    }
    return (long long)(r);
}

// 拡張ユークリッドの互除法: {gcd(a, b), x} s.t. ax ≡ gcd(a, b) (mod b)
constexpr std::pair<long long, long long> inv_gcd(long long a, long long b) {
    a = (a % b + b) % b;
    if (a == 0) return {b, 0};
    long long s = b, t = a;
    long long m0 = 0, m1 = 1;
    while (t) {
        long long u = s / t;
        s -= u * t;
        m0 -= u * m1;
        std::swap(s, t);
        std::swap(m0, m1);
    }
    if (m0 < 0) m0 += b / s;
    return {s, m0};
}

// 一次関数の床関数の総和: sum_{i=0}^{n-1} floor((ai + b) / m)
long long floor_sum(long long n, long long m, long long a, long long b) {
    assert(0 <= n);
    assert(1 <= m);
    long long ans = 0;
    if (a < 0) {
        long long a2 = (a % m + m) % m;
        ans -= n * (n - 1) / 2 * ((a2 - a) / m);
        a = a2;
    }
    if (b < 0) {
        long long b2 = (b % m + m) % m;
        ans -= n * ((b2 - b) / m);
        b = b2;
    }
    
    while (true) {
        if (a >= m) {
            ans += n * (n - 1) / 2 * (a / m);
            a %= m;
        }
        if (b >= m) {
            ans += n * (b / m);
            b %= m;
        }

        long long y_max = a * n + b;
        if (y_max < m) break;
        n = y_max / m;
        b = y_max % m;
        std::swap(m, a);
    }
    return ans;
}

// 中国剰余定理 (Chinese Remainder Theorem)
// x ≡ r_i (mod m_i) を満たす {rem, mod} を返す。解なしは {0, 0}
std::pair<long long, long long> crt(const std::vector<long long>& r, const std::vector<long long>& m) {
    assert(r.size() == m.size());
    int n = int(r.size());
    long long r0 = 0, m0 = 1;
    for (int i = 0; i < n; i++) {
        assert(1 <= m[i]);
        long long r1 = (r[i] % m[i] + m[i]) % m[i], m1 = m[i];
        if (m0 < m1) {
            std::swap(r0, r1);
            std::swap(m0, m1);
        }
        if (m0 % m1 == 0) {
            if (r0 % m1 != r1) return {0, 0};
            continue;
        }
        auto eg = inv_gcd(m0, m1);
        long long g = eg.first;
        long long im0 = eg.second;
        if ((r1 - r0) % g != 0) return {0, 0};
        long long u1 = (m1 / g);
        long long x = (r1 - r0) / g % u1 * im0 % u1;
        r0 += x * m0;
        m0 *= u1;
        if (r0 < 0) r0 += m0;
    }
    return {r0, m0};
}
