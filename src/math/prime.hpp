#pragma once

#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <chrono>
#include <random>

namespace internal {

// 巨大な数用の高速な冪乗算
constexpr unsigned long long pow_mod_u128(unsigned long long x, unsigned long long n, unsigned long long m) {
    unsigned long long r = 1;
    unsigned __int128 base = x % m;
    while (n > 0) {
        if (n & 1) r = (unsigned long long)((unsigned __int128)r * base % m);
        base = (unsigned __int128)base * base % m;
        n >>= 1;
    }
    return r;
}

// Miller-Rabin 素数判定法
constexpr bool miller_rabin(unsigned long long n) {
    if (n <= 1) return false;
    if (n == 2 || n == 3 || n == 5 || n == 7) return true;
    if (n % 2 == 0 || n % 3 == 0 || n % 5 == 0 || n % 7 == 0) return false;
    if (n < 121) return true;

    unsigned long long d = n - 1;
    int s = 0;
    while ((d & 1) == 0) {
        d >>= 1;
        s++;
    }

    constexpr unsigned long long bases[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
    for (unsigned long long a : bases) {
        if (a % n == 0) continue;
        unsigned long long x = pow_mod_u128(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool composite = true;
        for (int r = 1; r < s; r++) {
            x = (unsigned long long)((unsigned __int128)x * x % n);
            if (x == n - 1) {
                composite = false;
                break;
            }
        }
        if (composite) return false;
    }
    return true;
}

// 修正版 Pollard's rho アルゴリズム
long long pollard_rho(long long n) {
    if (n % 2 == 0) return 2;
    if (miller_rabin(n)) return n;
    
    // 撃墜ケースを回避するために実行時乱数を使用
    static std::mt19937_64 rng(std::chrono::steady_clock::now().time_since_epoch().count());
    
    while (true) {
        long long step = rng() % (n - 1) + 1;
        auto f = [&](long long x) {
            return (long long)((unsigned __int128)x * x % n + step) % n;
        };

        long long x = 0, y = 0, prd = 1;
        long long d = 1;
        
        // 128ステップごとに一括して gcd を取る（ボトルネック軽減）
        for (int i = 1; ; i++) {
            x = f(x);
            y = f(f(y));
            long long diff = std::abs(x - y);
            if (diff == 0) break;
            
            long long next_prd = (long long)((unsigned __int128)prd * diff % n);
            if (next_prd == 0) {
                // prd * diff が n の倍数になった場合、
                // 潰れる前の prd または diff から非自明な因数を回収する
                d = std::gcd(prd, n);
                if (d > 1 && d < n) return d;
                d = std::gcd(diff, n);
                if (d > 1 && d < n) return d;
                break;
            }
            prd = next_prd;
            
            if (i % 128 == 0) {
                d = std::gcd(prd, n);
                if (d > 1 && d < n) return d;
                prd = 1;
            }
        }
        
        d = std::gcd(prd, n);
        if (d > 1 && d < n) return d;
    }
}

void factorize_inner(long long n, std::vector<long long>& res) {
    if (n == 1) return;
    if (miller_rabin(n)) {
        res.push_back(n);
        return;
    }
    long long d = pollard_rho(n);
    factorize_inner(d, res);
    factorize_inner(n / d, res);
}

} // namespace internal


// ==================== 1. 巨大な数用 (Pollard's rho / Miller-Rabin) ====================

constexpr bool is_prime(long long n) {
    return internal::miller_rabin(n);
}

std::vector<long long> factorize(long long n) {
    std::vector<long long> res;
    if (n <= 1) return res;
    
    // 小さな素数で予め割っておくことで、ポラード・ローの呼び出し回数を劇的に減らす
    for (long long p : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47}) {
        while (n % p == 0) {
            res.push_back(p);
            n /= p;
        }
    }
    
    internal::factorize_inner(n, res);
    std::sort(res.begin(), res.end());
    return res;
}

std::vector<long long> divisors(long long n) {
    if (n <= 0) return {};
    auto facts = factorize(n);
    std::vector<std::pair<long long, int>> pf;
    for (long long p : facts) {
        if (pf.empty() || pf.back().first != p) {
            pf.push_back({p, 1});
        } else {
            pf.back().second++;
        }
    }
    std::vector<long long> res;
    auto dfs = [&](auto self, int idx, long long val) -> void {
        if (idx == (int)pf.size()) {
            res.push_back(val);
            return;
        }
        long long mul = 1;
        for (int i = 0; i <= pf[idx].second; i++) {
            self(self, idx + 1, val * mul);
            if (i < pf[idx].second) mul *= pf[idx].first;
        }
    };
    dfs(dfs, 0, 1);
    std::sort(res.begin(), res.end());
    return res;
}


// ==================== 2. 動的エラトステネスの篩 ====================

struct sieve {
  public:
    sieve() = default;

    void extend(int n) {
        if (n <= _max_n) return;
        int next_n = _max_n;
        while (next_n < n) next_n *= 2;

        int old_n = _max_n;
        _max_n = next_n;

        _min_factor.resize(_max_n + 1);
        for (int i = old_n + 1; i <= _max_n; i++) {
            _min_factor[i] = i;
        }

        for (int i = 2; i * i <= _max_n; i++) {
            if (_min_factor[i] == i) {
                int start = ((old_n + i) / i) * i;
                start = std::max(start, i * i);
                for (int j = start; j <= _max_n; j += i) {
                    if (_min_factor[j] == j) {
                        _min_factor[j] = i;
                    }
                }
            }
        }

        _mobius.resize(_max_n + 1);
        for (int i = old_n + 1; i <= _max_n; i++) {
            int p = _min_factor[i];
            if (i == p) {
                _mobius[i] = -1;
            } else if ((i / p) % p == 0) {
                _mobius[i] = 0;
            } else {
                _mobius[i] = -_mobius[i / p];
            }
        }
    }

    std::vector<int> factorize(int x) {
        assert(0 <= x && x <= _max_n);
        std::vector<int> res;
        while (x > 1) {
            res.push_back(_min_factor[x]);
            x /= _min_factor[x];
        }
        return res;
    }

    std::vector<int> divisors(int x) {
        assert(0 <= x && x <= _max_n);
        if (x <= 0) return {};
        std::vector<int> res = {1};
        while (x > 1) {
            int p = _min_factor[x];
            int cnt = 0;
            while (x % p == 0) {
                cnt++;
                x /= p;
            }
            int sz = int(res.size());
            int mul = 1;
            for (int i = 0; i < cnt; i++) {
                mul *= p;
                for (int j = 0; j < sz; j++) {
                    res.push_back(res[j] * mul);
                }
            }
        }
        std::sort(res.begin(), res.end());
        return res;
    }

    int mobius(int x) {
        assert(0 <= x && x <= _max_n);
        return _mobius[x];
    }

    bool is_prime(int x) {
        assert(0 <= x && x <= _max_n);
        return x >= 2 && _min_factor[x] == x;
    }

  private:
    int _max_n = 1;
    std::vector<int> _min_factor = {0, 1};
    std::vector<int> _mobius = {0, 1};
};

inline sieve default_sieve;

inline std::vector<int> factorize_sieve(int n) {
    default_sieve.extend(n);
    return default_sieve.factorize(n);
}

inline std::vector<int> divisors_sieve(int n) {
    default_sieve.extend(n);
    return default_sieve.divisors(n);
}

inline int mobius_sieve(int n) {
    default_sieve.extend(n);
    return default_sieve.mobius(n);
}

inline bool is_prime_sieve(int n) {
    default_sieve.extend(n);
    return default_sieve.is_prime(n);
}