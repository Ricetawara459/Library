#pragma once

#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <cassert>

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

long long pollard_rho(long long n) {
    if (n % 2 == 0) return 2;
    if (miller_rabin(n)) return n;
    
    long long step = 0;
    auto f = [&](long long x) {
        return (long long)((unsigned __int128)x * x % n + step) % n;
    };

    while (true) {
        step++;
        long long x = 2, y = 2, d = 1;
        long long prd = 1;
        int i = 0;
        while (d == 1) {
            x = f(x);
            y = f(f(y));
            long long diff = std::abs(x - y);
            if (diff == 0) break;
            prd = (long long)((unsigned __int128)prd * diff % n);
            if (prd == 0) break;
            i++;
            if (i == 127) {
                d = std::gcd(prd, n);
                if (d > 1) return d;
                i = 0;
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

    // 必要に応じてサイズを倍々に拡張する
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

        // 既存の素数および新規素数で拡張区間を塗りつぶす
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
    std::vector<int> _min_factor = {0, 1}; // private隠蔽
    std::vector<int> _mobius = {0, 1};     // private隠蔽
};

// ヘッダ内で共有されるグローバルなインスタンス
inline sieve default_sieve;

// インスタンスを意識せず直接呼べるトップレベル関数群（自動拡張付き）
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