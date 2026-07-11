#ifndef LIBRARY_MATH_FACTORIAL_HPP
#define LIBRARY_MATH_FACTORIAL_HPP

#line 1 "src/math/factorial.hpp"


#include <vector>

template <class mint>
struct factorial {
  public:
    factorial() : _n(0), _fact(1, mint(1)), _inv_fact(1, mint(1)) {}
    explicit factorial(int n) : factorial() { extend(n); }

    /// n! まで前計算を拡張する。
    void extend(int n) {
        if (n <= _n) return;
        int next_n = _n == 0 ? 1 : _n;
        while (next_n < n) next_n *= 2;
        
        _fact.resize(next_n + 1);
        _inv_fact.resize(next_n + 1);
        
        // 新しく増えた区間の階乗を計算
        for (int i = _n + 1; i <= next_n; i++) {
            _fact[i] = _fact[i - 1] * i;
        }
        // 新しい末尾の逆元から逆順に埋める
        _inv_fact[next_n] = _fact[next_n].inv();
        for (int i = next_n - 1; i > _n; i--) {
            _inv_fact[i] = _inv_fact[i + 1] * (i + 1);
        }
        _n = next_n;
    }

    /// n! を返す。n < 0 なら 0。
    mint fact(int n) {
        if (n < 0) return mint(0);
        extend(n);
        return _fact[n];
    }

    /// 1 / n! を返す。n < 0 なら 0。
    mint fact_inv(int n) {
        if (n < 0) return mint(0);
        extend(n);
        return _inv_fact[n];
    }

    /// 順列 nPr を返す。不正な範囲なら 0。
    mint perm(int n, int r) {
        if (r < 0 || n < r) return mint(0);
        extend(n);
        return _fact[n] * _inv_fact[n - r];
    }

    /// 組合せ nCr を返す。不正な範囲なら 0。
    mint comb(int n, int r) {
        if (r < 0 || n < r) return mint(0);
        extend(n);
        return _fact[n] * _inv_fact[r] * _inv_fact[n - r];
    }

    /// 重複組合せ nHr を返す。不正な範囲なら 0。
    mint homo(int n, int r) {
        if (n < 0 || r < 0) return mint(0);
        if (n == 0 && r == 0) return mint(1);
        return comb(n + r - 1, r);
    }

  private:
    int _n;
    std::vector<mint> _fact;     // privateに隠蔽
    std::vector<mint> _inv_fact; // privateに隠蔽
};

#endif  // LIBRARY_MATH_FACTORIAL_HPP
