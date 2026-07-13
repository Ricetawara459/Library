#ifndef LIBRARY_MATH_BINOMIAL_PREFIX_SUM_HPP
#define LIBRARY_MATH_BINOMIAL_PREFIX_SUM_HPP

#pragma once

#include <cassert>
#include "factorial.hpp"

template <class mint>
struct binomial_prefix_sum {
  public:
    binomial_prefix_sum() : binomial_prefix_sum(0, 0) {}
    /// 初期状態 (n, r) を作る。0 <= r <= n。値は sum_{i=0}^{r-1} C(n,i)。
    binomial_prefix_sum(int n, int r) { reset(n, r); }

    /// 状態を (n, r) に作り直す。0 <= r <= n。
    void reset(int n, int r) {
        assert(0 <= r && r <= n);
        _n = n;
        _r = r;
        fac.extend(n);
        _comb = fac.comb(n, r);
        _sum = 0;
        for (int i = 0; i < r; i++) _sum += fac.comb(n, i);
    }

    /// 現在の n を返す。
    int n() const { return _n; }
    /// 現在の r を返す。
    int r() const { return _r; }
    /// 現在の C(n,r) を返す。
    mint comb() const { return _comb; }
    /// 現在の sum_{i=0}^{r-1} C(n,i) を返す。
    mint get() const { return _sum; }

    /// r を 1 増やす。F(n,r+1) = F(n,r) + C(n,r)。
    void inc_r() {
        assert(_r < _n);
        _sum += _comb;
        _comb *= mint(_n - _r);
        _r++;
        _comb /= mint(_r);
    }

    /// r を 1 減らす。F(n,r-1) = F(n,r) - C(n,r-1)。
    void dec_r() {
        assert(_r > 0);
        _comb *= mint(_r);
        _r--;
        _comb /= mint(_n - _r);
        _sum -= _comb;
    }

    /// n を 1 増やす。F(n+1,r) = 2F(n,r) - C(n,r-1)。
    void inc_n() {
        mint prev = (_r == 0 ? mint(0) : _comb * mint(_r) / mint(_n - _r + 1));
        _sum += _sum - prev;
        _n++;
        _comb *= mint(_n);
        _comb /= mint(_n - _r);
    }

    /// n を 1 減らす。r <= n-1 が必要。
    void dec_n() {
        assert(_n > 0);
        assert(_r < _n);
        mint prev = (_r == 0 ? mint(0) : _comb * mint(_r) / mint(_n));
        _sum = (_sum + prev) / mint(2);
        _comb *= mint(_n - _r);
        _comb /= mint(_n);
        _n--;
    }

    /// 隣接移動を繰り返して (n, r) に移動する。移動距離に比例する。
    void move(int n, int r) {
        assert(0 <= r && r <= n);
        while (_r > r) dec_r();
        while (_n < n) inc_n();
        while (_r < r) inc_r();
        while (_n > n) dec_n();
    }

  private:
    int _n = 0;
    int _r = 0;
    mint _sum = 0;   // sum_{i=0}^{r-1} C(n,i)
    mint _comb = 1;  // C(n,r)
    factorial<mint> fac;
};

#endif  // LIBRARY_MATH_BINOMIAL_PREFIX_SUM_HPP
