#ifndef LIBRARY_DATA_STRUCTURE_LAZYSEGTREE_HPP
#define LIBRARY_DATA_STRUCTURE_LAZYSEGTREE_HPP

#line 1 "src/data-structure/lazysegtree.hpp"

#include <vector>
#include <cassert>

template <class S,
          S (*op)(S, S),
          S (*e)(),
          class F,
          S (*mapping)(F, S),
          F (*composition)(F, F),
          F (*id)()>
struct lazysegtree {
  public:
    lazysegtree() : lazysegtree(0) {}
    /// 長さ n、全要素 e() の遅延セグメント木を作る。
    explicit lazysegtree(int n) : lazysegtree(std::vector<S>(n, e())) {}
    /// 配列 v から遅延セグメント木を作る。
    explicit lazysegtree(const std::vector<S>& v) : _n(int(v.size())) {
        _log = 0;
        size = 1;
        while (size < _n) {
            size <<= 1;
            _log++;
        }
        d = std::vector<S>(2 * size, e());
        lz = std::vector<F>(size, id());
        for (int i = 0; i < _n; i++) d[size + i] = v[i];
        for (int i = size - 1; i >= 1; i--) {
            update(i);
        }
    }

    /// 一点更新。a[p] = x。
    void set(int p, S x) {
        assert(0 <= p && p < _n);
        p += size;
        for (int i = _log; i >= 1; i--) push(p >> i);
        d[p] = x;
        for (int i = 1; i <= _log; i++) update(p >> i);
    }

    /// 遅延を反映して a[p] を返す。
    S get(int p) {
        assert(0 <= p && p < _n);
        p += size;
        for (int i = _log; i >= 1; i--) push(p >> i);
        return d[p];
    }

    /// 半開区間 [l, r) を op で畳み込む。
    S prod(int l, int r) {
        assert(0 <= l && l <= r && r <= _n);
        if (l == r) return e();

        l += size;
        r += size;

        for (int i = _log; i >= 1; i--) {
            if (((l >> i) << i) != l) push(l >> i);
            if (((r >> i) << i) != r) push((r - 1) >> i);
        }

        S sml = e(), smr = e();
        while (l < r) {
            if (l & 1) sml = op(sml, d[l++]);
            if (r & 1) smr = op(d[--r], smr);
            l >>= 1;
            r >>= 1;
        }

        return op(sml, smr);
    }

    /// 全区間 [0, n) の畳み込み結果を返す。
    S all_prod() const { return d[1]; }

    /// 一点に作用 f を適用する。
    void apply(int p, F f) {
        assert(0 <= p && p < _n);
        p += size;
        for (int i = _log; i >= 1; i--) push(p >> i);
        d[p] = mapping(f, d[p]);
        for (int i = 1; i <= _log; i++) update(p >> i);
    }
    /// 半開区間 [l, r) に作用 f を適用する。
    void apply(int l, int r, F f) {
        assert(0 <= l && l <= r && r <= _n);
        if (l == r) return;

        l += size;
        r += size;

        for (int i = _log; i >= 1; i--) {
            if (((l >> i) << i) != l) push(l >> i);
            if (((r >> i) << i) != r) push((r - 1) >> i);
        }

        {
            int l2 = l, r2 = r;
            while (l < r) {
                if (l & 1) all_apply(l++, f);
                if (r & 1) all_apply(--r, f);
                l >>= 1;
                r >>= 1;
            }
            l = l2;
            r = r2;
        }

        for (int i = 1; i <= _log; i++) {
            if (((l >> i) << i) != l) update(l >> i);
            if (((r >> i) << i) != r) update((r - 1) >> i);
        }
    }

  private:
    int _n, size, _log;
    std::vector<S> d;
    std::vector<F> lz;

    void update(int k) { d[k] = op(d[2 * k], d[2 * k + 1]); }
    void all_apply(int k, F f) {
        d[k] = mapping(f, d[k]);
        if (k < size) lz[k] = composition(f, lz[k]);
    }
    void push(int k) {
        all_apply(2 * k, lz[k]);
        all_apply(2 * k + 1, lz[k]);
        lz[k] = id();
    }
};

#endif  // LIBRARY_DATA_STRUCTURE_LAZYSEGTREE_HPP
