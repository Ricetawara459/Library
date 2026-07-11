#ifndef LIBRARY_DATA_STRUCTURE_SEGTREE_HPP
#define LIBRARY_DATA_STRUCTURE_SEGTREE_HPP

#line 1 "src/data-structure/segtree.hpp"


#ifndef LIBRARY_DATA_STRUCTURE_SEGMENT_TREE_HPP
#define LIBRARY_DATA_STRUCTURE_SEGMENT_TREE_HPP

#line 1 "src/data-structure/segment_tree.hpp"


#include <vector>
#include <cassert>

template <class S, S (*op)(S, S), S (*e)()> 
struct segtree {
  public:
    segtree() : segtree(0) {}
    /// 長さ n、全要素 e() のセグメント木を作る。
    explicit segtree(int n) : segtree(std::vector<S>(n, e())) {}
    /// 配列 v からセグメント木を作る。
    explicit segtree(const std::vector<S>& v) : _n(int(v.size())) {
        _log = 0;
        size = 1;
        while (size < _n) {
            size <<= 1;
            _log++;
        }

        d = std::vector<S>(2 * size, e());
        for (int i = 0; i < _n; i++) d[size + i] = v[i];
        for (int i = size - 1; i >= 1; i--) {
            update(i);
        }
    }

    /// 一点更新。a[p] = x。
    void set(int p, S x) {
        assert(0 <= p && p < _n);
        p += size;
        d[p] = x;
        for (int i = 1; i <= _log; i++) update(p >> i);
    }

    /// a[p] を返す。
    S get(int p) const {
        assert(0 <= p && p < _n);
        return d[p + size];
    }

    /// 半開区間 [l, r) を op で畳み込む。
    S prod(int l, int r) const {
        assert(0 <= l && l <= r && r <= _n);
        S sml = e(), smr = e();
        l += size;
        r += size;

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

    /// l から右に伸ばして、条件 f が true である最大の右端を返す。
    template <bool (*f)(S)> int max_right(int l) const {
        return max_right(l, [](S x) { return f(x); });
    }
    /// l から右に伸ばして、条件 f が true である最大の右端を返す。
    template <class F> int max_right(int l, F f) const {
        assert(0 <= l && l <= _n);
        assert(f(e()));
        if (l == _n) return _n;
        l += size;
        S sm = e();
        do {
            while (l % 2 == 0) l >>= 1;
            if (!f(op(sm, d[l]))) {
                while (l < size) {
                    l = (2 * l);
                    if (f(op(sm, d[l]))) {
                        sm = op(sm, d[l]);
                        l++;
                    }
                }
                return l - size;
            }
            sm = op(sm, d[l]);
            l++;
        } while ((l & -l) != l);
        return _n;
    }

    /// r から左に伸ばして、条件 f が true である最小の左端を返す。
    template <bool (*f)(S)> int min_left(int r) const {
        return min_left(r, [](S x) { return f(x); });
    }
    /// r から左に伸ばして、条件 f が true である最小の左端を返す。
    template <class F> int min_left(int r, F f) const {
        assert(0 <= r && r <= _n);
        assert(f(e()));
        if (r == 0) return 0;
        r += size;
        S sm = e();
        do {
            r--;
            while (r > 1 && (r % 2)) r >>= 1;
            if (!f(op(d[r], sm))) {
                while (r < size) {
                    r = (2 * r + 1);
                    if (f(op(d[r], sm))) {
                        sm = op(d[r], sm);
                        r--;
                    }
                }
                return r + 1 - size;
            }
            sm = op(d[r], sm);
        } while ((r & -r) != r);
        return 0;
    }

  private:
    int _n, size, _log;
    std::vector<S> d;

    void update(int k) { d[k] = op(d[2 * k], d[2 * k + 1]); }
};

#endif  // LIBRARY_DATA_STRUCTURE_SEGMENT_TREE_HPP

#endif  // LIBRARY_DATA_STRUCTURE_SEGTREE_HPP
