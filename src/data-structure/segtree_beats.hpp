#ifndef LIBRARY_DATA_STRUCTURE_SEGTREE_BEATS_HPP
#define LIBRARY_DATA_STRUCTURE_SEGTREE_BEATS_HPP

#pragma once

#include <algorithm>
#include <cassert>
#include <vector>

// 抽象化 Segment Tree Beats。
//
// 使い方:
//   S: 区間情報。op(S, S) で左右の子をマージできるようにする。
//   F: 区間に作用させるクエリ情報。
//   mapping(f, x): ノード x 全体に作用 f をまとめて適用できるなら x を更新して true。
//                  まとめて適用できないなら x を変更せず false。false の時だけ子へ降りる。
//   composition(f, g): 既存の遅延 g の後に新しい作用 f を行う合成、つまり f(g(x))。
//   id(): 何もしない作用。
//
// 例: range chmin/chmax/add/sum は snippets/segtree_beats_chmin_chmax_sum.cpp を参照。
// 注意: Beats の償却計算量は S と mapping の設計に依存する。
template <class S,
          S (*op)(S, S),
          S (*e)(),
          class F,
          bool (*mapping)(F, S&),
          F (*composition)(F, F),
          F (*id)()>
struct segtree_beats {
  public:
    segtree_beats() : segtree_beats(0) {}
    /// 長さ n、全要素 e() の抽象化 Segment Tree Beats を作る。
    explicit segtree_beats(int n) : segtree_beats(std::vector<S>(n, e())) {}
    /// 配列 v から抽象化 Segment Tree Beats を作る。
    explicit segtree_beats(const std::vector<S>& v) : _n(int(v.size())) {
        size = 1;
        while (size < _n) size <<= 1;
        d = std::vector<S>(2 * size, e());
        lz = std::vector<F>(size, id());
        for (int i = 0; i < _n; i++) d[size + i] = v[i];
        for (int i = size - 1; i >= 1; i--) update(i);
    }

    /// 一点更新。a[p] = x。
    void set(int p, S x) {
        assert(0 <= p && p < _n);
        set_rec(p, x, 1, 0, size);
    }

    /// 遅延を反映して a[p] を返す。
    S get(int p) {
        assert(0 <= p && p < _n);
        return prod(p, p + 1);
    }

    /// 半開区間 [l, r) を op で畳み込む。
    S prod(int l, int r) {
        assert(0 <= l && l <= r && r <= _n);
        return prod_rec(l, r, 1, 0, size);
    }

    /// 全区間 [0, n) の畳み込み結果を返す。
    S all_prod() { return prod(0, _n); }

    /// 一点に作用 f を適用する。
    void apply(int p, F f) {
        assert(0 <= p && p < _n);
        apply(p, p + 1, f);
    }

    /// 半開区間 [l, r) に作用 f を適用する。mapping が失敗したノードだけ子へ降りる。
    void apply(int l, int r, F f) {
        assert(0 <= l && l <= r && r <= _n);
        if (l == r) return;
        apply_rec(l, r, f, 1, 0, size);
    }

    /// l から右に伸ばして、条件 g が true である最大の右端を返す。
    template <bool (*g)(S)> int max_right(int l) {
        return max_right(l, [](S x) { return g(x); });
    }
    /// l から右に伸ばして、条件 g が true である最大の右端を返す。
    template <class G> int max_right(int l, G g) {
        assert(0 <= l && l <= _n);
        assert(g(e()));
        if (l == _n) return _n;
        S sm = e();
        int res = max_right_rec(l, g, sm, 1, 0, size);
        return std::min(res, _n);
    }

    /// r から左に伸ばして、条件 g が true である最小の左端を返す。
    template <bool (*g)(S)> int min_left(int r) {
        return min_left(r, [](S x) { return g(x); });
    }
    /// r から左に伸ばして、条件 g が true である最小の左端を返す。
    template <class G> int min_left(int r, G g) {
        assert(0 <= r && r <= _n);
        assert(g(e()));
        if (r == 0) return 0;
        S sm = e();
        return min_left_rec(r, g, sm, 1, 0, size);
    }

  private:
    int _n, size;
    std::vector<S> d;
    std::vector<F> lz;

    void update(int k) { d[k] = op(d[2 * k], d[2 * k + 1]); }

    bool all_apply(int k, F f) {
        if (!mapping(f, d[k])) return false;
        if (k < size) lz[k] = composition(f, lz[k]);
        return true;
    }

    void push(int k) {
        if (k >= size) return;
        bool left_ok = all_apply(2 * k, lz[k]);
        bool right_ok = all_apply(2 * k + 1, lz[k]);
        if (!left_ok || !right_ok) assert(false);
        lz[k] = id();
    }

    void set_rec(int p, S x, int k, int l, int r) {
        if (r - l == 1) {
            d[k] = x;
            return;
        }
        push(k);
        int m = (l + r) >> 1;
        if (p < m) set_rec(p, x, 2 * k, l, m);
        else set_rec(p, x, 2 * k + 1, m, r);
        update(k);
    }

    S prod_rec(int a, int b, int k, int l, int r) {
        if (r <= a || b <= l) return e();
        if (a <= l && r <= b) return d[k];
        push(k);
        int m = (l + r) >> 1;
        return op(prod_rec(a, b, 2 * k, l, m),
                  prod_rec(a, b, 2 * k + 1, m, r));
    }

    void apply_rec(int a, int b, F f, int k, int l, int r) {
        if (r <= a || b <= l) return;
        if (a <= l && r <= b && all_apply(k, f)) return;
        if (r - l == 1) {
            assert(false);
            return;
        }
        push(k);
        int m = (l + r) >> 1;
        apply_rec(a, b, f, 2 * k, l, m);
        apply_rec(a, b, f, 2 * k + 1, m, r);
        update(k);
    }

    template <class G>
    int max_right_rec(int a, G& g, S& sm, int k, int l, int r) {
        if (r <= a) return _n;
        if (a <= l) {
            S nxt = op(sm, d[k]);
            if (g(nxt)) {
                sm = nxt;
                return _n;
            }
            if (r - l == 1) return l;
        }
        push(k);
        int m = (l + r) >> 1;
        int res = max_right_rec(a, g, sm, 2 * k, l, m);
        if (res != _n) return res;
        return max_right_rec(a, g, sm, 2 * k + 1, m, r);
    }

    template <class G>
    int min_left_rec(int b, G& g, S& sm, int k, int l, int r) {
        if (b <= l) return 0;
        if (r <= b) {
            S nxt = op(d[k], sm);
            if (g(nxt)) {
                sm = nxt;
                return 0;
            }
            if (r - l == 1) return r;
        }
        push(k);
        int m = (l + r) >> 1;
        int res = min_left_rec(b, g, sm, 2 * k + 1, m, r);
        if (res != 0) return res;
        return min_left_rec(b, g, sm, 2 * k, l, m);
    }
};

#endif  // LIBRARY_DATA_STRUCTURE_SEGTREE_BEATS_HPP
