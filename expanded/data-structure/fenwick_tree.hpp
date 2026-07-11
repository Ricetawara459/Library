#ifndef LIBRARY_DATA_STRUCTURE_FENWICK_TREE_HPP
#define LIBRARY_DATA_STRUCTURE_FENWICK_TREE_HPP

#line 1 "src/data-structure/fenwick_tree.hpp"

#include <vector>
#include <cassert>

template <class T>
struct fenwick_tree {
  public:
    fenwick_tree() : _n(0) {}
    /// 長さ n の Fenwick Tree を作る。
    explicit fenwick_tree(int n) : _n(n), data(n, 0) {}

    /// 一点加算。a[p] += x。
    void add(int p, T x) {
        assert(0 <= p && p < _n);
        p++;
        while (p <= _n) {
            data[p - 1] += x;
            p += p & -p;
        }
    }

    /// 半開区間 [l, r) の総和を返す。
    T sum(int l, int r) const {
        assert(0 <= l && l <= r && r <= _n);
        return sum(r) - sum(l);
    }

  private:
    int _n;
    std::vector<T> data;

    T sum(int r) const {
        T s = 0;
        while (r > 0) {
            s += data[r - 1];
            r -= r & -r;
        }
        return s;
    }
};

#endif  // LIBRARY_DATA_STRUCTURE_FENWICK_TREE_HPP
