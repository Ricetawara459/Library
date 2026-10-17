// Source: https://github.com/Ricetawara459/Library
#ifndef LIBRARY_DATA_STRUCTURE_FENWICK_TREE_HPP
#define LIBRARY_DATA_STRUCTURE_FENWICK_TREE_HPP



#include <vector>
#include <cassert>
#include <ostream>

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

    /// デバッグ用。元の配列を [x, y, ...] 形式で出力する。T の << が必要。
    /// コピーした区間和から各要素を復元する。時間・追加領域 O(n)。
    friend std::ostream& operator<<(std::ostream& os, const fenwick_tree& tree) {
        auto values = tree.data;
        for (int i = tree._n; i > 0; i--) {
            int parent = i + (i & -i);
            if (parent <= tree._n) values[parent - 1] = values[parent - 1] - values[i - 1];
        }
        os << '[';
        for (int i = 0; i < tree._n; i++) {
            if (i) os << ", ";
            os << values[i];
        }
        return os << ']';
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
