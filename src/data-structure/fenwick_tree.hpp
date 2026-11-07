#ifndef LIBRARY_DATA_STRUCTURE_FENWICK_TREE_HPP
#define LIBRARY_DATA_STRUCTURE_FENWICK_TREE_HPP

#pragma once

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

    /// sum(0, i+1) >= x となる最小の i を返す。存在しなければ n。x <= 0 なら 0。
    /// 全要素が非負で、T が順序比較可能であること。計算量 O(log n)。
    int lower_bound(T x) const {
        if (!(T(0) < x)) return 0;
        int pos = 0, step = 1;
        while (step <= _n / 2) step *= 2;
        T acc = 0;
        for (; step > 0; step /= 2) {
            if (step <= _n - pos && acc + data[pos + step - 1] < x) {
                acc += data[pos + step - 1];
                pos += step;
            }
        }
        return pos;
    }

    /// sum(0, i+1) > x となる最小の i を返す。存在しなければ n。x < 0 なら 0。
    /// 全要素が非負で、T が順序比較可能であること。計算量 O(log n)。
    int upper_bound(T x) const {
        if (x < T(0)) return 0;
        int pos = 0, step = 1;
        while (step <= _n / 2) step *= 2;
        T acc = 0;
        for (; step > 0; step /= 2) {
            if (step <= _n - pos && !(x < acc + data[pos + step - 1])) {
                acc += data[pos + step - 1];
                pos += step;
            }
        }
        return pos;
    }

    /// 元の配列を空白区切りで出力する。末尾空白・改行なし。T の << が必要。
    /// コピーした区間和から各要素を復元する。時間・追加領域 O(n)。
    friend std::ostream& operator<<(std::ostream& os, const fenwick_tree& tree) {
        auto values = tree.data;
        for (int i = tree._n; i > 0; i--) {
            int parent = i + (i & -i);
            if (parent <= tree._n) values[parent - 1] = values[parent - 1] - values[i - 1];
        }
        for (int i = 0; i < tree._n; i++) {
            if (i) os << ' ';
            os << values[i];
        }
        return os;
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
