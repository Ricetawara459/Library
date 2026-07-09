#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <vector>

// Mo's algorithm のクエリ順序管理。区間は半開区間 [l, r)。
// 左右どちらから追加/削除しても同じ処理でよい場合は run(add, erase, output) を使う。
// 左端と右端で処理を分けたい場合は run(add_left, add_right, erase_left, erase_right, output) を使う。
struct mo {
  public:
    mo() = default;
    explicit mo(int n) : _n(n) {}

    void add_query(int l, int r) {
        assert(0 <= l && l <= r && r <= _n);
        left.push_back(l);
        right.push_back(r);
    }

    template <class AddLeft, class AddRight, class EraseLeft, class EraseRight, class Output>
    void run(AddLeft add_left, AddRight add_right, EraseLeft erase_left, EraseRight erase_right, Output output) const {
        int q = int(left.size());
        if (q == 0) return;
        int block = std::max(1, int(_n / std::max(1.0, std::sqrt(double(q)))));
        std::vector<int> ord(q);
        std::iota(ord.begin(), ord.end(), 0);
        std::sort(ord.begin(), ord.end(), [&](int a, int b) {
            int ab = left[a] / block;
            int bb = left[b] / block;
            if (ab != bb) return ab < bb;
            if (ab & 1) return right[a] > right[b];
            return right[a] < right[b];
        });

        int nl = 0, nr = 0;
        for (int id : ord) {
            while (left[id] < nl) add_left(--nl);
            while (nr < right[id]) add_right(nr++);
            while (nl < left[id]) erase_left(nl++);
            while (right[id] < nr) erase_right(--nr);
            output(id);
        }
    }

    template <class Add, class Erase, class Output>
    void run(Add add, Erase erase, Output output) const {
        run(add, add, erase, erase, output);
    }

    int size() const { return int(left.size()); }

  private:
    int _n = 0;
    std::vector<int> left, right;
};
