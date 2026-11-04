// Source: https://github.com/Ricetawara459/Library
#ifndef LIBRARY_MISC_MO_2D_HPP
#define LIBRARY_MISC_MO_2D_HPP



#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <vector>

// 独立した非負整数の2変数 (x, y) に対するMo。x <= y の制約はない。
// mo_algorithm_2d mo(nx, ny); mo.add_query(x, y);
// mo.run(inc_x, inc_y, dec_x, dec_y, output); 初期状態は (0, 0)。
struct mo_algorithm_2d {
  public:
    mo_algorithm_2d() = default;
    /// 座標の上限をそれぞれ nx, ny にする。
    mo_algorithm_2d(int nx, int ny) : _nx(nx), _ny(ny) {
        assert(nx >= 0 && ny >= 0);
    }

    /// (x, y) を追加する。0 <= x <= nx, 0 <= y <= ny。x <= y は不要。
    void add_query(int x, int y) {
        assert(0 <= x && x <= _nx && 0 <= y && y <= _ny);
        xs.push_back(x);
        ys.push_back(y);
    }

    /// inc_x(x++), inc_y(y++), dec_x(--x), dec_y(--y) を呼んでMo順に動かす。
    /// 呼び出し側は (0, 0) の状態を用意する。到達時に output(登録順の0-indexed番号) を呼ぶ。
    template <class IncX, class IncY, class DecX, class DecY, class Output>
    void run(IncX inc_x, IncY inc_y, DecX dec_x, DecY dec_y, Output output) const {
        int q = int(xs.size());
        if (q == 0) return;
        int block = std::max(1, int(std::sqrt(double(_nx) * _ny / q)));
        std::vector<int> ord(q);
        std::iota(ord.begin(), ord.end(), 0);
        std::sort(ord.begin(), ord.end(), [&](int a, int b) {
            int ab = xs[a] / block;
            int bb = xs[b] / block;
            if (ab != bb) return ab < bb;
            if (ab & 1) return ys[a] > ys[b];
            return ys[a] < ys[b];
        });

        int x = 0, y = 0;
        for (int id : ord) {
            while (xs[id] < x) dec_x(--x);
            while (y < ys[id]) inc_y(y++);
            while (x < xs[id]) inc_x(x++);
            while (ys[id] < y) dec_y(--y);
            output(id);
        }
    }

    /// 追加したクエリ数を返す。
    int size() const { return int(xs.size()); }

  private:
    int _nx = 0, _ny = 0;
    std::vector<int> xs, ys;
};

#endif  // LIBRARY_MISC_MO_2D_HPP
