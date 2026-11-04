// Source: https://github.com/Ricetawara459/Library
#ifndef LIBRARY_MISC_MO_HPP
#define LIBRARY_MISC_MO_HPP



#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <vector>

// Mo's algorithm のクエリ順序管理。非負整数の2変数 (x, y) を独立に動かす。
// 2変数クエリ: run_xy(inc_x, inc_y, dec_x, dec_y, output)。初期状態は (0, 0)。
// 半開区間 [l, r): run(add, erase, output) または左右別の5引数版を使う。
struct mo_algorithm {
  public:
    mo_algorithm() = default;
    /// 両座標の上限を n にする。区間クエリでは配列の長さを渡す。
    explicit mo_algorithm(int n) : mo_algorithm(n, n) {}
    /// 座標の上限をそれぞれ nx, ny にする。状態は (0, 0) から開始する。
    mo_algorithm(int nx, int ny) : _nx(nx), _ny(ny) {
        assert(nx >= 0 && ny >= 0);
    }

    /// (x, y) を追加する。0 <= x <= nx, 0 <= y <= ny。x <= y は不要。
    /// 区間クエリでは半開区間の端点 (l, r) を渡す。
    void add_query(int x, int y) {
        assert(0 <= x && x <= _nx && 0 <= y && y <= _ny);
        xs.push_back(x);
        ys.push_back(y);
    }

    /// 2変数をMo順に動かす。inc_x(x++), inc_y(y++), dec_x(--x), dec_y(--y) を呼ぶ。
    /// 呼び出し側は (0, 0) の状態を用意する。到達時に output(登録順の0-indexed番号) を呼ぶ。
    template <class IncX, class IncY, class DecX, class DecY, class Output>
    void run_xy(IncX inc_x, IncY inc_y, DecX dec_x, DecY dec_y, Output output) const {
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

    /// 半開区間 [l, r) 用。左追加・右追加・左削除・右削除の順に関数を渡す。l <= r が必要。
    template <class AddLeft, class AddRight, class EraseLeft, class EraseRight, class Output>
    void run(AddLeft add_left, AddRight add_right, EraseLeft erase_left, EraseRight erase_right, Output output) const {
        run_xy(erase_left, add_right, add_left, erase_right, output);
    }

    /// 半開区間 [l, r) 用。追加/削除が左右で同じ場合の簡易版。
    template <class Add, class Erase, class Output>
    void run(Add add, Erase erase, Output output) const {
        run(add, add, erase, erase, output);
    }

    /// 追加したクエリ数を返す。
    int size() const { return int(xs.size()); }

  private:
    int _nx = 0, _ny = 0;
    std::vector<int> xs, ys;
};

#endif  // LIBRARY_MISC_MO_HPP
