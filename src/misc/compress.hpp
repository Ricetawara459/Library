#pragma once

#include <vector>
#include <algorithm>
#include <cassert>

template <typename T>
struct compress {
  public:
    compress() = default;
    explicit compress(const std::vector<T>& a) { build(a); }

    // ベクターを元に座標圧縮のテーブルを構築
    void build(const std::vector<T>& a) {
        _xs = a;
        std::sort(_xs.begin(), _xs.end());
        _xs.erase(std::unique(_xs.begin(), _xs.end()), _xs.end());
    }

    // 元の値 x が圧縮後に何番目になるか (0-indexed)
    int operator()(const T& x) const {
        auto it = std::lower_bound(_xs.begin(), _xs.end(), x);
        assert(it != _xs.end() && *it == x);
        return std::distance(_xs.begin(), it);
    }

    // テーブルに存在しない値かもしれないが、仮に挿入した場合の lower_bound の位置を返す
    int get_id(const T& x) const {
        return std::distance(_xs.begin(), std::lower_bound(_xs.begin(), _xs.end(), x));
    }

    // 圧縮後のインデックス id から元の値を取得
    T uncompress(int id) const {
        assert(0 <= id && id < size());
        return _xs[id];
    }

    // [] 演算子でも元の値を取得できるようにする
    T operator[](int id) const { return uncompress(id); }

    // 圧縮後の種類数
    int size() const { return int(_xs.size()); }

  private:
    std::vector<T> _xs;
};