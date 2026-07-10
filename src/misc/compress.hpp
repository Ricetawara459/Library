#pragma once

#include <vector>
#include <algorithm>
#include <cassert>

template <typename T>
struct compress {
  public:
    compress() = default;
    explicit compress(const std::vector<T>& a) { build(a); }

    /// 配列 a の値を集めて座標圧縮を構築する。
    void build(const std::vector<T>& a) {
        _xs = a;
        std::sort(_xs.begin(), _xs.end());
        _xs.erase(std::unique(_xs.begin(), _xs.end()), _xs.end());
    }

    /// 値 x の圧縮後 ID を返す。x は存在すること。
    int operator()(const T& x) const {
        auto it = std::lower_bound(_xs.begin(), _xs.end(), x);
        assert(it != _xs.end() && *it == x);
        return std::distance(_xs.begin(), it);
    }

    /// 値 x の圧縮後 ID を返す。存在しない場合は挿入位置。
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

     /// 異なる値の個数を返す。
    int size() const { return int(_xs.size()); }

  private:
    std::vector<T> _xs;
};
