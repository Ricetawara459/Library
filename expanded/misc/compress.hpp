#ifndef LIBRARY_MISC_COMPRESS_HPP
#define LIBRARY_MISC_COMPRESS_HPP

#line 1 "src/misc/compress.hpp"


#include <algorithm>
#include <cassert>
#include <initializer_list>
#include <iterator>
#include <vector>

template <class T>
struct compress {
  public:
    compress() = default;
    explicit compress(const std::vector<T>& a) { build(a); }
    compress(std::initializer_list<T> a) {
        add(a.begin(), a.end());
        build();
    }

    /// 座標圧縮の候補に値 x を追加する。追加後は build() を呼ぶこと。
    void add(const T& x) {
        _xs.push_back(x);
        _built = false;
    }

    /// 区間 [first, last) の値を候補に追加する。追加後は build() を呼ぶこと。
    template <class InputIt>
    void add(InputIt first, InputIt last) {
        _xs.insert(_xs.end(), first, last);
        _built = false;
    }

    /// 配列 a の値を候補に追加する。追加後は build() を呼ぶこと。
    void add(const std::vector<T>& a) { add(a.begin(), a.end()); }

    /// 追加済みの候補をソート・重複除去し、座標圧縮を構築する。
    void build() {
        std::sort(_xs.begin(), _xs.end());
        _xs.erase(std::unique(_xs.begin(), _xs.end()), _xs.end());
        _built = true;
    }

    /// 候補を配列 a で置き換え、座標圧縮を構築する。
    void build(const std::vector<T>& a) {
        _xs = a;
        build();
    }

    /// 値 x の圧縮後 ID を返す。x は存在すること。
    int operator()(const T& x) const {
        assert(_built);
        auto it = std::lower_bound(_xs.begin(), _xs.end(), x);
        assert(it != _xs.end() && *it == x);
        return int(std::distance(_xs.begin(), it));
    }

    /// 値 x の圧縮後 ID を返す。存在しない場合は挿入位置。
    int get_id(const T& x) const {
        assert(_built);
        return int(std::distance(
            _xs.begin(), std::lower_bound(_xs.begin(), _xs.end(), x)));
    }

    /// 圧縮後ID id に対応する元の値を返す。
    const T& uncompress(int id) const {
        assert(_built);
        assert(0 <= id && id < int(_xs.size()));
        return _xs[id];
    }

    /// 圧縮後ID id に対応する元の値を返す。
    const T& operator[](int id) const { return uncompress(id); }

    /// 異なる値の個数を返す。build() 後に呼ぶこと。
    int size() const {
        assert(_built);
        return int(_xs.size());
    }

  private:
    std::vector<T> _xs;
    bool _built = false;
};

#endif  // LIBRARY_MISC_COMPRESS_HPP
