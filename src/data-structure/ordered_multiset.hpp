#ifndef LIBRARY_DATA_STRUCTURE_ORDERED_MULTISET_HPP
#define LIBRARY_DATA_STRUCTURE_ORDERED_MULTISET_HPP

#pragma once

#include <cassert>
#include <cstdint>
#include <functional>
#include <limits>
#include <utility>
#include "ordered_set.hpp"

/// 重複を含む順序統計付き多重集合。各操作は期待 O(log n)。
/// Compare の順序で位置を数える。pair も辞書順で扱える。
template <class T, class Compare = std::less<T>>
struct ordered_multiset {
  private:
    using key = std::pair<T, std::uint64_t>;
    struct key_compare {
        Compare cmp;
        bool operator()(const key& a, const key& b) const {
            if (cmp(a.first, b.first)) return true;
            if (cmp(b.first, a.first)) return false;
            return a.second < b.second;
        }
    };
    Compare cmp;
    ordered_set<key, key_compare> tree;
    std::uint64_t next_id = 0;

  public:
    explicit ordered_multiset(Compare compare = Compare()) : cmp(compare), tree(key_compare{compare}) {}
    /// 要素数を返す。重複を含む。
    int size() const { return tree.size(); }
    /// 空なら true。
    bool empty() const { return tree.empty(); }
    /// x を1個追加する。
    void insert(const T& x) {
        assert(next_id < std::numeric_limits<std::uint64_t>::max());
        tree.insert({x, next_id++});
    }
    /// x を1個削除する。存在した場合 true。
    bool erase(const T& x) {
        int i = lower_bound(x);
        if (i == size()) return false;
        auto found = tree.kth_smallest(i);
        if (cmp(x, found.first) || cmp(found.first, x)) return false;
        return tree.erase(found);
    }
    /// x の個数を返す。
    int count(const T& x) const { return upper_bound(x) - lower_bound(x); }
    /// x が存在するか返す。
    bool contains(const T& x) const { return count(x) != 0; }
    /// 重複を含めて0-indexed k番目の値を返す。0 <= k < size()。
    T kth_smallest(int k) const { return tree.kth_smallest(k).first; }
    /// x より前にある要素数を返す。0-indexedの挿入位置。
    int lower_bound(const T& x) const { return tree.lower_bound({x, 0}); }
    /// x 以下の要素数を返す。重複を含む。
    int upper_bound(const T& x) const {
        return tree.lower_bound({x, std::numeric_limits<std::uint64_t>::max()});
    }
};

#endif  // LIBRARY_DATA_STRUCTURE_ORDERED_MULTISET_HPP
