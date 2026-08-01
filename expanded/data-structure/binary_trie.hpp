#ifndef LIBRARY_DATA_STRUCTURE_BINARY_TRIE_HPP
#define LIBRARY_DATA_STRUCTURE_BINARY_TRIE_HPP

#line 1 "src/data-structure/binary_trie.hpp"


#include <array>
#include <cassert>
#include <limits>
#include <type_traits>
#include <vector>

/// 0 <= x < 2^(max_bit+1) の整数を管理するBinary Trie。重複する値も保持する。
/// 例: 0 <= x < 2^30 なら binary_trie<int, 29>。
template <class T = unsigned int,
          int max_bit = std::numeric_limits<std::make_unsigned_t<T>>::digits - 1,
          class Sum = long long>
struct binary_trie {
    static_assert(std::is_integral_v<T>, "T must be an integral type");
    using unsigned_type = std::make_unsigned_t<T>;
    using sum_type = Sum;
    static_assert(std::is_arithmetic_v<sum_type>);
    static_assert(0 <= max_bit &&
                  max_bit < std::numeric_limits<unsigned_type>::digits);

  private:
    struct node {
        int next[2] = {-1, -1};
        int count = 0;
        sum_type alternating_sum = 0;
    };

  public:
    /// 空のBinary Trieを作る。
    binary_trie() : _nodes(1) {}

    /// 頂点配列の容量を node_capacity まで予約する。
    void reserve(int node_capacity) {
        assert(node_capacity >= 0);
        _nodes.reserve(node_capacity);
    }

    /// 値 x を 1 個追加する。計算量 O(max_bit)。
    void insert(T x) {
        unsigned_type ux = to_unsigned(x);
        std::array<int, max_bit + 2> path{};
        int v = 0;
        path[0] = v;
        for (int bit = max_bit; bit >= 0; bit--) {
            int b = int(ux >> bit & 1);
            if (_nodes[v].next[b] == -1) {
                _nodes[v].next[b] = int(_nodes.size());
                _nodes.emplace_back();
            }
            v = _nodes[v].next[b];
            path[max_bit - bit + 1] = v;
        }
        for (int u : path) _nodes[u].count++;
        update_leaf(path.back(), x);
        for (int depth = max_bit; depth >= 0; depth--) pull(path[depth]);
        _total_sum += static_cast<sum_type>(x);
    }

    /// 値 x を 1 個削除する。x が存在すること。計算量 O(max_bit)。
    void erase(T x) {
        assert(count(x) > 0);
        unsigned_type ux = to_unsigned(x);
        std::array<int, max_bit + 2> path{};
        int v = 0;
        path[0] = v;
        for (int bit = max_bit; bit >= 0; bit--) {
            int b = int(ux >> bit & 1);
            v = _nodes[v].next[b];
            path[max_bit - bit + 1] = v;
        }
        for (int u : path) _nodes[u].count--;
        update_leaf(path.back(), x);
        for (int depth = max_bit; depth >= 0; depth--) pull(path[depth]);
        _total_sum -= static_cast<sum_type>(x);
    }

    /// 値 x の登録個数を返す。計算量 O(max_bit)。
    int count(T x) const {
        unsigned_type ux = to_unsigned(x);
        int v = 0;
        for (int bit = max_bit; bit >= 0; bit--) {
            int b = int(ux >> bit & 1);
            v = _nodes[v].next[b];
            if (v == -1 || _nodes[v].count == 0) return 0;
        }
        return _nodes[v].count;
    }

    /// 値 x が 1 個以上登録されていれば true。計算量 O(max_bit)。
    bool contains(T x) const { return count(x) > 0; }

    /// x とのxorが最小になる登録要素を返す。空でないこと。計算量 O(max_bit)。
    T min_xor_element(T x) const {
        return xor_extreme_element(x, false);
    }

    /// x とのxorが最大になる登録要素を返す。空でないこと。計算量 O(max_bit)。
    T max_xor_element(T x) const {
        return xor_extreme_element(x, true);
    }

    /// 登録要素 y に対する x xor y の最小値を返す。空でないこと。
    T min_xor(T x) const {
        return T(to_unsigned(x) ^ to_unsigned(min_xor_element(x)));
    }

    /// 登録要素 y に対する x xor y の最大値を返す。空でないこと。
    T max_xor(T x) const {
        return T(to_unsigned(x) ^ to_unsigned(max_xor_element(x)));
    }

    /// 0-indexedで k 番目に小さい登録要素を返す。重複を含む。
    T kth_smallest(int k) const {
        assert(0 <= k && k < size());
        int v = 0;
        unsigned_type value = 0;
        for (int bit = max_bit; bit >= 0; bit--) {
            int left = child_count(v, 0);
            int b = 0;
            if (k >= left) {
                k -= left;
                b = 1;
                value |= unsigned_type(1) << bit;
            }
            v = _nodes[v].next[b];
        }
        return T(value);
    }

    /// xとのxorが0-indexedでk番目に小さくなる登録要素を返す。重複を含む。
    T kth_xor_element(T x, int k) const {
        assert(0 <= k && k < size());
        unsigned_type ux = to_unsigned(x);
        unsigned_type value = 0;
        int v = 0;
        for (int bit = max_bit; bit >= 0; bit--) {
            int xb = int(ux >> bit & 1);
            int same_count = child_count(v, xb);
            int b = xb;
            if (k >= same_count) {
                k -= same_count;
                b ^= 1;
            }
            if (b == 1) value |= unsigned_type(1) << bit;
            v = _nodes[v].next[b];
        }
        return T(value);
    }

    /// (x xor y) < upper を満たす登録要素 y の個数を返す。重複を含む。
    int count_xor_less(T x, T upper) const {
        unsigned_type ux = to_unsigned(x);
        unsigned_type limit = to_unsigned(upper);
        int v = 0;
        int res = 0;
        for (int bit = max_bit; bit >= 0 && v != -1; bit--) {
            int xb = int(ux >> bit & 1);
            int limit_bit = int(limit >> bit & 1);
            if (limit_bit == 1) {
                res += child_count(v, xb);
                v = _nodes[v].next[xb ^ 1];
            } else {
                v = _nodes[v].next[xb];
            }
        }
        return res;
    }

    /// x 未満の登録要素数を返す。重複を含む。計算量 O(max_bit)。
    int count_less(T x) const {
        unsigned_type ux = to_unsigned(x);
        int v = 0;
        int res = 0;
        for (int bit = max_bit; bit >= 0 && v != -1; bit--) {
            int b = int(ux >> bit & 1);
            if (b == 1) res += child_count(v, 0);
            v = _nodes[v].next[b];
        }
        return res;
    }

    /// 登録要素数を返す。重複を含む。
    int size() const { return _nodes[0].count; }

    /// 登録要素が 1 個もなければ true。
    bool empty() const { return size() == 0; }

    /// 昇順に並べたときの1,3,5,...番目（1-indexed）の要素の和を返す。
    sum_type odd_position_sum() const {
        return (_total_sum + _nodes[0].alternating_sum) / 2;
    }

    /// 昇順に並べたときの2,4,6,...番目（1-indexed）の要素の和を返す。
    sum_type even_position_sum() const {
        return (_total_sum - _nodes[0].alternating_sum) / 2;
    }

    /// Binary Trie が確保している頂点数を返す。
    int node_count() const { return int(_nodes.size()); }

  private:
    static unsigned_type to_unsigned(T x) {
        if constexpr (std::is_signed_v<T>) assert(x >= 0);
        unsigned_type ux = static_cast<unsigned_type>(x);
        if constexpr (
            max_bit + 1 < std::numeric_limits<unsigned_type>::digits) {
            assert((ux >> (max_bit + 1)) == 0);
        }
        return ux;
    }

    int child_count(int v, int b) const {
        int to = _nodes[v].next[b];
        return to == -1 ? 0 : _nodes[to].count;
    }

    void update_leaf(int v, T x) {
        sum_type value = static_cast<sum_type>(x);
        _nodes[v].alternating_sum = _nodes[v].count % 2 == 0 ? 0 : value;
    }

    void pull(int v) {
        int left = _nodes[v].next[0];
        int right = _nodes[v].next[1];
        int left_count = left == -1 ? 0 : _nodes[left].count;
        sum_type left_sum =
            left == -1 ? sum_type(0) : _nodes[left].alternating_sum;
        sum_type right_sum =
            right == -1 ? sum_type(0) : _nodes[right].alternating_sum;
        _nodes[v].alternating_sum =
            left_sum + (left_count % 2 == 0 ? right_sum : -right_sum);
    }

    T xor_extreme_element(T x, bool maximize) const {
        assert(!empty());
        unsigned_type ux = to_unsigned(x);
        unsigned_type value = 0;
        int v = 0;
        for (int bit = max_bit; bit >= 0; bit--) {
            int xb = int(ux >> bit & 1);
            int b = xb ^ int(maximize);
            if (child_count(v, b) == 0) b ^= 1;
            if (b == 1) value |= unsigned_type(1) << bit;
            v = _nodes[v].next[b];
        }
        return T(value);
    }

    std::vector<node> _nodes;
    sum_type _total_sum = 0;
};

#endif  // LIBRARY_DATA_STRUCTURE_BINARY_TRIE_HPP
