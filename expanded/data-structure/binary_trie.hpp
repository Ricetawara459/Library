// Source: https://github.com/Ricetawara459/Library
#ifndef LIBRARY_DATA_STRUCTURE_BINARY_TRIE_HPP
#define LIBRARY_DATA_STRUCTURE_BINARY_TRIE_HPP



#include <cassert>
#include <limits>
#include <type_traits>
#include <vector>

/// 非負整数の多重集合。max_bit は最上位ビットの番号。
/// 例: binary_trie<int, 29> は 0 <= x < 2^30 を扱う。
/// xor_all(x) は全要素を y -> y xor x に変更する。以降の操作は変更後の値を扱う。
template <class T = unsigned int, int max_bit = std::numeric_limits<T>::digits - 1>
struct binary_trie {
    static_assert(std::is_integral_v<T> && !std::is_same_v<T, bool>);
    static_assert(0 <= max_bit && max_bit < std::numeric_limits<T>::digits);
    using unsigned_type = std::make_unsigned_t<T>;

  private:
    struct node { int next[2] = {-1, -1}; int count = 0; };
    std::vector<node> nodes{1};
    unsigned_type lazy = 0;

    static unsigned_type checked(T x) {
        if constexpr (std::is_signed_v<T>) assert(x >= 0);
        unsigned_type u = static_cast<unsigned_type>(x);
        if constexpr (max_bit + 1 < std::numeric_limits<unsigned_type>::digits)
            assert((u >> (max_bit + 1)) == 0);
        return u;
    }
    int child_count(int v, int b) const {
        int to = nodes[v].next[b];
        return to == -1 ? 0 : nodes[to].count;
    }

  public:
    binary_trie() = default;
    /// 頂点配列の容量を予約する。
    void reserve(int node_capacity) { assert(node_capacity >= 0); nodes.reserve(node_capacity); }
    /// 登録要素数を返す。重複を含む。
    int size() const { return nodes[0].count; }
    /// 空なら true。
    bool empty() const { return size() == 0; }
    /// 確保した頂点数を返す。削除後の頂点も含む。
    int node_count() const { return int(nodes.size()); }

    /// 全要素に x を xor する。O(1)。
    void xor_all(T x) { lazy ^= checked(x); }
    /// 現在の値 x を1個追加する。O(max_bit+1)。
    void insert(T x) {
        unsigned_type key = checked(x) ^ lazy;
        int v = 0;
        nodes[v].count++;
        for (int bit = max_bit; bit >= 0; bit--) {
            int b = int(key >> bit & 1);
            if (nodes[v].next[b] == -1) {
                nodes[v].next[b] = int(nodes.size());
                nodes.emplace_back();
            }
            v = nodes[v].next[b];
            nodes[v].count++;
        }
    }

    /// 現在の値 x を1個削除する。x が存在すること。O(max_bit+1)。
    void erase(T x) {
        assert(count(x) > 0);
        unsigned_type key = checked(x) ^ lazy;
        int v = 0;
        nodes[v].count--;
        for (int bit = max_bit; bit >= 0; bit--) {
            v = nodes[v].next[int(key >> bit & 1)];
            nodes[v].count--;
        }
    }

    /// 現在の値 x の個数を返す。O(max_bit+1)。
    int count(T x) const {
        unsigned_type key = checked(x) ^ lazy;
        int v = 0;
        for (int bit = max_bit; bit >= 0; bit--) {
            v = nodes[v].next[int(key >> bit & 1)];
            if (v == -1 || nodes[v].count == 0) return 0;
        }
        return nodes[v].count;
    }
    /// 現在の値 x が存在するか返す。
    bool contains(T x) const { return count(x) > 0; }

    /// 現在の要素のうち、0-indexed k 番目に小さい値を返す。
    T kth_smallest(int k) const {
        assert(0 <= k && k < size());
        unsigned_type key = 0;
        int v = 0;
        for (int bit = max_bit; bit >= 0; bit--) {
            int b = int(lazy >> bit & 1);
            int cnt = child_count(v, b);
            if (k >= cnt) { k -= cnt; b ^= 1; }
            if (b) key |= unsigned_type(1) << bit;
            v = nodes[v].next[b];
        }
        return T(key ^ lazy);
    }
    /// 最小要素を返す。空でないこと。
    T min_element() const { return kth_smallest(0); }
    /// 最大要素を返す。空でないこと。
    T max_element() const { return kth_smallest(size() - 1); }
    /// 現在の値が upper 未満の要素数を返す。
    /// upper が値域以上なら全要素、upper <= 0 なら 0。O(max_bit+1)。
    int count_less(T upper) const {
        if constexpr (std::is_signed_v<T>) if (upper <= 0) return 0;
        unsigned_type limit = static_cast<unsigned_type>(upper);
        if constexpr (max_bit + 1 < std::numeric_limits<unsigned_type>::digits)
            if ((limit >> (max_bit + 1)) != 0) return size();
        int v = 0, res = 0;
        for (int bit = max_bit; bit >= 0 && v != -1; bit--) {
            int b = int(lazy >> bit & 1);
            if (limit >> bit & 1) {
                res += child_count(v, b);
                b ^= 1;
            }
            v = nodes[v].next[b];
        }
        return res;
    }
};

#endif  // LIBRARY_DATA_STRUCTURE_BINARY_TRIE_HPP
