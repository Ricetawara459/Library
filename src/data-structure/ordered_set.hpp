#ifndef LIBRARY_DATA_STRUCTURE_ORDERED_SET_HPP
#define LIBRARY_DATA_STRUCTURE_ORDERED_SET_HPP

#pragma once

#include <cassert>
#include <cstdint>
#include <functional>
#include <vector>

/// 重複を持たない順序統計付き集合。各操作は期待 O(log n)。
/// Compare の順序で kth_smallest と lower_bound/upper_bound の位置を数える。
template <class T, class Compare = std::less<T>>
struct ordered_set {
  private:
    struct node {
        T key;
        std::uint64_t priority;
        int left = -1, right = -1, size = 1;
        node(const T& x, std::uint64_t p) : key(x), priority(p) {}
    };
    std::vector<node> nodes;
    std::vector<int> free_nodes;
    int root = -1;
    Compare cmp;
    std::uint64_t rng = 0x9e3779b97f4a7c15ULL;

    int size_of(int v) const { return v == -1 ? 0 : nodes[v].size; }
    void update(int v) { nodes[v].size = 1 + size_of(nodes[v].left) + size_of(nodes[v].right); }
    std::uint64_t next_random() {
        rng ^= rng << 7;
        rng ^= rng >> 9;
        return rng;
    }
    int make_node(const T& x) {
        auto p = next_random();
        if (free_nodes.empty()) {
            nodes.emplace_back(x, p);
            return int(nodes.size()) - 1;
        }
        int v = free_nodes.back();
        free_nodes.pop_back();
        nodes[v] = node(x, p);
        return v;
    }
    int rotate_right(int v) {
        int u = nodes[v].left;
        nodes[v].left = nodes[u].right;
        nodes[u].right = v;
        update(v); update(u);
        return u;
    }
    int rotate_left(int v) {
        int u = nodes[v].right;
        nodes[v].right = nodes[u].left;
        nodes[u].left = v;
        update(v); update(u);
        return u;
    }
    int insert_rec(int v, const T& x, bool& inserted) {
        if (v == -1) { inserted = true; return make_node(x); }
        if (cmp(x, nodes[v].key)) {
            int child = insert_rec(nodes[v].left, x, inserted);
            nodes[v].left = child;
            if (nodes[child].priority > nodes[v].priority) v = rotate_right(v);
        } else if (cmp(nodes[v].key, x)) {
            int child = insert_rec(nodes[v].right, x, inserted);
            nodes[v].right = child;
            if (nodes[child].priority > nodes[v].priority) v = rotate_left(v);
        }
        update(v);
        return v;
    }
    int merge(int a, int b) {
        if (a == -1) return b;
        if (b == -1) return a;
        if (nodes[a].priority > nodes[b].priority) {
            nodes[a].right = merge(nodes[a].right, b);
            update(a);
            return a;
        }
        nodes[b].left = merge(a, nodes[b].left);
        update(b);
        return b;
    }
    int erase_rec(int v, const T& x, bool& erased) {
        if (v == -1) return -1;
        if (cmp(x, nodes[v].key)) nodes[v].left = erase_rec(nodes[v].left, x, erased);
        else if (cmp(nodes[v].key, x)) nodes[v].right = erase_rec(nodes[v].right, x, erased);
        else {
            int result = merge(nodes[v].left, nodes[v].right);
            free_nodes.push_back(v);
            erased = true;
            return result;
        }
        update(v);
        return v;
    }

  public:
    explicit ordered_set(Compare compare = Compare()) : cmp(compare) {}
    /// 要素数を返す。
    int size() const { return size_of(root); }
    /// 空なら true。
    bool empty() const { return root == -1; }
    /// x を挿入する。新規に挿入した場合 true。
    bool insert(const T& x) { bool ok = false; root = insert_rec(root, x, ok); return ok; }
    /// x を削除する。存在した場合 true。
    bool erase(const T& x) { bool ok = false; root = erase_rec(root, x, ok); return ok; }
    /// x が存在するか返す。
    bool contains(const T& x) const {
        int v = root;
        while (v != -1) {
            if (cmp(x, nodes[v].key)) v = nodes[v].left;
            else if (cmp(nodes[v].key, x)) v = nodes[v].right;
            else return true;
        }
        return false;
    }
    /// 0-indexed で k 番目の要素を返す。0 <= k < size()。
    T kth_smallest(int k) const {
        assert(0 <= k && k < size());
        int v = root;
        while (true) {
            int left_size = size_of(nodes[v].left);
            if (k < left_size) v = nodes[v].left;
            else if (k == left_size) return nodes[v].key;
            else { k -= left_size + 1; v = nodes[v].right; }
        }
    }
    /// x 未満の要素数、つまり x の挿入位置を返す。値ではなく0-indexed位置。
    int lower_bound(const T& x) const {
        int v = root, result = 0;
        while (v != -1) {
            if (cmp(nodes[v].key, x)) {
                result += size_of(nodes[v].left) + 1;
                v = nodes[v].right;
            } else v = nodes[v].left;
        }
        return result;
    }
    /// x 以下の要素数、つまり x より大きい要素の先頭位置を返す。
    int upper_bound(const T& x) const {
        int v = root, result = 0;
        while (v != -1) {
            if (!cmp(x, nodes[v].key)) {
                result += size_of(nodes[v].left) + 1;
                v = nodes[v].right;
            } else v = nodes[v].left;
        }
        return result;
    }
};

#endif  // LIBRARY_DATA_STRUCTURE_ORDERED_SET_HPP
