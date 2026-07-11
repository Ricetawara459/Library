#ifndef LIBRARY_DATA_STRUCTURE_DSU_HPP
#define LIBRARY_DATA_STRUCTURE_DSU_HPP

#line 1 "src/data-structure/dsu.hpp"

#include <vector>
#include <numeric>
#include <algorithm>

struct dsu {
  public:
    dsu() : _n(0) {}
    /// n 頂点の Union-Find を作る。
    explicit dsu(int n) : _n(n), parent_or_size(n, -1) {}

    /// a と b を同じ集合にまとめ、代表元を返す。
    int merge(int a, int b) {
        int x = leader(a), y = leader(b);
        if (x == y) return x;
        if (-parent_or_size[x] < -parent_or_size[y]) std::swap(x, y);
        parent_or_size[x] += parent_or_size[y];
        parent_or_size[y] = x;
        return x;
    }

    /// a と b が同じ集合なら true。
    bool same(int a, int b) {
        return leader(a) == leader(b);
    }

    /// a が属する集合の代表元を返す。
    int leader(int a) {
        if (parent_or_size[a] < 0) return a;
        return parent_or_size[a] = leader(parent_or_size[a]);
    }

    /// a が属する集合のサイズを返す。
    int size(int a) {
        return -parent_or_size[leader(a)];
    }

    /// 全ての連結成分を頂点リストとして返す。
    std::vector<std::vector<int>> groups() {
        std::vector<int> leader_buf(_n), group_size(_n);
        for (int i = 0; i < _n; i++) {
            leader_buf[i] = leader(i);
            group_size[leader_buf[i]]++;
        }
        std::vector<std::vector<int>> result(_n);
        for (int i = 0; i < _n; i++) {
            result[i].reserve(group_size[i]);
        }
        for (int i = 0; i < _n; i++) {
            result[leader_buf[i]].push_back(i);
        }
        result.erase(
            std::remove_if(result.begin(), result.end(),
                           [](const std::vector<int>& v) { return v.empty(); }),
            result.end());
        return result;
    }

  private:
    int _n;
    // 負の値の場合：自身が根であり、その絶対値が連結成分のサイズ
    // 正の値の場合：親のインデックス
    std::vector<int> parent_or_size;
};

#endif  // LIBRARY_DATA_STRUCTURE_DSU_HPP
