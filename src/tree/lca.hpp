#ifndef LIBRARY_TREE_LCA_HPP
#define LIBRARY_TREE_LCA_HPP

#pragma once

#include <algorithm>
#include <cassert>
#include <queue>
#include <vector>

// ダブリングによる LCA。
// 使い方:
//   lca tree(n);
//   tree.add_edge(u, v, cost);  // cost は省略可
//   tree.build(root);
//   tree.query(u, v), tree.distance(u, v), tree.kth_ancestor(v, k)
struct lca {
  public:
    lca() = default;
    explicit lca(int n) : g(n), depth(n), dist(n), parent(n, -1), _n(n) {}

    /// 木の無向辺 u-v を追加する。w は距離用の重み。
    void add_edge(int u, int v, long long w = 1) {
        assert(0 <= u && u < _n);
        assert(0 <= v && v < _n);
        g[u].push_back({v, w});
        g[v].push_back({u, w});
    }

    /// root を根としてダブリング表を構築する。
    void build(int root = 0) {
        assert(0 <= root && root < _n);
        int lg = 1;
        while ((1 << lg) <= std::max(1, _n)) lg++;
        up.assign(lg, std::vector<int>(_n, -1));

        std::queue<int> que;
        que.push(root);
        parent[root] = -1;
        depth[root] = 0;
        dist[root] = 0;
        std::vector<int> seen(_n, 0);
        seen[root] = 1;

        while (!que.empty()) {
            int v = que.front();
            que.pop();
            up[0][v] = parent[v];
            for (auto [to, cost] : g[v]) {
                if (seen[to]) continue;
                seen[to] = 1;
                parent[to] = v;
                depth[to] = depth[v] + 1;
                dist[to] = dist[v] + cost;
                que.push(to);
            }
        }
        for (int k = 0; k + 1 < lg; k++) {
            for (int v = 0; v < _n; v++) {
                if (up[k][v] == -1) {
                    up[k + 1][v] = -1;
                } else {
                    up[k + 1][v] = up[k][up[k][v]];
                }
            }
        }
    }

    /// v から k 個上の祖先を返す。存在しなければ -1。
    int kth_ancestor(int v, long long k) const {
        assert(0 <= v && v < _n);
        for (int i = 0; i < int(up.size()) && v != -1; i++) {
            if ((k >> i) & 1) v = up[i][v];
        }
        return v;
    }

    /// u と v の LCA を返す。
    int query(int u, int v) const {
        assert(0 <= u && u < _n);
        assert(0 <= v && v < _n);
        if (depth[u] < depth[v]) std::swap(u, v);
        u = kth_ancestor(u, depth[u] - depth[v]);
        if (u == v) return u;
        for (int k = int(up.size()) - 1; k >= 0; k--) {
            if (up[k][u] != up[k][v]) {
                u = up[k][u];
                v = up[k][v];
            }
        }
        return parent[u];
    }

    int operator()(int u, int v) const {
        return query(u, v);
    }

    /// u-v パスの重み付き距離を返す。
    long long distance(int u, int v) const {
        int w = query(u, v);
        return dist[u] + dist[v] - 2 * dist[w];
    }

    /// from から to へ k 辺進んだ頂点を返す。範囲外なら -1。
    int jump(int from, int to, long long k) const {
        int w = query(from, to);
        long long d1 = depth[from] - depth[w];
        long long d2 = depth[to] - depth[w];
        if (k > d1 + d2) return -1;
        if (k <= d1) return kth_ancestor(from, k);
        return kth_ancestor(to, d1 + d2 - k);
    }

    /// 頂点数を返す。
    int size() const { return _n; }

    std::vector<std::vector<std::pair<int, long long>>> g;
    std::vector<int> depth;
    std::vector<long long> dist;
    std::vector<int> parent;

  private:
    int _n = 0;
    std::vector<std::vector<int>> up;
};

#endif  // LIBRARY_TREE_LCA_HPP
