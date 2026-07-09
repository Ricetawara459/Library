#pragma once

#include <algorithm>
#include <cassert>
#include <vector>

// Heavy-Light Decomposition。
// 各頂点 v は pos[v] によって [0, n) の位置へ写される。
// path_query(u, v, edge, f) はパスを O(log n) 個の半開区間 [l, r) に分けて f(l, r) を呼ぶ。
// edge == true のときは LCA の頂点を除外するので、辺の値を pos[child] に置いて扱える。
// path_query_ordered(u, v, edge, f) は u から v へ進む順に f(l, r, rev) を呼ぶ。
// rev == true の区間は [l, r) を右から左へ読む。非可換な演算や関数合成で使う。
struct hld {
    struct segment {
        int l, r;
        bool rev;
    };

  public:
    hld() = default;
    explicit hld(int n) : g(n), parent(n, -1), depth(n), heavy(n, -1), head(n), pos(n), rev(n), subtree_size(n), _n(n) {}

    void add_edge(int u, int v) {
        assert(0 <= u && u < _n);
        assert(0 <= v && v < _n);
        g[u].push_back(v);
        g[v].push_back(u);
    }

    void build(int root = 0) {
        assert(0 <= root && root < _n);
        int cur = 0;
        dfs_size(root, -1);
        dfs_hld(root, root, cur);
    }

    int lca(int u, int v) const {
        assert(0 <= u && u < _n);
        assert(0 <= v && v < _n);
        while (head[u] != head[v]) {
            if (depth[head[u]] > depth[head[v]]) {
                u = parent[head[u]];
            } else {
                v = parent[head[v]];
            }
        }
        return depth[u] < depth[v] ? u : v;
    }

    int distance(int u, int v) const {
        int w = lca(u, v);
        return depth[u] + depth[v] - 2 * depth[w];
    }

    template <class F>
    void path_query(int u, int v, bool edge, F f) const {
        assert(0 <= u && u < _n);
        assert(0 <= v && v < _n);
        while (head[u] != head[v]) {
            if (depth[head[u]] < depth[head[v]]) std::swap(u, v);
            f(pos[head[u]], pos[u] + 1);
            u = parent[head[u]];
        }
        if (depth[u] > depth[v]) std::swap(u, v);
        int l = pos[u] + (edge ? 1 : 0);
        int r = pos[v] + 1;
        if (l < r) f(l, r);
    }

    std::vector<std::pair<int, int>> path_segments(int u, int v, bool edge = false) const {
        std::vector<std::pair<int, int>> res;
        path_query(u, v, edge, [&](int l, int r) { res.push_back({l, r}); });
        return res;
    }

    template <class F>
    void path_query_ordered(int u, int v, bool edge, F f) const {
        for (auto seg : path_segments_ordered(u, v, edge)) {
            f(seg.l, seg.r, seg.rev);
        }
    }

    std::vector<segment> path_segments_ordered(int u, int v, bool edge = false) const {
        assert(0 <= u && u < _n);
        assert(0 <= v && v < _n);
        std::vector<segment> left, right;
        while (head[u] != head[v]) {
            if (depth[head[u]] >= depth[head[v]]) {
                left.push_back({pos[head[u]], pos[u] + 1, true});
                u = parent[head[u]];
            } else {
                right.push_back({pos[head[v]], pos[v] + 1, false});
                v = parent[head[v]];
            }
        }
        if (depth[u] >= depth[v]) {
            int l = pos[v] + (edge ? 1 : 0);
            int r = pos[u] + 1;
            if (l < r) left.push_back({l, r, true});
        } else {
            int l = pos[u] + (edge ? 1 : 0);
            int r = pos[v] + 1;
            if (l < r) left.push_back({l, r, false});
        }
        for (auto it = right.rbegin(); it != right.rend(); ++it) {
            left.push_back(*it);
        }
        return left;
    }

    std::pair<int, int> subtree_query(int v, bool edge = false) const {
        assert(0 <= v && v < _n);
        return {pos[v] + (edge ? 1 : 0), pos[v] + subtree_size[v]};
    }

    int operator[](int v) const {
        assert(0 <= v && v < _n);
        return pos[v];
    }

    int vertex_at(int p) const {
        assert(0 <= p && p < _n);
        return rev[p];
    }

    int size() const { return _n; }

    std::vector<std::vector<int>> g;
    std::vector<int> parent, depth, heavy, head, pos, rev, subtree_size;

  private:
    int _n = 0;

    void dfs_size(int v, int p) {
        parent[v] = p;
        subtree_size[v] = 1;
        int max_size = 0;
        for (int to : g[v]) {
            if (to == p) continue;
            depth[to] = depth[v] + 1;
            dfs_size(to, v);
            subtree_size[v] += subtree_size[to];
            if (subtree_size[to] > max_size) {
                max_size = subtree_size[to];
                heavy[v] = to;
            }
        }
    }

    void dfs_hld(int v, int h, int& cur) {
        head[v] = h;
        pos[v] = cur;
        rev[cur] = v;
        cur++;
        if (heavy[v] != -1) dfs_hld(heavy[v], h, cur);
        for (int to : g[v]) {
            if (to == parent[v] || to == heavy[v]) continue;
            dfs_hld(to, to, cur);
        }
    }
};
