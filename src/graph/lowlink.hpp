#pragma once

#include <algorithm>
#include <cassert>
#include <vector>

// 無向グラフの LowLink。
// build() 後に以下が使える。
//   articulation: 関節点の一覧
//   bridges: 橋の一覧。各辺は (min(u, v), max(u, v)) の形で入る
// 多重辺は辺 ID で区別して扱う。
struct lowlink {
    struct edge {
        int to;
        int id;
    };

    lowlink() = default;
    explicit lowlink(int n) : _n(n), g(n), ord(n, -1), low(n, -1) {}

    /// 無向辺 u-v を追加する。多重辺も区別される。
    void add_edge(int u, int v) {
        assert(0 <= u && u < _n);
        assert(0 <= v && v < _n);
        int id = int(edges.size());
        edges.push_back({u, v});
        g[u].push_back({v, id});
        g[v].push_back({u, id});
    }

    /// ord/low を計算し、関節点 articulation と橋 bridges を作る。
    void build() {
        articulation.clear();
        bridges.clear();
        std::fill(ord.begin(), ord.end(), -1);
        std::fill(low.begin(), low.end(), -1);
        int timer = 0;
        for (int v = 0; v < _n; v++) {
            if (ord[v] == -1) dfs(v, -1, timer);
        }
        std::sort(articulation.begin(), articulation.end());
        articulation.erase(std::unique(articulation.begin(), articulation.end()), articulation.end());
        std::sort(bridges.begin(), bridges.end());
    }

    /// id 番目に追加した辺が橋なら true。
    bool is_bridge_id(int id) const {
        assert(0 <= id && id < int(edges.size()));
        auto [u, v] = edges[id];
        if (ord[u] > ord[v]) std::swap(u, v);
        return ord[u] < low[v];
    }

    /// 頂点数を返す。
    int size() const { return _n; }

    int _n = 0;
    std::vector<std::vector<edge>> g;
    std::vector<std::pair<int, int>> edges;
    std::vector<int> ord, low;
    std::vector<int> articulation;
    std::vector<std::pair<int, int>> bridges;

  private:
    void dfs(int v, int parent_edge, int& timer) {
        ord[v] = low[v] = timer++;
        int child_count = 0;
        bool is_articulation = false;

        for (auto e : g[v]) {
            if (e.id == parent_edge) continue;
            int to = e.to;
            if (ord[to] == -1) {
                child_count++;
                dfs(to, e.id, timer);
                low[v] = std::min(low[v], low[to]);
                if (parent_edge != -1 && ord[v] <= low[to]) is_articulation = true;
                if (ord[v] < low[to]) {
                    bridges.push_back(std::minmax(v, to));
                }
            } else {
                low[v] = std::min(low[v], ord[to]);
            }
        }

        if (parent_edge == -1 && child_count >= 2) is_articulation = true;
        if (is_articulation) articulation.push_back(v);
    }
};
