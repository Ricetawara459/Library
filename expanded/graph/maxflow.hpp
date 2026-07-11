#ifndef LIBRARY_GRAPH_MAXFLOW_HPP
#define LIBRARY_GRAPH_MAXFLOW_HPP

#line 1 "src/graph/maxflow.hpp"

#include <vector>
#include <queue>
#include <algorithm>
#include <cassert>
#include <limits>

template <class Cap>
struct mf_graph {
  public:
    mf_graph() : _n(0) {}
    explicit mf_graph(int n) : _n(n), g(n) {}

    /// from から to へ容量 cap の有向辺を追加し、辺番号を返す。
    int add_edge(int from, int to, Cap cap) {
        assert(0 <= from && from < _n);
        assert(0 <= to && to < _n);
        assert(0 <= cap);
        int m = int(pos.size());
        pos.push_back({from, int(g[from].size())});
        int from_id = int(g[from].size());
        int to_id = int(g[to].size());
        if (from == to) to_id++;
        g[from].push_back(_edge{to, to_id, cap});
        g[to].push_back(_edge{from, from_id, 0});
        return m;
    }

    struct edge {
        int from;
        int to;
        Cap cap;
        Cap flow;
    };

    /// i 番目に追加した辺の現在の容量・流量を返す。
    edge get_edge(int i) const {
        int m = int(pos.size());
        assert(0 <= i && i < m);
        auto _e = g[pos[i].first][pos[i].second];
        auto _re = g[_e.to][_e.rev];
        return edge{pos[i].first, _e.to, _e.cap + _re.cap, _re.cap};
    }

    /// 追加した全ての辺の現在状態を返す。
    std::vector<edge> edges() const {
        int m = int(pos.size());
        std::vector<edge> result;
        for (int i = 0; i < m; i++) {
            result.push_back(get_edge(i));
        }
        return result;
    }

    /// i 番目の辺の容量と流量を直接変更する。
    void change_edge(int i, Cap new_cap, Cap new_flow) {
        int m = int(pos.size());
        assert(0 <= i && i < m);
        assert(0 <= new_flow && new_flow <= new_cap);
        auto& _e = g[pos[i].first][pos[i].second];
        auto& _re = g[_e.to][_e.rev];
        _e.cap = new_cap - new_flow;
        _re.cap = new_flow;
    }

    /// s から t への最大流を流せるだけ流す。
    Cap flow(int s, int t) {
        return flow(s, t, std::numeric_limits<Cap>::max());
    }
    /// 流量上限 flow_limit まで、s から t への最大流を流す。
    Cap flow(int s, int t, Cap flow_limit) {
        assert(0 <= s && s < _n);
        assert(0 <= t && t < _n);
        assert(s != t);

        std::vector<int> level(_n), iter(_n);
        
        auto bfs = [&]() {
            std::fill(level.begin(), level.end(), -1);
            level[s] = 0;
            std::queue<int> que;
            que.push(s);
            while (!que.empty()) {
                int v = que.front();
                que.pop();
                for (const auto& e : g[v]) {
                    if (e.cap == 0 || level[e.to] >= 0) continue;
                    level[e.to] = level[v] + 1;
                    if (e.to == t) return;
                    que.push(e.to);
                }
            }
        };

        auto dfs = [&](auto self, int v, Cap up) -> Cap {
            if (v == s) return up;
            Cap res = 0;
            int level_v = level[v];
            for (int& i = iter[v]; i < int(g[v].size()); i++) {
                _edge& e = g[v][i];
                if (level_v <= level[e.to] || g[e.to][e.rev].cap == 0) continue;
                Cap d = self(self, e.to, std::min(up - res, g[e.to][e.rev].cap));
                if (d <= 0) continue;
                g[e.to][e.rev].cap -= d;
                g[v][i].cap += d;
                res += d;
                if (res == up) return res;
            }
            level[v] = _n;
            return res;
        };

        Cap flow = 0;
        while (flow < flow_limit) {
            bfs();
            if (level[t] == -1) break;
            std::fill(iter.begin(), iter.end(), 0);
            while (flow < flow_limit) {
                Cap f = dfs(dfs, t, flow_limit - flow);
                if (!f) break;
                flow += f;
            }
        }
        return flow;
    }

    /// 最後の flow 後の残余グラフで、s から到達可能な頂点集合を返す。
    std::vector<bool> min_cut(int s) {
        std::vector<bool> visited(_n, false);
        std::queue<int> que;
        que.push(s);
        visited[s] = true;
        while (!que.empty()) {
            int p = que.front();
            que.pop();
            for (const auto& e : g[p]) {
                if (e.cap > 0 && !visited[e.to]) {
                    visited[e.to] = true;
                    que.push(e.to);
                }
            }
        }
        return visited;
    }

  private:
    int _n;
    struct _edge {
        int to;
        int rev;
        Cap cap;
    };
    std::pair<int, int> pos_t;
    std::vector<std::pair<int, int>> pos;
    std::vector<std::vector<_edge>> g;
};

#endif  // LIBRARY_GRAPH_MAXFLOW_HPP
