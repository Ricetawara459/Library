#ifndef LIBRARY_GRAPH_MINCOSTFLOW_HPP
#define LIBRARY_GRAPH_MINCOSTFLOW_HPP

#line 1 "src/graph/mincostflow.hpp"

#include <vector>
#include <queue>
#include <algorithm>
#include <cassert>
#include <limits>
#include <utility>

template <class Cap, class Cost>
struct mcf_graph {
  public:
    mcf_graph() : _n(0) {}
    explicit mcf_graph(int n) : _n(n), g(n) {}

    /// from から to へ容量 cap・費用 cost の有向辺を追加し、辺番号を返す。
    int add_edge(int from, int to, Cap cap, Cost cost) {
        assert(0 <= from && from < _n);
        assert(0 <= to && to < _n);
        assert(0 <= cap);
        assert(0 <= cost);
        int m = int(pos.size());
        pos.push_back({from, int(g[from].size())});
        int from_id = int(g[from].size());
        int to_id = int(g[to].size());
        if (from == to) to_id++;
        g[from].push_back(_edge{to, to_id, cap, cost});
        g[to].push_back(_edge{from, from_id, 0, -cost});
        return m;
    }

    struct edge {
        int from;
        int to;
        Cap cap;
        Cap flow;
        Cost cost;
    };

    /// i 番目に追加した辺の現在の容量・流量・費用を返す。
    edge get_edge(int i) const {
        int m = int(pos.size());
        assert(0 <= i && i < m);
        auto _e = g[pos[i].first][pos[i].second];
        auto _re = g[_e.to][_e.rev];
        return edge{pos[i].first, _e.to, _e.cap + _re.cap, _re.cap, _e.cost};
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

    /// s から t へ流せるだけ流し、{流量, 最小費用} を返す。
    std::pair<Cap, Cost> flow(int s, int t) {
        return flow(s, t, std::numeric_limits<Cap>::max());
    }
    
    /// 流量上限 flow_limit まで流し、{流量, 最小費用} を返す。
    std::pair<Cap, Cost> flow(int s, int t, Cap flow_limit) {
        return slope(s, t, flow_limit).back();
    }

    /// 流量を増やしたときの {流量, 最小費用} の折れ線を返す。
    std::vector<std::pair<Cap, Cost>> slope(int s, int t) {
        return slope(s, t, std::numeric_limits<Cap>::max());
    }

    /// 流量上限つきで、{流量, 最小費用} の折れ線を返す。
    std::vector<std::pair<Cap, Cost>> slope(int s, int t, Cap flow_limit) {
        assert(0 <= s && s < _n);
        assert(0 <= t && t < _n);
        assert(s != t);

        std::vector<Cost> dual(_n, 0), dist(_n);
        std::vector<int> pv(_n), pe(_n);
        std::vector<bool> vis(_n);

        auto dual_ref = [&]() {
            std::fill(dist.begin(), dist.end(), std::numeric_limits<Cost>::max());
            std::fill(vis.begin(), vis.end(), false);
            dist[s] = 0;
            
            using P = std::pair<Cost, int>;
            std::priority_queue<P, std::vector<P>, std::greater<P>> que;
            que.push({0, s});
            
            while (!que.empty()) {
                int v = que.top().second;
                que.pop();
                if (vis[v]) continue;
                vis[v] = true;
                if (v == t) break;
                
                for (int i = 0; i < int(g[v].size()); i++) {
                    auto e = g[v][i];
                    if (vis[e.to] || e.cap == 0) continue;
                    Cost cost = e.cost - dual[e.to] + dual[v];
                    if (dist[e.to] - dist[v] > cost) {
                        dist[e.to] = dist[v] + cost;
                        pv[e.to] = v;
                        pe[e.to] = i;
                        que.push({dist[e.to], e.to});
                    }
                }
            }
            if (!vis[t]) return false;

            for (int v = 0; v < _n; v++) {
                if (!vis[v]) continue;
                dual[v] -= dist[t] - dist[v];
            }
            return true;
        };

        Cap flow = 0;
        Cost cost = 0;
        Cost prev_cost_per_flow = -1;
        std::vector<std::pair<Cap, Cost>> result;
        result.push_back({flow, cost});

        while (flow < flow_limit) {
            if (!dual_ref()) break;
            Cap c = flow_limit - flow;
            for (int v = t; v != s; v = pv[v]) {
                c = std::min(c, g[pv[v]][pe[v]].cap);
            }
            for (int v = t; v != s; v = pv[v]) {
                auto& e = g[pv[v]][pe[v]];
                e.cap -= c;
                g[v][e.rev].cap += c;
            }
            Cost d = -dual[s];
            flow += c;
            cost += c * d;
            if (prev_cost_per_flow == d) {
                result.pop_back();
            }
            result.push_back({flow, cost});
            prev_cost_per_flow = d;
        }
        return result;
    }

  private:
    int _n;
    struct _edge {
        int to;
        int rev;
        Cap cap;
        Cost cost;
    };
    std::vector<std::pair<int, int>> pos;
    std::vector<std::vector<_edge>> g;
};

#endif  // LIBRARY_GRAPH_MINCOSTFLOW_HPP
