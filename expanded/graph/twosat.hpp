#line 2 "src/graph/twosat.hpp"

#include <vector>
#include <cassert>
#line 2 "src/graph/scc.hpp"

#line 4 "src/graph/scc.hpp"
#include <algorithm>
#include <utility>

struct scc_graph {
  public:
    scc_graph() : _n(0) {}
    explicit scc_graph(int n) : _n(n) {}

    void add_edge(int from, int to) {
        edges.push_back({from, to});
    }

    std::vector<std::vector<int>> scc() {
        std::vector<int> start(_n + 1, 0);
        for (const auto& e : edges) {
            start[e.from + 1]++;
        }
        for (int i = 0; i < _n; i++) {
            start[i + 1] += start[i];
        }
        auto counter = start;
        std::vector<int> elist(edges.size());
        for (const auto& e : edges) {
            elist[counter[e.from]++] = e.to;
        }

        int now_ord = 0, group_num = 0;
        std::vector<int> low(_n), ord(_n, -1), ids(_n);
        std::vector<int> visited;
        visited.reserve(_n);

        // 非再帰 DFS 用のスタック: {頂点, 次に探索する辺のインデックス}
        std::vector<std::pair<int, int>> call_stack;
        call_stack.reserve(_n);

        for (int i = 0; i < _n; i++) {
            if (ord[i] != -1) continue;
            
            call_stack.push_back({i, start[i]});
            ord[i] = low[i] = now_ord++;
            visited.push_back(i);

            while (!call_stack.empty()) {
                int u = call_stack.back().first;
                int& edge_idx = call_stack.back().second;

                if (edge_idx < start[u + 1]) {
                    int to = elist[edge_idx];
                    edge_idx++;
                    if (ord[to] == -1) {
                        ord[to] = low[to] = now_ord++;
                        visited.push_back(to);
                        call_stack.push_back({to, start[to]});
                    } else {
                        low[u] = std::min(low[u], ord[to]);
                    }
                } else {
                    call_stack.pop_back();
                    if (!call_stack.empty()) {
                        int parent = call_stack.back().first;
                        low[parent] = std::min(low[parent], low[u]);
                    }
                    if (low[u] == ord[u]) {
                        while (true) {
                            int v = visited.back();
                            visited.pop_back();
                            ord[v] = _n;
                            ids[v] = group_num;
                            if (v == u) break;
                        }
                        group_num++;
                    }
                }
            }
        }

        for (int i = 0; i < _n; i++) {
            ids[i] = group_num - 1 - ids[i];
        }

        std::vector<std::vector<int>> groups(group_num);
        std::vector<int> counts(group_num);
        for (int i = 0; i < _n; i++) {
            counts[ids[i]]++;
        }
        for (int i = 0; i < group_num; i++) {
            groups[i].reserve(counts[i]);
        }
        for (int i = 0; i < _n; i++) {
            groups[ids[i]].push_back(i);
        }
        return groups;
    }

  private:
    int _n;
    struct edge {
        int from;
        int to;
    };
    std::vector<edge> edges;
};
#line 6 "src/graph/twosat.hpp"

struct twosat {
  public:
    twosat() : _n(0), _has_values(false) {}
    explicit twosat(int n) : _n(n), _has_values(false), _answer(n), scc(2 * n) {}

    // (x_i == f) v (x_j == g) というクローズを追加する
    void add_clause(int i, bool f, int j, bool g) {
        assert(0 <= i && i < _n);
        assert(0 <= j && j < _n);
        // ¬x_i => x_j
        scc.add_edge(2 * i + (f ? 0 : 1), 2 * j + (g ? 1 : 0));
        // ¬x_j => x_i
        scc.add_edge(2 * j + (g ? 0 : 1), 2 * i + (f ? 1 : 0));
    }

    bool satisfiable() {
        auto scc_res = scc.scc();
        std::vector<int> id(2 * _n);
        for (int i = 0; i < int(scc_res.size()); i++) {
            for (int v : scc_res[i]) {
                id[v] = i;
            }
        }
        for (int i = 0; i < _n; i++) {
            if (id[2 * i] == id[2 * i + 1]) {
                return false;
            }
            _answer[i] = id[2 * i] > id[2 * i + 1];
        }
        _has_values = true;
        return true;
    }

    std::vector<bool> answer() {
        assert(_has_values);
        return _answer;
    }

  private:
    int _n;
    bool _has_values;
    std::vector<bool> _answer;
    scc_graph scc; // 変数名衝突を避けるため g から scc に変更
};
