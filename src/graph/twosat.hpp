#pragma once

#include <vector>
#include <cassert>
#include "scc.hpp"

struct twosat {
  public:
    twosat() : _n(0), _has_values(false) {}
    explicit twosat(int n) : _n(n), _has_values(false), _answer(n), g(2 * n) {}

    // (x_i == f) v (x_j == g) というクローズを追加する
    void add_clause(int i, bool f, int j, bool g) {
        assert(0 <= i && i < _n);
        assert(0 <= j && j < _n);
        // ¬x_i => x_j  (f が偽なら真、真なら偽に変換して命題の辺を張る)
        g.add_edge(2 * i + (f ? 0 : 1), 2 * j + (g ? 1 : 0));
        // ¬x_j => x_i
        g.add_edge(2 * j + (g ? 0 : 1), 2 * i + (f ? 1 : 0));
    }

    bool satisfiable() {
        auto scc = g.scc();
        std::vector<int> id(2 * _n);
        for (int i = 0; i < int(scc.size()); i++) {
            for (int v : scc[i]) {
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
    scc_graph g;
};