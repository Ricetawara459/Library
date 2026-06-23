#pragma once

#include <vector>
#include <cassert>
#include "scc.hpp"

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
            // 【修正点】id[偽] < id[真] のときに true と判定する
            _answer[i] = id[2 * i] < id[2 * i + 1];
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
    scc_graph scc;
};