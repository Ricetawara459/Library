#define PROBLEM "https://judge.yosupo.jp/problem/two_sat"

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "../../src/graph/twosat.hpp"

using namespace std;

int main() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);

    string p, cnf;
    int n, m;
    if (!(cin >> p >> cnf >> n >> m)) return 0;

    twosat ts(n);

    for (int k = 0; k < m; k++) {
        int i, j, zero;
        cin >> i >> j >> zero;
        
        // 1-indexed で、負数は否定を表す
        bool f = (i > 0);
        bool g = (j > 0);
        
        int var_i = abs(i) - 1;
        int var_j = abs(j) - 1;
        
        ts.add_clause(var_i, f, var_j, g);
    }

    if (!ts.satisfiable()) {
        cout << "s UNSATISFIABLE\n";
    } else {
        cout << "s SATISFIABLE\n";
        cout << "v";
        auto ans = ts.answer();
        for (int i = 0; i < n; i++) {
            // 真なら i+1、偽なら -(i+1) を出力
            if (ans[i]) {
                cout << " " << (i + 1);
            } else {
                cout << " " << -(i + 1);
            }
        }
        cout << " 0\n";
    }

    return 0;
}