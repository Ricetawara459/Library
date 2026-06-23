#define PROBLEM "https://judge.yosupo.jp/problem/scc"

#include <iostream>
#include <vector>
#include "../../src/graph/scc.hpp"

using namespace std;

int main() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);

    int n, m;
    if (!(cin >> n >> m)) return 0;

    scc_graph g(n);
    for (int i = 0; i < m; i++) {
        int u, v;
        cin >> u >> v;
        g.add_edge(u, v);
    }

    auto groups = g.scc();

    cout << groups.size() << "\n";
    for (const auto& group : groups) {
        cout << group.size();
        for (int v : group) {
            cout << " " << v;
        }
        cout << "\n";
    }

    return 0;
}