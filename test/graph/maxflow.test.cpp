#define PROBLEM "https://onlinejudge.u-aizu.ac.jp/courses/library/5/GRL/all/GRL_6_A"

#include <iostream>
#include "../../src/graph/maxflow.hpp"

using namespace std;

int main() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);

    int n, m;
    if (!(cin >> n >> m)) return 0;

    mf_graph<int> graph(n);

    for (int i = 0; i < m; i++) {
        int u, v, c;
        cin >> u >> v >> c;
        graph.add_edge(u, v, c);
    }

    // 頂点 0 から 頂点 n-1 への最大流
    cout << graph.flow(0, n - 1) << "\n";

    return 0;
}