#define PROBLEM "https://onlinejudge.u-aizu.ac.jp/courses/library/5/GRL/all/GRL_6_B"

#include <iostream>
#include "../../src/graph/mincostflow.hpp"

using namespace std;

int main() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);

    int n, m;
    long long f;
    if (!(cin >> n >> m >> f)) return 0;

    mcf_graph<long long, long long> graph(n);

    for (int i = 0; i < m; i++) {
        int u, v;
        long long cap, cost;
        cin >> u >> v >> cap >> cost;
        graph.add_edge(u, v, cap, cost);
    }

    // 頂点 0 から 頂点 n-1 へ流量 f の最小費用流を流す
    auto res = graph.flow(0, n - 1, f);

    // 要求された流量 f を流し切れなかった場合は -1
    if (res.first < f) {
        cout << -1 << "\n";
    } else {
        cout << res.second << "\n";
    }

    return 0;
}