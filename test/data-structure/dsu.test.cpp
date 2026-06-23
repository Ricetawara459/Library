#define PROBLEM "https://judge.yosupo.jp/problem/unionfind"

#include <iostream>
#include "../../src/data-structure/dsu.hpp"

using namespace std;

int main() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);

    int n, q;
    if (!(cin >> n >> q)) return 0;

    dsu d(n);

    for (int i = 0; i < q; i++) {
        int t, u, v;
        cin >> t >> u >> v;
        if (t == 0) {
            d.merge(u, v);
        } else {
            cout << (d.same(u, v) ? 1 : 0) << "\n";
        }
    }

    return 0;
}