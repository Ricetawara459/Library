#define PROBLEM "https://judge.yosupo.jp/problem/point_add_range_sum"

#include <iostream>
#include <vector>
#include "../../src/data-structure/fenwick_tree.hpp"

using namespace std;

int main() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);

    int n, q;
    if (!(cin >> n >> q)) return 0;

    fenwick_tree<long long> fw(n);
    for (int i = 0; i < n; i++) {
        long long a;
        cin >> a;
        fw.add(i, a);
    }

    for (int i = 0; i < q; i++) {
        int type;
        cin >> type;
        if (type == 0) {
            int p;
            long long x;
            cin >> p >> x;
            fw.add(p, x);
        } else {
            int l, r;
            cin >> l >> r;
            cout << fw.sum(l, r) << "\n";
        }
    }

    return 0;
}