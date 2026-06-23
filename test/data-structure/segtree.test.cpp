#define PROBLEM "https://judge.yosupo.jp/problem/point_add_range_sum"

#include <iostream>
#include <vector>
#include "../../src/data-structure/segment_tree.hpp" // 作成したセグ木をインクルード

using namespace std;

// セグ木に持たせる演算と単位元
long long op(long long a, long long b) { return a + b; }
long long e() { return 0LL; }

int main() {
    // 高速化
    cin.tie(nullptr);
    ios::sync_with_stdio(false);

    int n, q;
    if (!(cin >> n >> q)) return 0;

    vector<long long> a(n);
    for (int i = 0; i < n; i++) {
        cin >> a[i];
    }

    // セグ木の初期化
    segtree<long long, op, e> seg(a);

    for (int i = 0; i < q; i++) {
        int type;
        cin >> type;
        if (type == 0) {
            int p;
            long long x;
            cin >> p >> x;
            // 点加算
            seg.set(p, seg.get(p) + x);
        } else {
            int l, r;
            cin >> l >> r;
            // 区間和取得
            cout << seg.prod(l, r) << "\n";
        }
    }

    return 0;
}