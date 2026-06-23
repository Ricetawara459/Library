#define PROBLEM "https://judge.yosupo.jp/problem/range_affine_range_sum"

#include <iostream>
#include <vector>
#include "../../src/data-structure/lazysegtree.hpp"
#include "../../src/math/modint.hpp"

using namespace std;

using mint = modint998244353;

// セグ木のモノイド（値、区間の長さ）
struct S {
    mint val;
    int size;
};

// 作用素（b, c） -> f(x) = b * x + c
struct F {
    mint b, c;
};

S op(S a, S b) {
    return S{a.val + b.val, a.size + b.size};
}
S e() {
    return S{0, 0};
}
S mapping(F f, S x) {
    return S{f.b * x.val + f.c * x.size, x.size};
}
F composition(F f, F g) {
    // f(g(x)) = f(b_g * x + c_g) = b_f * (b_g * x + c_g) + c_f
    return F{f.b * g.b, f.b * g.c + f.c};
}
F id() {
    return F{1, 0};
}

int main() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);

    int n, q;
    if (!(cin >> n >> q)) return 0;

    vector<S> a(n);
    for (int i = 0; i < n; i++) {
        long long v;
        cin >> v;
        a[i] = S{mint(v), 1};
    }

    lazysegtree<S, op, e, F, mapping, composition, id> seg(a);

    for (int i = 0; i < q; i++) {
        int type;
        cin >> type;
        if (type == 0) {
            int l, r;
            long long b, c;
            cin >> l >> r >> b >> c;
            seg.apply(l, r, F{mint(b), mint(c)});
        } else {
            int l, r;
            cin >> l >> r;
            // 【修正】末尾に () を追加して整数値を取得
            cout << seg.prod(l, r).val.val() << "\n";
        }
    }

    return 0;
}