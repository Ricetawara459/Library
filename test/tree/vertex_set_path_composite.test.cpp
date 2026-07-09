#define PROBLEM "https://judge.yosupo.jp/problem/vertex_set_path_composite"

#include <iostream>
#include <vector>

#include "../../src/data-structure/segment_tree.hpp"
#include "../../src/math/modint.hpp"
#include "../../src/tree/hld.hpp"

using namespace std;
using mint = modint998244353;

struct affine {
    mint a, b;
};

affine compose(affine g, affine f) {
    return {g.a * f.a, g.a * f.b + g.b};
}

struct node {
    affine lr, rl;
};

node op(node l, node r) {
    return {
        compose(r.lr, l.lr),
        compose(l.rl, r.rl),
    };
}

node e() {
    affine id{1, 0};
    return {id, id};
}

node make_node(mint a, mint b) {
    affine f{a, b};
    return {f, f};
}

int main() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);

    int n, q;
    cin >> n >> q;

    vector<mint> a(n), b(n);
    for (int i = 0; i < n; i++) {
        int c, d;
        cin >> c >> d;
        a[i] = c;
        b[i] = d;
    }

    hld tree(n);
    for (int i = 0; i < n - 1; i++) {
        int u, v;
        cin >> u >> v;
        tree.add_edge(u, v);
    }
    tree.build(0);

    vector<node> init(n);
    for (int v = 0; v < n; v++) {
        init[tree[v]] = make_node(a[v], b[v]);
    }
    segtree<node, op, e> seg(init);

    while (q--) {
        int type;
        cin >> type;
        if (type == 0) {
            int p, c, d;
            cin >> p >> c >> d;
            seg.set(tree[p], make_node(c, d));
        } else {
            int u, v, x;
            cin >> u >> v >> x;
            affine f{1, 0};
            tree.path_query_ordered(u, v, false, [&](int l, int r, bool rev) {
                node cur = seg.prod(l, r);
                f = compose(rev ? cur.rl : cur.lr, f);
            });
            cout << (f.a * x + f.b).val() << '\n';
        }
    }
}
