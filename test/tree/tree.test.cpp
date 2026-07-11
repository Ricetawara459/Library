#include <cassert>
#include <vector>

#include "../../src/tree/lca.hpp"
#include "../../src/tree/hld.hpp"
#include "../../src/data-structure/segtree.hpp"

int op_sum(int a, int b) { return a + b; }
int e_sum() { return 0; }

int main() {
    lca lc(7);
    hld hd(7);
    std::vector<std::pair<int, int>> edges = {
        {0, 1}, {0, 2}, {1, 3}, {1, 4}, {2, 5}, {2, 6},
    };
    for (auto [u, v] : edges) {
        lc.add_edge(u, v);
        hd.add_edge(u, v);
    }
    lc.build(0);
    hd.build(0);

    assert(lc.query(3, 4) == 1);
    assert(lc.query(3, 6) == 0);
    assert(lc.distance(3, 6) == 4);
    assert(lc.jump(3, 6, 2) == 0);
    assert(lc.kth_ancestor(6, 2) == 0);

    assert(hd.lca(3, 4) == 1);
    assert(hd.distance(3, 6) == 4);

    std::vector<int> base(7);
    for (int v = 0; v < 7; v++) base[hd[v]] = v + 1;
    segtree<int, op_sum, e_sum> seg(base);

    int path_sum = 0;
    hd.path_query(3, 6, false, [&](int l, int r) {
        path_sum += seg.prod(l, r);
    });
    assert(path_sum == 4 + 2 + 1 + 3 + 7);

    auto [l, r] = hd.subtree_query(1);
    assert(seg.prod(l, r) == 2 + 4 + 5);
}
