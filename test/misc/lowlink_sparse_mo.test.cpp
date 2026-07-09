#include <algorithm>
#include <cassert>
#include <vector>

#include "../../src/graph/lowlink.hpp"
#include "../../src/data-structure/sparse_table.hpp"
#include "../../src/misc/mo.hpp"

int op_min(int a, int b) { return std::min(a, b); }

int main() {
    lowlink ll(5);
    ll.add_edge(0, 1);
    ll.add_edge(1, 2);
    ll.add_edge(2, 0);
    ll.add_edge(1, 3);
    ll.add_edge(3, 4);
    ll.build();
    assert((ll.articulation == std::vector<int>{1, 3}));
    assert((ll.bridges == std::vector<std::pair<int, int>>{{1, 3}, {3, 4}}));

    std::vector<int> a = {5, 2, 6, 3, 1, 4};
    sparse_table<int, op_min> st(a);
    assert(st.prod(0, 6) == 1);
    assert(st.prod(1, 4) == 2);
    assert(st.prod(2, 4) == 3);

    mo solver(int(a.size()));
    solver.add_query(0, 3);
    solver.add_query(2, 6);
    solver.add_query(1, 5);
    std::vector<int> ans(3);
    int sum = 0;
    solver.run(
        [&](int i) { sum += a[i]; },
        [&](int i) { sum -= a[i]; },
        [&](int id) { ans[id] = sum; }
    );
    assert((ans == std::vector<int>{13, 14, 12}));
}
