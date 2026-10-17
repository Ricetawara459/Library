#include <algorithm>
#include <cassert>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "../../src/data-structure/fenwick_tree.hpp"
#include "../../src/data-structure/lazysegtree.hpp"
#include "../../src/data-structure/segtree.hpp"
#include "../../src/data-structure/segtree_beats.hpp"
#include "../../src/data-structure/sparse_table.hpp"
#include "../../src/math/modint.hpp"
#include "../../snippets/segtree_beats_chmin_chmax_sum.cpp"

using value_type = long long;

value_type add_values(value_type a, value_type b) { return a + b; }
value_type zero() { return 0; }
value_type min_value(value_type a, value_type b) { return std::min(a, b); }

struct node {
    value_type sum;
    int len;

    friend std::ostream& operator<<(std::ostream& os, const node& x) {
        return os << x.sum;
    }
};

struct affine {
    value_type mul, add;
};

node merge(node a, node b) { return {a.sum + b.sum, a.len + b.len}; }
node empty_node() { return {0, 0}; }
int mapping_calls = 0;
node map_node(affine f, node x) {
    mapping_calls++;
    return {f.mul * x.sum + f.add * x.len, x.len};
}
bool map_beats(affine f, node& x) {
    x = map_node(f, x);
    return true;
}
affine compose(affine f, affine g) { return {f.mul * g.mul, f.mul * g.add + f.add}; }
affine identity() { return {1, 0}; }

using sum_tree = segtree<value_type, add_values, zero>;
using lazy_tree = lazysegtree<node, merge, empty_node, affine, map_node, compose, identity>;
using beats_tree = segtree_beats<node, merge, empty_node, affine, map_beats, compose, identity>;

template <class T>
std::string dump(const T& value) {
    std::ostringstream out;
    out << value;
    return out.str();
}

std::string expected_output(const std::vector<value_type>& a) {
    std::ostringstream out;
    out << '[';
    for (int i = 0; i < int(a.size()); i++) {
        if (i) out << ", ";
        out << a[i];
    }
    out << ']';
    return out.str();
}

int main() {
    assert(dump(sum_tree{}) == "[]");
    assert(dump(lazy_tree{}) == "[]");
    assert(dump(beats_tree{}) == "[]");
    assert(dump(fenwick_tree<value_type>{}) == "[]");
    assert(dump(sparse_table<value_type, min_value>{}) == "[]");
    assert(dump(beats_chmin_chmax_sum::range_chmin_chmax_add_sum{}) == "[]");

    std::mt19937 rng(459);
    for (int n : {0, 1, 2, 3, 5, 16, 31, 257}) {
        std::vector<value_type> a(n);
        std::vector<node> initial(n);
        fenwick_tree<value_type> bit(n);
        for (int i = 0; i < n; i++) {
            a[i] = int(rng() % 101) - 50;
            initial[i] = {a[i], 1};
            bit.add(i, a[i]);
        }
        sum_tree seg(a);
        lazy_tree lazy(initial);
        beats_tree beats(initial);
        sparse_table<value_type, min_value> table(a);

        for (int q = 0; q < 60; q++) {
            int l = int(rng() % (n + 1)), r = int(rng() % (n + 1));
            if (l > r) std::swap(l, r);
            affine f{int(rng() % 3) - 1, int(rng() % 21) - 10};
            lazy.apply(l, r, f);
            beats.apply(l, r, f);
            for (int i = l; i < r; i++) {
                value_type next = f.mul * a[i] + f.add;
                bit.add(i, next - a[i]);
                a[i] = next;
                seg.set(i, next);
            }
            if (n && q % 3 == 0) {
                int p = int(rng() % n);
                value_type next = int(rng() % 101) - 50;
                bit.add(p, next - a[p]);
                a[p] = next;
                seg.set(p, next);
                lazy.set(p, {next, 1});
                beats.set(p, {next, 1});
            }
            table.build(a);
            std::string expected = expected_output(a);
            assert(dump(seg) == expected);
            assert(dump(bit) == expected);
            assert(dump(table) == expected);
            mapping_calls = 0;
            assert(dump(lazy) == expected);
            assert(mapping_calls <= 4 * n);
            mapping_calls = 0;
            assert(dump(beats) == expected);
            assert(mapping_calls <= 4 * n);
            assert(dump(lazy) == expected);
            assert(dump(beats) == expected);

            value_type sum = std::accumulate(a.begin(), a.end(), 0LL);
            assert(seg.all_prod() == sum);
            assert(bit.sum(0, n) == sum);
            assert(lazy.all_prod().sum == sum);
            assert(beats.all_prod().sum == sum);
            if (n) assert(table.prod(0, n) == *std::min_element(a.begin(), a.end()));
        }
        std::ostringstream chained;
        chained << seg << '!' << bit;
        assert(chained.str() == expected_output(a) + '!' + expected_output(a));
    }

    for (int n : {1, 5, 16, 31}) {
        std::vector<value_type> a(n);
        for (auto& x : a) x = int(rng() % 101) - 50;
        beats_chmin_chmax_sum::range_chmin_chmax_add_sum seg(a);
        for (int q = 0; q < 300; q++) {
            int l = int(rng() % (n + 1)), r = int(rng() % (n + 1));
            if (l > r) std::swap(l, r);
            value_type x = int(rng() % 101) - 50;
            int type = int(rng() % 3);
            if (type == 0) seg.range_add(l, r, x);
            if (type == 1) seg.range_chmin(l, r, x);
            if (type == 2) seg.range_chmax(l, r, x);
            for (int i = l; i < r; i++) {
                if (type == 0) a[i] += x;
                if (type == 1) a[i] = std::min(a[i], x);
                if (type == 2) a[i] = std::max(a[i], x);
            }
            if (q % 17 == 0) {
                assert(dump(seg) == expected_output(a));
                assert(dump(seg.seg) == expected_output(a));
            }
            assert(seg.all_sum() == std::accumulate(a.begin(), a.end(), 0LL));
        }
        assert(dump(seg) == expected_output(a));
    }

    fenwick_tree<modint998244353> modular(3);
    modular.add(0, -1);
    modular.add(1, 2);
    assert(dump(modular) == "[998244352, 2, 0]");
    assert(modular.sum(0, 3) == modint998244353(1));
}
