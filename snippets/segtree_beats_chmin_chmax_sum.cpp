// Range Chmin Chmax Add Range Sum 用の Segment Tree Beats チートシート。
// src/data-structure/segtree_beats.hpp を include して使う。

#include <algorithm>
#include <limits>
#include <vector>

using ll = long long;

namespace beats_chmin_chmax_sum {
constexpr ll INF = std::numeric_limits<ll>::max() / 4;

struct S {
    ll sum;
    ll max1, max2;
    int max_count;
    ll min1, min2;
    int min_count;
    int len;
};

S make_node(ll x) { return {x, x, -INF, 1, x, INF, 1, 1}; }

S e() { return {0, -INF, -INF, 0, INF, INF, 0, 0}; }

S op(S a, S b) {
    if (a.len == 0) return b;
    if (b.len == 0) return a;
    S c{};
    c.sum = a.sum + b.sum;
    c.len = a.len + b.len;

    if (a.max1 == b.max1) {
        c.max1 = a.max1;
        c.max2 = std::max(a.max2, b.max2);
        c.max_count = a.max_count + b.max_count;
    } else if (a.max1 > b.max1) {
        c.max1 = a.max1;
        c.max2 = std::max(a.max2, b.max1);
        c.max_count = a.max_count;
    } else {
        c.max1 = b.max1;
        c.max2 = std::max(a.max1, b.max2);
        c.max_count = b.max_count;
    }

    if (a.min1 == b.min1) {
        c.min1 = a.min1;
        c.min2 = std::min(a.min2, b.min2);
        c.min_count = a.min_count + b.min_count;
    } else if (a.min1 < b.min1) {
        c.min1 = a.min1;
        c.min2 = std::min(a.min2, b.min1);
        c.min_count = a.min_count;
    } else {
        c.min1 = b.min1;
        c.min2 = std::min(a.min1, b.min2);
        c.min_count = b.min_count;
    }
    return c;
}

struct F {
    ll lower, upper, add;
};

ll clamp_value(ll x, ll lower, ll upper) {
    return std::min(std::max(x, lower), upper);
}

ll shift_lower(ll x, ll add) {
    return x == -INF ? -INF : x + add;
}

ll shift_upper(ll x, ll add) {
    return x == INF ? INF : x + add;
}

F id() { return {-INF, INF, 0}; }

F composition(F f, F g) {
    return {clamp_value(shift_lower(g.lower, f.add), f.lower, f.upper),
            clamp_value(shift_upper(g.upper, f.add), f.lower, f.upper),
            g.add + f.add};
}

void apply_add(S& x, ll a) {
    if (x.len == 0 || a == 0) return;
    x.sum += a * x.len;
    x.max1 += a;
    if (x.max2 != -INF) x.max2 += a;
    x.min1 += a;
    if (x.min2 != INF) x.min2 += a;
}

bool apply_chmin(S& x, ll a) {
    if (x.len == 0 || x.max1 <= a) return true;
    if (a <= x.max2) return false;
    x.sum += (a - x.max1) * x.max_count;
    if (x.min1 == x.max1) x.min1 = a;
    else if (x.min2 == x.max1) x.min2 = a;
    x.max1 = a;
    return true;
}

bool apply_chmax(S& x, ll a) {
    if (x.len == 0 || a <= x.min1) return true;
    if (x.min2 <= a) return false;
    x.sum += (a - x.min1) * x.min_count;
    if (x.max1 == x.min1) x.max1 = a;
    else if (x.max2 == x.min1) x.max2 = a;
    x.min1 = a;
    return true;
}

bool mapping(F f, S& x) {
    S y = x;
    apply_add(y, f.add);
    if (!apply_chmax(y, f.lower)) return false;
    if (!apply_chmin(y, f.upper)) return false;
    x = y;
    return true;
}

F add_query(ll x) { return {-INF, INF, x}; }
F chmin_query(ll x) { return {-INF, x, 0}; }
F chmax_query(ll x) { return {x, INF, 0}; }

using segtree = segtree_beats<S, op, e, F, mapping, composition, id>;

std::vector<S> make_nodes(const std::vector<ll>& a) {
    std::vector<S> v;
    v.reserve(a.size());
    for (ll x : a) v.push_back(make_node(x));
    return v;
}

struct range_chmin_chmax_add_sum {
    segtree seg;

    range_chmin_chmax_add_sum() = default;
    explicit range_chmin_chmax_add_sum(const std::vector<ll>& a)
        : seg(make_nodes(a)) {}

    // a[p] = x。
    void set(int p, ll x) { seg.set(p, make_node(x)); }
    // a[p] を返す。
    ll get(int p) { return seg.get(p).sum; }
    // 半開区間 [l, r) に x を加算する。
    void range_add(int l, int r, ll x) { seg.apply(l, r, add_query(x)); }
    // 半開区間 [l, r) の各値を min(a[i], x) にする。
    void range_chmin(int l, int r, ll x) { seg.apply(l, r, chmin_query(x)); }
    // 半開区間 [l, r) の各値を max(a[i], x) にする。
    void range_chmax(int l, int r, ll x) { seg.apply(l, r, chmax_query(x)); }
    // 半開区間 [l, r) の総和を返す。
    ll range_sum(int l, int r) { return seg.prod(l, r).sum; }
    // 全体の総和を返す。
    ll all_sum() { return seg.all_prod().sum; }
};
}  // namespace beats_chmin_chmax_sum
