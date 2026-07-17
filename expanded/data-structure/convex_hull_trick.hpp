#ifndef LIBRARY_DATA_STRUCTURE_CONVEX_HULL_TRICK_HPP
#define LIBRARY_DATA_STRUCTURE_CONVEX_HULL_TRICK_HPP

#line 1 "src/data-structure/convex_hull_trick.hpp"


#include <cassert>
#include <limits>
#include <set>

namespace cht_internal {

template <class T>
struct line_container {
  private:
    struct line {
        mutable T a, b; // y = ax + b
        mutable T p;    // 次の直線に負け始める最初の x

        bool operator<(const line& other) const { return a < other.a; }
        bool operator<(T x) const { return p < x; }
    };

    static constexpr T inf = std::numeric_limits<T>::max();

    std::multiset<line, std::less<>> lines;

    // floor(a / b)。b != 0。
    static T floor_div(T a, T b) {
        assert(b != 0);
        if (b < 0) a = -a, b = -b;
        if (a >= 0) return a / b;
        return -((-a + b - 1) / b);
    }

    // x の直線と y の直線の交点を見て、x.p を更新する。
    bool intersect(typename std::multiset<line, std::less<>>::iterator x,
                   typename std::multiset<line, std::less<>>::iterator y) {
        if (y == lines.end()) {
            x->p = inf;
            return false;
        }
        if (x->a == y->a) {
            x->p = (x->b > y->b ? inf : -inf);
        } else {
            x->p = floor_div(y->b - x->b, x->a - y->a);
        }
        return x->p >= y->p;
    }

  public:
    line_container() = default;

    void add(T a, T b) {
        auto z = lines.insert({a, b, 0});
        auto y = z++;
        auto x = y;
        while (intersect(y, z)) z = lines.erase(z);
        if (x != lines.begin() && intersect(--x, y)) {
            intersect(x, y = lines.erase(y));
        }
        while ((y = x) != lines.begin() && (--x)->p >= y->p) {
            intersect(x, lines.erase(y));
        }
    }

    bool empty() const { return lines.empty(); }
    int size() const { return int(lines.size()); }

    T query(T x) const {
        assert(!lines.empty());
        auto l = *lines.lower_bound(x);
        return l.a * x + l.b;
    }
};

} // namespace cht_internal

template <class T>
struct cht_max {
  public:
    cht_max() = default;

    /// 最大値クエリ用に直線 y = ax + b を追加する。任意の順番で追加してよい。
    void add(T a, T b) { hull.add(a, b); }

    /// x における最大値 max(ax+b) を返す。
    T query(T x) const { return hull.query(x); }

    /// 直線が 1 本以上入っているかを返す。
    bool empty() const { return hull.empty(); }
    /// 入っている直線の本数を返す。不要になった直線は内部で消える。
    int size() const { return hull.size(); }

  private:
    cht_internal::line_container<T> hull;
};

template <class T>
struct cht_min {
  public:
    cht_min() = default;

    /// 最小値クエリ用に直線 y = ax + b を追加する。任意の順番で追加してよい。
    void add(T a, T b) { hull.add(-a, -b); }

    /// x における最小値 min(ax+b) を返す。
    T query(T x) const { return -hull.query(x); }

    /// 直線が 1 本以上入っているかを返す。
    bool empty() const { return hull.empty(); }
    /// 入っている直線の本数を返す。不要になった直線は内部で消える。
    int size() const { return hull.size(); }

  private:
    cht_internal::line_container<T> hull;
};

#endif  // LIBRARY_DATA_STRUCTURE_CONVEX_HULL_TRICK_HPP
