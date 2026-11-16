// Source: https://github.com/Ricetawara459/Library
#ifndef LIBRARY_GEOMETRY_CONVEX_HULL_HPP
#define LIBRARY_GEOMETRY_CONVEX_HULL_HPP



#include <algorithm>
#include <vector>
#ifndef LIBRARY_GEOMETRY_GEOMETRY_TEMPLATE_HPP
#define LIBRARY_GEOMETRY_GEOMETRY_TEMPLATE_HPP



#include <cmath>
#include <iostream>

using Real = double;
const Real EPS = 1e-9;

inline int sign(Real r) {
    return r < -EPS ? -1 : (r > EPS ? 1 : 0);
}

struct Point {
    Real x, y;
    Point operator+(const Point& p) const { return {x + p.x, y + p.y}; }
    Point operator-(const Point& p) const { return {x - p.x, y - p.y}; }
    friend std::ostream& operator<<(std::ostream& os, const Point& p) {
        return os << p.x << ' ' << p.y;
    }
};

// 内積 (Dot Product)
Real dot(const Point& a, const Point& b) {
    return a.x * b.x + a.y * b.y;
}

// 外積 (Cross Product)
Real cross(const Point& a, const Point& b) {
    return a.x * b.y - a.y * b.x;
}

// ccw の戻り値の定義
// 状態の区別がつけば値自体は何でも良いですが、標準的な割り当てです
enum {
    COUNTER_CLOCKWISE = 1,  // 反時計回り
    CLOCKWISE = -1,         // 時計回り
    ONLINE_BACK = 2,        // c, a, b の順で同一直線上
    ONLINE_FRONT = -2,      // a, b, c の順で同一直線上
    ON_SEGMENT = 0          // a, c, b の順で線分上（端点含む）
};

int ccw(Point a, Point b, Point c) {
    b = b - a; 
    c = c - a;
    
    // 外積の符号で反時計回り・時計回りを判定
    if (sign(cross(b, c)) > 0) return COUNTER_CLOCKWISE;
    if (sign(cross(b, c)) < 0) return CLOCKWISE;
    
    // 同一直線上にある場合、内積の符号と長さで判定
    if (sign(dot(b, c)) < 0) return ONLINE_BACK;
    if (sign(dot(b, b) - dot(c, c)) < 0) return ONLINE_FRONT;
    
    return ON_SEGMENT;
}

#endif  // LIBRARY_GEOMETRY_GEOMETRY_TEMPLATE_HPP


/// 点集合の凸包を反時計回りに返す。始点は (x, y) が辞書順最小の点。
/// 重複点は除く。include_collinear=true なら辺上の途中の点も含める。
/// 全点が一直線上なら、false では両端のみ、true では辞書順に全点を返す。
/// 回転方向の判定には geometry_template.hpp の EPS を使う。計算量 O(n log n)。
std::vector<Point> convex_hull(std::vector<Point> points, bool include_collinear = false) {
    std::sort(points.begin(), points.end(), [](const Point& a, const Point& b) {
        return a.x != b.x ? a.x < b.x : a.y < b.y;
    });
    points.erase(std::unique(points.begin(), points.end(), [](const Point& a, const Point& b) {
        return a.x == b.x && a.y == b.y;
    }), points.end());
    int n = int(points.size());
    if (n <= 2) return points;
    auto turn = [](const Point& a, const Point& b, const Point& c) {
        return sign(cross(b - a, c - a));
    };
    if (include_collinear) {
        bool line = true;
        for (int i = 1; i + 1 < n; i++) {
            if (turn(points.front(), points.back(), points[i]) != 0) { line = false; break; }
        }
        if (line) return points;
    }
    std::vector<Point> hull;
    hull.reserve(2 * n);
    auto add = [&](const Point& p, int lower_size) {
        while (int(hull.size()) >= lower_size + 2) {
            int t = turn(hull[hull.size() - 2], hull.back(), p);
            if (include_collinear ? t >= 0 : t > 0) break;
            hull.pop_back();
        }
        hull.push_back(p);
    };
    for (const auto& p : points) add(p, 0);
    int lower_size = int(hull.size()) - 1;
    for (int i = n - 2; i >= 0; i--) add(points[i], lower_size);
    hull.pop_back();
    return hull;
}

#endif  // LIBRARY_GEOMETRY_CONVEX_HULL_HPP
