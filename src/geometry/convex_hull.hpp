#ifndef LIBRARY_GEOMETRY_CONVEX_HULL_HPP
#define LIBRARY_GEOMETRY_CONVEX_HULL_HPP

#pragma once

#include <algorithm>
#include <vector>
#include "geometry_template.hpp"

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
