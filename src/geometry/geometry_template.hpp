#pragma once

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
