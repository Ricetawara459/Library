#include <algorithm>
#include <cassert>
#include <random>
#include <vector>
#include "../../src/geometry/convex_hull.hpp"

bool same(Point a, Point b) { return a.x == b.x && a.y == b.y; }
long long area(Point a, Point b, Point c) { return static_cast<long long>(cross(b - a, c - a)); }
std::vector<Point> reference(std::vector<Point> p) {
    std::sort(p.begin(), p.end(), [](Point a, Point b) { return a.x != b.x ? a.x < b.x : a.y < b.y; });
    p.erase(std::unique(p.begin(), p.end(), same), p.end());
    if (p.size() <= 2) return p;
    std::vector<Point> ans;
    int at = 0;
    do {
        ans.push_back(p[at]);
        int next = (at + 1) % p.size();
        for (int i = 0; i < int(p.size()); i++) {
            auto t = area(p[at], p[next], p[i]);
            auto d = [](Point a, Point b) { return dot(a - b, a - b); };
            if (t < 0 || (t == 0 && d(p[at], p[i]) > d(p[at], p[next]))) next = i;
        }
        at = next;
    } while (at != 0);
    return ans;
}

int main() {
    assert(convex_hull({}).empty());
    auto one = convex_hull({{2, 3}, {2, 3}});
    assert(one.size() == 1 && same(one[0], {2, 3}));
    std::vector<Point> line{{2, 2}, {0, 0}, {1, 1}, {1, 1}, {-1, -1}};
    auto ends = convex_hull(line);
    auto all = convex_hull(line, true);
    assert(ends.size() == 2 && same(ends[0], {-1, -1}) && same(ends[1], {2, 2}));
    assert(all.size() == 4 && same(all[0], {-1, -1}) && same(all[3], {2, 2}));
    std::vector<Point> square{{0,0},{1,0},{2,0},{2,1},{2,2},{1,2},{0,2},{0,1},{1,1},{0,0}};
    assert(convex_hull(square).size() == 4);
    assert(convex_hull(square, true).size() == 8);
    std::vector<Point> almost_line{{0, 0}, {1, 0}, {0.5, 1e-10}};
    assert(convex_hull(almost_line).size() == 2);
    assert(convex_hull(almost_line, true).size() == 3);
    std::mt19937 rng(459);
    for (int tc = 0; tc < 5000; tc++) {
        std::vector<Point> p(rng() % 35);
        for (auto& x : p) x = {double(int(rng() % 15) - 7), double(int(rng() % 15) - 7)};
        auto actual = convex_hull(p), expected = reference(p);
        assert(actual.size() == expected.size());
        for (int i = 0; i < int(actual.size()); i++) assert(same(actual[i], expected[i]));
        auto boundary = convex_hull(p, true);
        if (actual.size() >= 3) {
            for (int i = 0; i < int(boundary.size()); i++) {
                int j = (i + 1) % boundary.size();
                for (auto x : p) assert(area(boundary[i], boundary[j], x) >= 0);
            }
            for (auto x : p) {
                bool on = false;
                for (int i = 0; i < int(actual.size()); i++) {
                    Point a = actual[i], b = actual[(i + 1) % actual.size()];
                    if (area(a, b, x) == 0 && dot(x - a, x - b) <= 0) on = true;
                }
                if (on) assert(std::find_if(boundary.begin(), boundary.end(), [&](Point y) { return same(x, y); }) != boundary.end());
            }
        }
    }
}
