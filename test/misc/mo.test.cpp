#include <algorithm>
#include <array>
#include <cassert>
#include <deque>
#include <random>
#include <utility>
#include <vector>

#include "../../src/misc/mo.hpp"
#include "../../src/misc/mo_2d.hpp"

using query = std::pair<int, int>;

std::array<int, 4> test_points(int nx, int ny, const std::vector<query>& queries) {
    mo_algorithm_2d solver(nx, ny);
    for (auto [x, y] : queries) solver.add_query(x, y);
    assert(solver.size() == int(queries.size()));

    std::vector<std::vector<int>> a(nx, std::vector<int>(ny));
    std::vector<std::vector<int>> prefix(nx + 1, std::vector<int>(ny + 1));
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            a[i][j] = (i * 13 + j * 7 + i * j) % 23 - 11;
            prefix[i + 1][j + 1] = a[i][j] + prefix[i][j + 1]
                                  + prefix[i + 1][j] - prefix[i][j];
        }
    }

    std::array<int, 4> moves{};
    for (int repeat = 0; repeat < 2; repeat++) {
        int x = 0, y = 0, sum = 0;
        std::vector<int> seen(queries.size());
        solver.run(
            [&](int i) {
                assert(i == x && x < nx);
                for (int j = 0; j < y; j++) sum += a[i][j];
                x++;
                moves[0]++;
            },
            [&](int j) {
                assert(j == y && y < ny);
                for (int i = 0; i < x; i++) sum += a[i][j];
                y++;
                moves[1]++;
            },
            [&](int i) {
                assert(x > 0 && i == x - 1);
                x--;
                for (int j = 0; j < y; j++) sum -= a[i][j];
                moves[2]++;
            },
            [&](int j) {
                assert(y > 0 && j == y - 1);
                y--;
                for (int i = 0; i < x; i++) sum -= a[i][j];
                moves[3]++;
            },
            [&](int id) {
                assert(0 <= id && id < int(queries.size()));
                assert(seen[id] == 0);
                seen[id]++;
                assert(queries[id] == query(x, y));
                assert(sum == prefix[x][y]);
            }
        );
        for (int count : seen) assert(count == 1);
    }
    return moves;
}

void test_intervals(int n, const std::vector<query>& queries) {
    mo_algorithm solver(n);
    for (auto [l, r] : queries) solver.add_query(l, r);
    int left = 0, right = 0;
    std::deque<int> values;
    solver.run(
        [&](int i) {
            left--;
            assert(i == left);
            values.push_front(i);
        },
        [&](int i) {
            assert(i == right);
            right++;
            values.push_back(i);
        },
        [&](int i) {
            assert(i == left && !values.empty() && values.front() == i);
            left++;
            values.pop_front();
        },
        [&](int i) {
            right--;
            assert(i == right && !values.empty() && values.back() == i);
            values.pop_back();
        },
        [&](int id) {
            assert(queries[id] == query(left, right));
            assert(int(values.size()) == right - left);
            for (int i = left; i < right; i++) assert(values[i - left] == i);
        }
    );

    std::vector<int> active(n);
    solver.run(
        [&](int i) { assert(active[i] == 0); active[i]++; },
        [&](int i) { assert(active[i] == 1); active[i]--; },
        [&](int id) {
            auto [l, r] = queries[id];
            for (int i = 0; i < n; i++) assert(active[i] == int(l <= i && i < r));
        }
    );
}

int main() {
    mo_algorithm empty;
    auto unexpected = [](int) { assert(false); };
    empty.run(unexpected, unexpected, unexpected);
    mo_algorithm_2d empty_2d;
    empty_2d.run(unexpected, unexpected, unexpected, unexpected, unexpected);
    empty_2d.add_query(0, 0);
    empty_2d.run(unexpected, unexpected, unexpected, unexpected, [](int id) { assert(id == 0); });

    auto moves = test_points(5, 7, {{0, 0}, {5, 0}, {5, 7}, {0, 7}, {3, 1}, {1, 6}, {4, 2}});
    for (int count : moves) assert(count > 0);
    for (auto [nx, ny] : std::vector<query>{{0, 0}, {0, 7}, {7, 0}, {1, 1}, {3, 17}, {19, 2}, {32, 32}}) {
        std::vector<query> queries;
        for (int x = 0; x <= nx; x++) {
            for (int y = 0; y <= ny; y++) queries.emplace_back(x, y);
        }
        test_points(nx, ny, queries);
    }

    std::mt19937 rng(459);
    for (int tc = 0; tc < 100; tc++) {
        int nx = int(rng() % 30), ny = int(rng() % 30);
        std::vector<query> queries;
        for (int i = 0; i < 100; i++) queries.emplace_back(rng() % (nx + 1), rng() % (ny + 1));
        std::shuffle(queries.begin(), queries.end(), rng);
        test_points(nx, ny, queries);
        test_points(nx, ny, {});

        int n = int(rng() % 40);
        queries = {{0, 0}, {n, n}, {0, n}};
        for (int i = 0; i < 100; i++) {
            int l = int(rng() % (n + 1)), r = int(rng() % (n + 1));
            if (l > r) std::swap(l, r);
            queries.emplace_back(l, r);
        }
        test_intervals(n, queries);
        test_intervals(n, {});
    }
}
