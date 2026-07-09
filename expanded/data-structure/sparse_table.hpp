#line 2 "src/data-structure/sparse_table.hpp"

#include <cassert>
#include <vector>

// Sparse Table。min/max/gcd などの冪等な演算に使う。
// prod(l, r) で半開区間 [l, r) を取得する。空区間は扱わないので l < r が必要。
template <class S, S (*op)(S, S)>
struct sparse_table {
  public:
    sparse_table() = default;
    explicit sparse_table(const std::vector<S>& v) {
        build(v);
    }

    void build(const std::vector<S>& v) {
        _n = int(v.size());
        log_table.assign(_n + 1, 0);
        for (int i = 2; i <= _n; i++) log_table[i] = log_table[i >> 1] + 1;

        table.assign(log_table[_n] + 1, std::vector<S>(_n));
        if (_n == 0) return;
        table[0] = v;
        for (int k = 1; k < int(table.size()); k++) {
            int len = 1 << k;
            for (int i = 0; i + len <= _n; i++) {
                table[k][i] = op(table[k - 1][i], table[k - 1][i + (len >> 1)]);
            }
        }
    }

    S prod(int l, int r) const {
        assert(0 <= l && l < r && r <= _n);
        int k = log_table[r - l];
        return op(table[k][l], table[k][r - (1 << k)]);
    }

    int size() const { return _n; }

  private:
    int _n = 0;
    std::vector<int> log_table;
    std::vector<std::vector<S>> table;
};
