// Source: https://github.com/Ricetawara459/Library
#ifndef LIBRARY_DATA_STRUCTURE_WAVELET_MATRIX_HPP
#define LIBRARY_DATA_STRUCTURE_WAVELET_MATRIX_HPP



#include <algorithm>
#include <cassert>
#include <cstdint>
#include <optional>
#include <vector>

/// 静的配列の Wavelet Matrix。値は内部で座標圧縮する。
/// 区間は [l, r)、順位 k は 0-indexed。構築 O(n log n)、各クエリ O(log n)。
template <class T>
struct wavelet_matrix {
  private:
    int n = 0, levels = 1;
    std::vector<T> values;
    std::vector<std::vector<std::uint64_t>> bits;
    std::vector<std::vector<int>> sums;
    std::vector<int> zeros;

    int rank1(int level, int pos) const {
        int block = pos >> 6, offset = pos & 63;
        int result = sums[level][block];
        if (offset) result += __builtin_popcountll(bits[level][block] & ((std::uint64_t(1) << offset) - 1));
        return result;
    }
    void check_range(int l, int r) const { assert(0 <= l && l <= r && r <= n); }
    int less_id(int l, int r, int upper) const {
        if (upper <= 0) return 0;
        if (upper >= int(values.size())) return r - l;
        int result = 0;
        for (int level = levels - 1; level >= 0; level--) {
            int a = rank1(level, l), b = rank1(level, r);
            if ((upper >> level) & 1) {
                result += (r - l) - (b - a);
                l = zeros[level] + a; r = zeros[level] + b;
            } else {
                l -= a; r -= b;
            }
        }
        return result;
    }

  public:
    wavelet_matrix() { build({}); }
    /// 配列 v から構築する。値の型 T は比較可能であること。
    explicit wavelet_matrix(const std::vector<T>& v) { build(v); }

    /// 配列 v で再構築する。
    void build(const std::vector<T>& v) {
        n = int(v.size());
        values = v;
        std::sort(values.begin(), values.end());
        values.erase(std::unique(values.begin(), values.end()), values.end());
        levels = 1;
        while ((std::uint64_t(1) << levels) < values.size()) levels++;
        bits.assign(levels, std::vector<std::uint64_t>((n + 63) / 64));
        sums.assign(levels, std::vector<int>((n + 63) / 64 + 1));
        zeros.assign(levels, 0);
        std::vector<int> a(n), next(n);
        for (int i = 0; i < n; i++) a[i] = int(std::lower_bound(values.begin(), values.end(), v[i]) - values.begin());
        for (int level = levels - 1; level >= 0; level--) {
            for (int i = 0; i < n; i++) {
                if ((a[i] >> level) & 1) bits[level][i >> 6] |= std::uint64_t(1) << (i & 63);
                else zeros[level]++;
            }
            for (int b = 0; b < int(bits[level].size()); b++)
                sums[level][b + 1] = sums[level][b] + __builtin_popcountll(bits[level][b]);
            int z = 0, o = zeros[level];
            for (int x : a) next[((x >> level) & 1) ? o++ : z++] = x;
            a.swap(next);
        }
    }

    /// 配列の長さを返す。
    int size() const { return n; }
    /// a[k] を返す。0 <= k < size()。
    T access(int k) const {
        assert(0 <= k && k < n);
        int id = 0;
        for (int level = levels - 1; level >= 0; level--) {
            int b = int(bits[level][k >> 6] >> (k & 63) & 1);
            int count = rank1(level, k);
            if (b) { id |= 1 << level; k = zeros[level] + count; }
            else k -= count;
        }
        return values[id];
    }
    /// a[k] を返す。
    T operator[](int k) const { return access(k); }

    /// [l, r) 内の値 x の個数を返す。
    int count(int l, int r, const T& x) const {
        check_range(l, r);
        auto it = std::lower_bound(values.begin(), values.end(), x);
        if (it == values.end() || !(*it == x)) return 0;
        int id = int(it - values.begin());
        return less_id(l, r, id + 1) - less_id(l, r, id);
    }
    /// [0, r) 内の値 x の個数を返す。
    int rank(const T& x, int r) const { return count(0, r, x); }

    /// [l, r) 内の小さい方から k 番目の値を返す。0 <= k < r-l。
    T kth_smallest(int l, int r, int k) const {
        check_range(l, r);
        assert(0 <= k && k < r - l);
        int id = 0;
        for (int level = levels - 1; level >= 0; level--) {
            int a = rank1(level, l), b = rank1(level, r);
            int count = (r - l) - (b - a);
            if (k < count) { l -= a; r -= b; }
            else { k -= count; id |= 1 << level; l = zeros[level] + a; r = zeros[level] + b; }
        }
        return values[id];
    }
    /// [l, r) 内の大きい方から k 番目の値を返す。0 <= k < r-l。
    T kth_largest(int l, int r, int k) const {
        check_range(l, r);
        assert(0 <= k && k < r - l);
        return kth_smallest(l, r, r - l - 1 - k);
    }
    /// [l, r) 内で upper 未満の値の個数を返す。
    int range_freq(int l, int r, const T& upper) const {
        check_range(l, r);
        int id = int(std::lower_bound(values.begin(), values.end(), upper) - values.begin());
        return less_id(l, r, id);
    }
    /// [l, r) 内で lower <= 値 < upper を満たす個数を返す。lower <= upper。
    int range_freq(int l, int r, const T& lower, const T& upper) const {
        assert(!(upper < lower));
        return range_freq(l, r, upper) - range_freq(l, r, lower);
    }
    /// [l, r) 内で upper 未満の最大値。存在しなければ nullopt。
    std::optional<T> prev_value(int l, int r, const T& upper) const {
        int count = range_freq(l, r, upper);
        if (count == 0) return std::nullopt;
        return kth_smallest(l, r, count - 1);
    }
    /// [l, r) 内で lower 以上の最小値。存在しなければ nullopt。
    std::optional<T> next_value(int l, int r, const T& lower) const {
        int count = range_freq(l, r, lower);
        if (count == r - l) return std::nullopt;
        return kth_smallest(l, r, count);
    }
};

#endif  // LIBRARY_DATA_STRUCTURE_WAVELET_MATRIX_HPP
