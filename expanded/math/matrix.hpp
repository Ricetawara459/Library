#ifndef LIBRARY_MATH_MATRIX_HPP
#define LIBRARY_MATH_MATRIX_HPP

#line 1 "src/math/matrix.hpp"


#include <vector>
#include <cassert>
#include <iostream>

template <class T>
struct matrix {
  public:
    int h, w;
    std::vector<std::vector<T>> data;

    matrix() : h(0), w(0) {}
    matrix(int h, int w, T val = T(0)) : h(h), w(w), data(h, std::vector<T>(w, val)) {}
    matrix(const std::vector<std::vector<T>>& d) : h(d.size()), w(d.empty() ? 0 : d[0].size()), data(d) {}

    /// n x n の単位行列を返す。
    static matrix identity(int n) {
        matrix res(n, n);
        for (int i = 0; i < n; i++) res[i][i] = T(1);
        return res;
    }

    /// 行数を返す。
    int height() const { return h; }
    /// 列数を返す。
    int width() const { return w; }

    const std::vector<T>& operator[](int i) const { return data[i]; }
    std::vector<T>& operator[](int i) { return data[i]; }

    matrix& operator+=(const matrix& other) {
        assert(h == other.h && w == other.w);
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                data[i][j] += other.data[i][j];
            }
        }
        return *this;
    }

    matrix& operator-=(const matrix& other) {
        assert(h == other.h && w == other.w);
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                data[i][j] -= other.data[i][j];
            }
        }
        return *this;
    }

    matrix& operator*=(const matrix& other) {
        assert(w == other.h);
        matrix res(h, other.w);
        for (int i = 0; i < h; i++) {
            for (int k = 0; k < w; k++) {
                for (int j = 0; j < other.w; j++) {
                    res.data[i][j] += data[i][k] * other.data[k][j];
                }
            }
        }
        return *this = res;
    }

    matrix operator+(const matrix& other) const { return matrix(*this) += other; }
    matrix operator-(const matrix& other) const { return matrix(*this) -= other; }
    matrix operator*(const matrix& other) const { return matrix(*this) *= other; }

    /// 行列式を返す。T は割り算ができる体であること。計算量 O(n^3)。
    T determinant() const {
        assert(h == w);
        matrix a = *this;
        T det = T(1);
        for (int i = 0; i < h; i++) {
            int pivot = -1;
            for (int j = i; j < h; j++) {
                if (a[j][i] != T(0)) {
                    pivot = j;
                    break;
                }
            }
            if (pivot == -1) return T(0);
            if (pivot != i) {
                std::swap(a[pivot], a[i]);
                det = -det;
            }
            det *= a[i][i];
            T inv = T(1) / a[i][i];
            for (int j = i + 1; j < h; j++) {
                T factor = a[j][i] * inv;
                if (factor == T(0)) continue;
                for (int k = i; k < h; k++) {
                    a[j][k] -= factor * a[i][k];
                }
            }
        }
        return det;
    }

    /// 正方行列の k 乗を返す。k >= 0。
    matrix pow(long long k) const {
        assert(h == w);
        matrix res = identity(h);
        matrix base = *this;
        while (k > 0) {
            if (k & 1) res *= base;
            base *= base;
            k >>= 1;
        }
        return res;
    }

    // デバッグ出力用演算子のオーバーロード
    friend std::ostream& operator<<(std::ostream& os, const matrix& mat) {
        for (int i = 0; i < mat.h; i++) {
            for (int j = 0; j < mat.w; j++) {
                os << mat.data[i][j] << (j + 1 == mat.w ? "" : " ");
            }
            if (i + 1 < mat.h) os << "\n";
        }
        return os;
    }
};

#endif  // LIBRARY_MATH_MATRIX_HPP
