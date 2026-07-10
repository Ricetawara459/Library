#pragma once

#include <vector>
#include <string>
#include <algorithm>

/// 各 i について s[0..] と s[i..] の LCP 長を返す。
template <class T>
std::vector<int> z_algorithm(const std::vector<T>& s) {
    int n = int(s.size());
    if (n == 0) return {};
    std::vector<int> z(n);
    z[0] = n;
    int i = 1, j = 0;
    while (i < n) {
        while (i + j < n && s[j] == s[i + j]) j++;
        z[i] = j;
        if (j == 0) {
            i++;
            continue;
        }
        int k = 1;
        while (i + k < n && k + z[k] < j) {
            z[i + k] = z[k];
            k++;
        }
        i += k;
        j -= k;
    }
    return z;
}

/// 文字列 s の Z-array を返す。
std::vector<int> z_algorithm(const std::string& s) {
    int n = int(s.size());
    std::vector<int> s2(n);
    for (int i = 0; i < n; i++) {
        s2[i] = s[i];
    }
    return z_algorithm(s2);
}
