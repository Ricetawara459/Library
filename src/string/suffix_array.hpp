#pragma once

#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <cassert>

template <class T>
std::vector<int> suffix_array(const std::vector<T>& s) {
    int n = int(s.size());
    if (n == 0) return {};
    
    // 最初のランク付けのために座標圧縮
    std::vector<T> vals = s;
    std::sort(vals.begin(), vals.end());
    vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
    
    std::vector<int> sa(n), rank(n), next_rank(n), sa_tmp(n);
    for (int i = 0; i < n; i++) {
        sa[i] = i;
        rank[i] = std::lower_bound(vals.begin(), vals.end(), s[i]) - vals.begin();
    }
    
    int max_val = std::max(256, n);
    std::vector<int> count(max_val + 1, 0);
    
    for (int k = 1; k < n; k <<= 1) {
        int p = 0;
        for (int i = n - k; i < n; i++) sa_tmp[p++] = i;
        for (int i = 0; i < n; i++) if (sa[i] >= k) sa_tmp[p++] = sa[i] - k;
        
        std::fill(count.begin(), count.end(), 0);
        for (int i = 0; i < n; i++) count[rank[i]]++;
        for (int i = 1; i <= max_val; i++) count[i] += count[i - 1];
        for (int i = n - 1; i >= 0; i--) sa[--count[rank[sa_tmp[i]]]] = sa_tmp[i];
        
        next_rank[sa[0]] = 0;
        int r = 0;
        for (int i = 1; i < n; i++) {
            int i1 = sa[i], i2 = sa[i - 1];
            bool same = (rank[i1] == rank[i2] && 
                         (i1 + k < n ? rank[i1 + k] : -1) == (i2 + k < n ? rank[i2 + k] : -1));
            next_rank[i1] = same ? r : ++r;
        }
        rank = next_rank;
        if (r == n - 1) break;
    }
    return sa;
}

std::vector<int> suffix_array(const std::string& s) {
    int n = int(s.size());
    std::vector<int> s2(n);
    for (int i = 0; i < n; i++) s2[i] = (unsigned char)s[i];
    return suffix_array(s2);
}

// Kasai's algorithm による O(N) の LCP Array 構築
template <class T>
std::vector<int> lcp_array(const std::vector<T>& s, const std::vector<int>& sa) {
    int n = int(s.size());
    if (n <= 1) return {};
    std::vector<int> rnk(n);
    for (int i = 0; i < n; i++) {
        rnk[sa[i]] = i;
    }
    std::vector<int> lcp(n - 1);
    int h = 0;
    for (int i = 0; i < n; i++) {
        if (rnk[i] == 0) continue;
        int j = sa[rnk[i] - 1];
        if (h > 0) h--;
        while (i + h < n && j + h < n && s[i + h] == s[j + h]) h++;
        lcp[rnk[i] - 1] = h;
    }
    return lcp;
}

std::vector<int> lcp_array(const std::string& s, const std::vector<int>& sa) {
    int n = int(s.size());
    std::vector<int> s2(n);
    for (int i = 0; i < n; i++) s2[i] = (unsigned char)s[i];
    return lcp_array(s2, sa);
}