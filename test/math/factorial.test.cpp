#define PROBLEM "https://yukicoder.me/problems/no/117"

#include <iostream>
#include <cstdio>
#include "../../src/math/modint.hpp"
#include "../../src/math/factorial.hpp"

using namespace std;

using mint = modint1000000007;

// グローバル領域に宣言して、型指定なしで使えるようにする
factorial<mint> fc;

int main() {
    int t;
    if (scanf("%d\n", &t) != 1) return 0;

    while (t--) {
        char type;
        int n, r;
        if (scanf("%c(%d,%d)\n", &type, &n, &r) != 3) break;

        if (type == 'C') {
            printf("%d\n", fc.comb(n, r).val());
        } else if (type == 'P') {
            printf("%d\n", fc.perm(n, r).val());
        } else if (type == 'H') {
            printf("%d\n", fc.homo(n, r).val());
        }
    }

    return 0;
}