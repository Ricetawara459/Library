#define PROBLEM "https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=ALDS1_4_C"

#include <iostream>
#include <string>

#include "../../src/string/trie.hpp"

int main() {
    std::cin.tie(nullptr);
    std::ios::sync_with_stdio(false);

    int q;
    std::cin >> q;
    trie<> dict;
    while (q--) {
        std::string command, word;
        std::cin >> command >> word;
        if (command == "insert") {
            if (!dict.contains(word)) dict.insert(word);
        } else {
            std::cout << (dict.contains(word) ? "yes\n" : "no\n");
        }
    }
    return 0;
}
