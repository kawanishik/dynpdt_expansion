#include "map_akr.hpp"
#include "plain_bonsai_trie_akr.hpp"
#include "plain_bonsai_nlm_akr.hpp"

#include <vector>
#include <string>
#include <cstdio>
#include <algorithm>

using namespace poplar;

int main() {
    std::vector<std::string> keys = {
        "aaaaaaaaa",
        "aaaaaabbb",
        "aaaaaaccc",
        "aaabbbbbb",
        "aaabbbccc",
        "bbbbbbbbb",
        "bbbcccccc",
        "bbbccccccddd"
    };

    map_akr<plain_bonsai_trie_akr<>, plain_bonsai_nlm_akr<int>> mp;
    for(auto it=keys.begin(); it != keys.end(); it++) {
        mp.update(*it);
    }

    std::vector<std::string> restore_keys = mp.all_key_restore();
    // std::vector<std::string> restore_keys = mp.all_key_restore_simple();
    std::sort(restore_keys.begin(), restore_keys.end());

    // 文字列が全て、復元できているのかを確認する
    if(restore_keys.size() != keys.size()) {
        std::cout << "size is different" << std::endl;
        return 0;
    }

    for(uint64_t i=0; i < keys.size(); i++) {
        if(keys[i] != restore_keys[i]) {
            std::cout << "failed to compare key" << std::endl;
            std::cout << "original_key: " << keys.size() << std::endl;
            std::cout << "restore_key : " << restore_keys[i] << std::endl;
            return 0;
        }
    }

    std::cout << "ok." << std::endl;

    return 0;
}