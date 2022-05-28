/*
文字列を復元し、新しい辞書に追加する際に使用するテストファイル
確かめる用の関数を使用して、動作を確認している
*/

#include "map_dr.hpp"
#include "plain_bonsa_trie_dr.hpp"
#include "plain_bonsai_nlm_dr.hpp"
#include <cstdio>
#include <vector>
#include <string>

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
    map_dr<plain_bonsai_trie_dr<>, plain_bonsai_nlm_dr<int>> map(8);    // 元の辞書

    // 辞書に追加
    for(std::string& key : keys) {
        int* ptr = map.update(key);
        *ptr = 1;
    }

    std::vector<std::vector<std::pair<uint64_t, uint64_t>>> children;   // 子ノードの集合
    std::vector<uint64_t> cnt_leaf_per_node;                            // それぞれのノード以降の葉の数

    map.compute_node_connect_and_blanch_num(children, cnt_leaf_per_node);

    map_dr<plain_bonsai_trie_dr<>, plain_bonsai_nlm_dr<int>> new_map(8);
    uint64_t node_id = 1;   // 比較対象のノード番号
    std::string store_string = "";
    map.test_CPD_order_and_insert_new_dic(new_map, children, node_id, cnt_leaf_per_node, 0, store_string);
    
    bool test_check = true;
    for(int i=0; i < keys.size(); i++) {
        const int *ptr = map.find(keys[i]);
        if(not (ptr != nullptr and *ptr == 1)) {
            std::cout << "search_failed : " << i << "," << keys[i] << std::endl;
            test_check = false;
            return 0;
        }
    }

    std::cout << "serch_ok." << std::endl;

    std::vector<std::string> keys_exp = {
        "aaaaaaccc",
        "aaaaaabbb",
        "aaabbbccc",
        "aaabbbbbb",
        "bbbccccccddd",
        "bbbcccccc",
        "bbbbbbbbb",
        "aaaaaaaaa"
    };

    std::cout << "比較するキー" << std::endl;
    for(auto key : keys_exp) std::cout << key << std::endl;
}