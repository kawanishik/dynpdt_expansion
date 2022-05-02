/*
文字列を復元し、新しい辞書に追加する際に使用するテストファイル
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
    map.require_centroid_path_order_and_insert_dictionary_using_backward_cumulative(new_map, children, node_id, cnt_leaf_per_node, 0);
    

    std::cout << "ok." << std::endl;
}