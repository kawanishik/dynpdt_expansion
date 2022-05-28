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
    std::vector<std::string> restore_keys_order;                                // 復元した文字列順を保存
    std::vector<std::pair<std::string, uint64_t>> not_leaf_node_compare_string; // 葉ノード以外の文字列を復元する際に比較した文字列とprefixの長さ
    map.test_CPD_order_and_insert_new_dic(new_map, children, node_id, cnt_leaf_per_node, 0, store_string, restore_keys_order, not_leaf_node_compare_string);
    
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

    std::vector<std::string> restore_keys_order_exp = {
        "aaaaaaccc",
        "aaaaaabbb",
        "aaabbbccc",
        "aaabbbbbb",
        "bbbccccccddd",
        "bbbcccccc",
        "bbbbbbbbb",
        "aaaaaaaaa"
    };

    std::vector<std::pair<std::string, uint64_t>> not_leaf_node_compare_string_exp {
        {"aaabbbbbb", 4},
        {"bbbcccccc", 4},
        {"bbbbbbbbb", 1}
    };

    // 文字列の復元順のテスト
    bool check_restore_key_order = true;
    if(restore_keys_order.size() != restore_keys_order_exp.size()) check_restore_key_order = false;
    if(check_restore_key_order) {
        for(int i=0; i < restore_keys_order_exp.size(); i++) {
            if(restore_keys_order[i] != restore_keys_order_exp[i]) check_restore_key_order = false;
        }
    }
    
    // prefixを使用して復元した文字列、抽出に使用した文字列と長さのテスト
    bool check_not_leaf_node_compare_string = true;
    if(not_leaf_node_compare_string.size() != not_leaf_node_compare_string_exp.size()) check_not_leaf_node_compare_string = false;
    if(check_not_leaf_node_compare_string) {
        for(int i=0; i < not_leaf_node_compare_string_exp.size(); i++) {
            if(not_leaf_node_compare_string[i].first != not_leaf_node_compare_string_exp[i].first) check_not_leaf_node_compare_string = false;
            if(not_leaf_node_compare_string[i].second != not_leaf_node_compare_string_exp[i].second) check_not_leaf_node_compare_string = false;
        }
    }

    // テスト結果の出力
    std::cout << "文字列の復元順のテスト : " << (check_restore_key_order ? "ok." : "failed.") << std::endl;
    std::cout << "部分文字列のテスト     : " << (check_not_leaf_node_compare_string ? "ok." : "failed.") << std::endl;
}