/*
再帰関数を使用せずに作成した関数をテストするためのファイル
再帰関数を使用して求めたものと比較する
*/

#include "map_dr.hpp"
#include "plain_bonsa_trie_dr.hpp"
#include "plain_bonsai_nlm_dr.hpp"
#include <cstdio>
#include <vector>
#include <string>
#include <stack>

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

    map.compute_node_connect_and_blanch_num(children, cnt_leaf_per_node);// ノード関係の取得

    // 再帰関数を使用しないCPD順の求め方
    std::stack<uint64_t> cpd_ord = map.require_centroid_path_order_not_using_recursion(children, cnt_leaf_per_node);

    std::vector<std::string> restore_keys_order_exp = {
        "aaaaaabbb",
        "aaaaaaccc",
        "aaabbbccc",
        "aaabbbbbb",
        "bbbccccccddd",
        "bbbcccccc",
        "bbbbbbbbb",
        "aaaaaaaaa"
    };

    bool check_restore_key_order = true;
    uint64_t i = 0;
    if(cpd_ord.size() != restore_keys_order_exp.size()) check_restore_key_order = false;
    while(!cpd_ord.empty()) {
        uint64_t node_id = cpd_ord.top();
        cpd_ord.pop();
        std::string restore_key = map.restore_insert_string(node_id);
        if(restore_key != restore_keys_order_exp[i]) {
            check_restore_key_order = false;
            break;
        }
        i++;
    }

    std::cout << "文字列の復元順のテスト : " << (check_restore_key_order ? "ok." : "failed.") << std::endl;

    return 0;
}