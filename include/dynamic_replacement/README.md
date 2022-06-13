# 動的に入れ替える手法を試すフォルダ

## それぞれのファイルの説明
1. map_dr.cpp, plain_bonsai_trie_dr.hpp, plain_bonsai_nlm_dr.hpp
    - それぞれの使用用途は、元のファイルと同様
1. testファイルについて
    - test_compute_node_connect_and_blanch_num.cpp
        - 関数compute_node_connect_and_blanch_numを確かめるためのファイル
    - test_restore_string.cpp
        - 文字列が、復元されているかを確認するためのファイル
        - どの部分をテストしているのかについて
            - 全文字列の復元をしているので、文字列順が想定通りであるのか
            - 文字列を復元する際にprefix部分を抜き出しているのだが、比較する文字列と長さが想定通りであるのか
    - test_CPD_order.cpp
        - 再帰関数を使用せずにCPD順を求めたものを再帰関数を使用したものと同じ結果になるのかを確かめたもの
        - compute_node_connect_and_blanch_num関数でノードの関係を調べる部分までは同じ