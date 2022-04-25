# 変更点

## dynamic_replacement関数
- 引数のvector<uint64_t> boxは、テストする際に使用(test_require_centroid_path_order_and_insert_dictionary)
    - 本番環境で実行する際には、このファイル内で使用しているboxはすべて、コメントアウトします

## compute_node_connect_and_blanch_num関数
- vector<uint64_t> blanch_num_expect_zeroの消去
    - 元々は、ゼロ分岐以外の累積和を計算するために使用していた
    - 累積和自体は、別の関数内で求めています

## require_centroid_path_order_and_insert_dictionary関数
- 変数の多くは、構造体にまとめました
- 累積和を計算する部分の追加
    - 0分岐以外の和をvariable.cumulative_sumに加算することにより実現
- map<> start_posについて
    - compute_node...関数を使用して、vector<> childrenは計算されていると想定
    - その際、図1のようなDynPDTは、図2のようなDynPDTとして格納しているので、children内の値は以下のようになっていて、lambda_(=32)を超える値が多々、存在しています
        - 図2のようにしている理由は、累積和などを求める際に、ダミーノードの個数をカウントしないようにするためです
    - そのため、今回は、vector<> start_pos(lambda_)のような形ではなく、map関数を使用しています

| ![Test1](img/origin.png) | ![Test1](img/children.png) |
| :---: | :---: |
| 図1 | 図2 |

- match_per_leaf_numをソートしている理由について
    - match_per_leaf_numの値が、下のようになっていると想定
    - ex. ( {0, 8}, {3, 7}, {13, 15} ) [first:分岐位置, second:葉の数]
    - 処理したい順としては、13 → 3 → 0分岐の順
    - これを実現するためにソートしています

- ソート以降の処理
    - for分で分岐後のノード数が多い順に処理していく
    - 初めに、分岐後のノードをchildrenより取得（同時に、子ノード以下の葉の数も取得）
    - 葉の数が多い順にソート
    - for分で多い順に処理

## HL分解について
- なんとなくの理解はできました
- しかし、どの部分に使用するのかが思いついていないです
- 現在は、新しい辞書に追加する際に、使用できるのではないかと考えている段階です

## 実行方法(それぞれのmapの挙動を確認する)
> git clone https://github.com/kawanishik/dynpdt_expansion.git  
> cd dynpdt_expansion  
> mkdir build && cd build  
> cmake ..  
> cd sample  
> make  
> ./sample ['マップ名']  

## 関数のテストを実行する際
> cd dynpdt_expansion  
> cd include/dynamic_replacement  
> mkdir build && cd build  
> cmake ..  
> make  
> ./test_compute_node_connect_and_blanch_num or ./test_require_centroid_path_order_and_insert_dictionary  