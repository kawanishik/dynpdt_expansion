#include <iostream>
#include <map>
#include <tuple>
#include <chrono>
#include <fstream>
#include <random>
#include <poplar.hpp>
#include <algorithm>
#include <string_view>
#include <unistd.h>

namespace {

const int search_cnt = 10;
const int search_get_string = 1000000;

class Stopwatch {
  using hrc = std::chrono::high_resolution_clock;
  hrc::time_point start_;
 public:
  Stopwatch() : start_(hrc::now()) {}
  auto time_process() const {
    return hrc::now() - start_;
  }
  double get_sec() const {
    return std::chrono::duration<double>(time_process()).count();
  }
  double get_milli_sec() const {
    return std::chrono::duration<double, std::milli>(time_process()).count();
  }
  double get_micro_sec() const {
    return std::chrono::duration<double, std::micro>(time_process()).count();
  }
};

template <class Process>
double milli_sec_in(Process process) {
  Stopwatch sw;
  process();
  return sw.get_milli_sec();
}

inline uint64_t get_process_size() {
    FILE* fp = std::fopen("/proc/self/statm", "r");
    uint64_t dummy(0), vm(0);
    std::fscanf(fp, "%ld %ld ", &dummy, &vm);  // get resident (see procfs)
    std::fclose(fp);
    return vm * ::getpagesize();
}

std::vector<std::string> MakePartialTimeKeysets(const std::vector<std::string>& keys, uint64_t size) {
    std::vector<std::string> keysets(search_get_string);
    std::random_device rnd;
    for(int i=0; i < search_get_string; i++) {
        int v = rnd() % size;
        // keysets.push_back((keys[v]));
        keysets[i] = keys[v];
    }
    // std::sort(time_keysets[j].begin(), time_keysets[j].end());
    return keysets;
}

template <typename MAP>
double CalcSearchSpeed(const MAP &dyn_, const std::vector<std::string>& keys, uint64_t size) {
    std::cout << "--- CalcSearchSpeed ---" << std::endl;
    double time = 0.0;
    for(int i=0; i < search_cnt; i++) {
        std::vector<std::string> keysets = MakePartialTimeKeysets(keys, size);
        Stopwatch sw;
        for (int j = 0; j < search_get_string; j++) {
            // const int *ptr = dyn_.find(time_keysets[j]);
            const int *ptr = dyn_.find(keysets[j]);
            if (not (ptr != nullptr and *ptr == 1)) {
                // std::cerr << "Failed to find " << time_keysets[j] << std::endl;
                std::cerr << "Failed to find " << keysets[j] << std::endl;
                return -1;
            }
        }
        time += sw.get_milli_sec();
        // time += sw.get_sec();
    }
    return time / (double)(search_cnt);
}

std::pair<std::string, std::vector<std::pair<uint64_t, uint64_t>>> MakeSequentialKeysets(const std::vector<std::string>& keys, uint64_t size) {
    std::string sequential_key = "";
    std::vector<std::pair<uint64_t, uint64_t>> sequential_set(search_get_string);
    std::random_device rnd;
    for(int i=0; i < search_get_string; i++) {
        int v = rnd() % size;
        std::string key = keys[v];
        // std::cout << "key : " << key << std::endl;
        // sleep(1);
        sequential_set[i] = std::pair<uint64_t, uint64_t>{sequential_key.size(), key.size()};
        sequential_key += key;
        // std::cout << "sequential_key : " << sequential_key.size() << std::endl;
    }
    return {sequential_key, sequential_set};
}

// 検索するデータセットをシーケンシャルに配置
template <typename MAP>
double CalcSearchSpeedSequence(const MAP &dyn_, const std::vector<std::string>& keys, uint64_t size) {
    std::cout << "--- CalcSearchSpeedSequence ---" << std::endl;
    double time = 0.0;
    for(int i=0; i < search_cnt; i++) {
        auto [sequential_key, sequential_set] = MakeSequentialKeysets(keys, size);
        Stopwatch sw;
        for(int j=0; j < search_get_string; j++) {
            // std::cout << std::string_view{sequential_key}.substr(sequential_set[i].first, sequential_set[i].second) << std::endl;
            const int *ptr = dyn_.find(std::string_view{sequential_key}.substr(sequential_set[i].first, sequential_set[i].second));
            if (not (ptr != nullptr and *ptr == 1)) {
                std::cerr << "Failed to find " << std::string_view{sequential_key}.substr(sequential_set[i].first, sequential_set[i].second) << std::endl;
                return -1;
            }
        }
        time += sw.get_milli_sec();
    }
    return time / (double)(search_cnt);
}

template <typename MAP>
double AllDatasetSearchSpeed(const MAP &dyn_, const std::vector<std::string>& test_keys) {
    std::cout << "--- AllDatasetSearchSpeed ---" << std::endl;
    Stopwatch sw;
    for(int i=0; i < int(test_keys.size()); i++) {
        const int *ptr = dyn_.find(test_keys[i]);
        if (not (ptr != nullptr and *ptr == 1)) {
            std::cerr << "Failed to find " << test_keys[i] << std::endl;
            return -1;
        }
    }
    return sw.get_milli_sec();
}

template <typename MAP>
double CalcRandomFileSearchSpeed(const MAP &dyn_, const std::vector<std::vector<std::string>>& random_test_keys) {
    std::cout << "--- CalcRandomFileSearchSpeed ---" << std::endl;
    double time = 0.0;
    for(int i=0; i < search_cnt; i++) {
        Stopwatch sw;
        for(int j=0; j < search_get_string; j++) {
            const int *ptr = dyn_.find(random_test_keys[i][j]);
            if (not (ptr != nullptr and *ptr == 1)) {
                std::cerr << "Failed to find " << random_test_keys[i][j] << std::endl;
                return -1;
            }
        }
        time += sw.get_milli_sec();
    }
    return time / (double)(search_cnt);
}

void FileRead(std::vector<std::string>& keys, std::vector<std::string>& test_keys, std::vector<std::vector<std::string>>& random_test_keys) {
    std::string input_name = "../../../dataset/Titles-enwiki.txt";
    // std::string input_name = "../../../dataset/DS5";
    // std::string input_name = "../../../dataset/GeoNames.txt";
    // std::string input_name = "../../../dataset/AOL.txt";
    // std::string input_name = "../../../dataset/in-2004.txt";
    // std::string input_name = "../../../dataset/uk-2005.txt";
    // std::string input_name = "../../../dataset/webbase-2001.txt";

    // std::string input_name = "../../../dataset/enwiki_CP.txt";
    // std::string input_name = "../../../dataset/DS5_CP.txt";
    // std::string input_name = "../../../dataset/GeoNames_CP.txt";
    // std::string input_name = "../../../dataset/AOL_CP.txt";

    // std::string input_name = "../../../dataset/DS5_partial";
    // std::string input_name = "../../../dataset/tmp/enwiki_shuf.txt";
    // std::string input_name = "../../../dataset/tmp/DS5.txt";
    // std::string input_name = "../../../dataset/enwiki-20150205.line";
    // std::string input_name = "../../../dataset/wordnet-3.0-word";
    // std::string input_name = "../../../dataset/s9.txt"; // 家で実験するとき
    // std::string input_name = "../../../dataset/s14.txt"; // 大学で実験するとき

    std::ifstream ifs(input_name);
    if (!ifs) {
        std::cerr << "File not found input file: "<< input_name << std::endl;
        exit(0);
    }
    std::cout << "dataset : " << input_name.substr(17) << std::endl;
    for (std::string s; std::getline(ifs, s); ) {
        keys.push_back(s);
    }

    std::string test_name = "../../../dataset/test_enwiki.txt";
    // std::string test_name = "../../../dataset/test_DS5.txt";
    // std::string test_name = "../../../dataset/test_GeoNames.txt";
    // std::string test_name = "../../../dataset/test_AOL.txt";
    std::ifstream tfs(test_name);
    if (!tfs) {
        std::cerr << "File not found input file: "<< test_name << std::endl;
        exit(0);
    }
    std::cout << "dataset : " << test_name.substr(17) << std::endl;
    for (std::string s; std::getline(tfs, s); ) {
        test_keys.push_back(s);
    }

    // test_name = "../../../dataset/random_test_enwiki.txt";
    // // test_name = "../../../dataset/random_test_DS5.txt";
    // // test_name = "../../../dataset/random_test_GeoNames.txt";
    // // test_name = "../../../dataset/random_test_AOL.txt";
    // std::ifstream tfs2(test_name);
    // if (!tfs2) {
    //     std::cerr << "File not found input file: "<< test_name << std::endl;
    //     exit(0);
    // }
    // std::cout << "dataset : " << test_name.substr(17) << std::endl;
    // int cnt = 0;
    // int i = 0;
    // random_test_keys.resize(search_cnt);
    // for (std::string s; std::getline(tfs2, s); ) {
    //     random_test_keys[i].push_back(s);
    //     cnt++;
    //     if(cnt % search_get_string == 0) i++;
    // }
}
} // namespace

template<class It, class StringDictionary>
void insert_by_centroid_path_order(It begin, It end, int depth, StringDictionary& dict) {
    assert(end-begin > 0);
    if (end-begin == 1) {
        int* ptr = dict.update(*begin);
        *ptr = 1;
        return;
    }
    std::vector<std::tuple<It,It,char>> ranges;
    auto from = begin;
    auto to = begin;
    // leaf character is '\0'
    if (from->length() == depth) {
        ranges.emplace_back(from, ++to, '\0');
        from = to;
    }
    while (from != end) {
        assert(from->length() > depth);
        char c = (*from)[depth];
        while (to != end and (*to)[depth] == c) {
            ++to;
        }
        ranges.emplace_back(from, to, c);
        from = to;
    }
    std::sort(ranges.begin(), ranges.end(), [](auto l, auto r) {
        auto [fl,tl,cl] = l;
        auto [fr,tr,cr] = r;
        return tl-fl > tr-fr;
    });
    for (auto [f,t,c] : ranges) {
        insert_by_centroid_path_order(f, t, depth+1, dict);
    }
}

template<class It>
void require_by_centroid_path_order(It begin, It end, int depth, std::vector<std::string>& keys) {
    assert(end-begin > 0);
    if (end-begin == 1) {
        keys.push_back(*begin);
        return;
    }
    std::vector<std::tuple<It,It,char>> ranges;
    auto from = begin;
    auto to = begin;
    // leaf character is '\0'
    if (from->length() == depth) {
        ranges.emplace_back(from, ++to, '\0');
        from = to;
    }
    while (from != end) {
        assert(from->length() > depth);
        char c = (*from)[depth];
        while (to != end and (*to)[depth] == c) {
            ++to;
        }
        ranges.emplace_back(from, to, c);
        from = to;
    }
    std::sort(ranges.begin(), ranges.end(), [](auto l, auto r) {
        auto [fl,tl,cl] = l;
        auto [fr,tr,cr] = r;
        return tl-fl > tr-fr;
    });
    for (auto [f,t,c] : ranges) {
        require_by_centroid_path_order(f, t, depth+1, keys);
    }
}

void remake_dataset(std::vector<std::string>& keys) {
    std::sort(keys.begin(), keys.end());
    std::vector<std::string> tmp_keys = keys;
    keys.clear();
    require_by_centroid_path_order(tmp_keys.begin(), tmp_keys.end(), 0, keys);
    // for(int i=0; i < 10; i++) std::cout << keys[i] << std::endl;
}

template<class Map>
void bench(Map& map, std::vector<std::string>& keys, std::vector<std::string>& test_keys, std::vector<std::vector<std::string>>& random_test_keys) {
    // uint32_t capa_bits = 16;
    // uint64_t lambda = 32;
    // Map map{capa_bits, lambda};
    // Map map;

    // キー追加
    const auto input_num_keys = static_cast<int>(keys.size());
    std::cout << "input_num_keys : " << input_num_keys << std::endl;
    const auto test_num_keys = static_cast<int>(test_keys.size());
    std::cout << "test_num_keys : " << test_num_keys << std::endl;
    const auto random_test_num_keys = static_cast<int>(random_test_keys.size());
    std::cout << "random_test_num_keys : " << random_test_num_keys << std::endl;
    auto begin_size = get_process_size();

    Stopwatch sw;
    for(int i=0; i < input_num_keys; i++) {
        int* ptr = map.update(keys[i]);
        *ptr = 1;
    }

    auto ram_size = get_process_size() - begin_size;
    // auto time = sw.get_milli_sec();
    auto time = sw.get_sec();

    std::cout << "----------------------" << std::endl;
    // 追加されたキーが全て辞書に存在するのかをチェック
    // bool test_check = true;
    // map.reset_cnt_hash();
    // for(int i=0; i < input_num_keys; i++) {
    //     // std::cout << "*** key" << i << " : " << keys[i] << std::endl;
    //     const int *ptr = map.find(keys[i]);
    //     if(not (ptr != nullptr and *ptr == 1)) {
    //         std::cout << "search_failed : " << keys[i] << std::endl;
    //         test_check = false;
    //         return;
    //     }
    // }
    // std::cout << (test_check ? "ok." : "failed.") << std::endl;

    // std::cout << "Build time(m): " << time / 1000.0 << std::endl;
    std::cout << "Build time(m): " << time << std::endl;
    std::cout << "ram_size : " << ram_size << std::endl;
    std::cout << "capa_size : " << map.capa_size() << std::endl;

    // map.reset_cnt_hash();
    // auto search_time = AllDatasetSearchSpeed(map, test_keys);
    // auto search_time = CalcSearchSpeed(map, keys, input_num_keys);           // 実験に使用する測定方法
    // auto search_time = CalcRandomFileSearchSpeed(map, random_test_keys);
    // auto search_time = CalcSearchSpeedSequence(map, keys, input_num_keys); // string_viewで検索する際にエラー
    // map.show_cnt_hash();
    // std::cout << "time_search : " << search_time << std::endl;
}

void write_file(std::vector<uint64_t>& box) {
    std::ofstream of;
    std::string filename = "../../../result/result.txt";
    of.open(filename, std::ios::out);
    for(auto b : box) {
        of << b << std::endl;
    }
    of.close();
}

// 作成したCP入れ替えを複数回試す
template<class Map>
void multi_CP_swap(Map& map, std::vector<std::string>& keys, std::vector<std::string>& test_keys) {
    std::cout << "--- multi_CP_swap ---" << std::endl;
    // uint64_t input_num_keys = keys.size();
    // map.reset_cnt_hash();
    // double search_time = AllDatasetSearchSpeed(map, test_keys);
    // // search_time = CalcSearchSpeed(map, keys, input_num_keys);
    // map.show_cnt_hash();
    // std::cout << "time_search : " << search_time << std::endl;
    for(int i=0; i < 1; i++) { // 複数回呼び出して、どのようになるのかを検証する際に使用
        Stopwatch sw;
        // map.call_topo();
        // map.call_restore_string_CP(); // 全ての文字列を復元してCP順に並べる
        map.dynamic_replacement();
        double time = sw.get_milli_sec();
        std::cout << "time : " << time / 1000.0 << std::endl;
        if(i == -1) {
            map.reset_cnt_hash();
            double search_time = AllDatasetSearchSpeed(map, test_keys);
            // search_time = CalcSearchSpeed(map, keys, input_num_keys);
            map.show_cnt_hash();
            std::cout << "time_search : " << search_time << std::endl;
        }
    }
}

int main() {

    std::vector<std::string> keys;                              // 辞書構築用の配列
    std::vector<std::string> test_keys;                         // 全てのキーに対する検索時間を測定する際に使用
    std::vector<std::vector<std::string>> random_test_keys;     // ランダムに取り出したキーを検索する際に使用
    
    std::cout << "normal, initial_CPD, remake_CPD, check, check_CPD, dr" << std::endl;
    std::cout << "上記より選択してください : ";
    std::string input_name;
    std::cin >> input_name;

    
    FileRead(keys, test_keys, random_test_keys);

    if(input_name == "normal") { // 何も変更していない通常の辞書
        poplar::plain_bonsai_map<int> map;
        bench(map, keys, test_keys, random_test_keys);
    } else if(input_name == "remake_CPD") { // keysに入れ替えた文字列を一度保存してから、辞書に追加していく
        remake_dataset(keys);
        poplar::plain_bonsai_map<int> map;
        bench(map, keys, test_keys, random_test_keys);
    } else if(input_name == "check") {
        poplar::plain_bonsai_map_check<int> map;
        bench(map, keys, test_keys, random_test_keys);
        // map.call_topo();
        // multi_CP_swap(map, keys, test_keys);
    } else if(input_name == "check_CPD") {
        auto begin_size = get_process_size();
        poplar::plain_bonsai_map_check<int> map;
        Stopwatch sw;
        std::sort(keys.begin(), keys.end());
        insert_by_centroid_path_order(keys.begin(), keys.end(), 0, map);
        auto ram_size = get_process_size() - begin_size;
        auto time = sw.get_sec();
        std::cout << "Build time(m): " << time << std::endl;
        std::cout << "ram_size : " << ram_size << std::endl;
        std::cout << "capa_size : " << map.capa_size() << std::endl;
        map.call_topo();
    } else if(input_name == "initial_CPD") { // CPD順に入れ替えてから、辞書に格納
        const auto input_num_keys = static_cast<int>(keys.size());
        std::cout << "input_num_keys : " << input_num_keys << std::endl;
        const auto test_num_keys = static_cast<int>(test_keys.size());
        std::cout << "test_num_keys : " << test_num_keys << std::endl;
        const auto random_test_num_keys = static_cast<int>(random_test_keys.size());
        std::cout << "random_test_num_keys : " << random_test_num_keys << std::endl;

        auto begin_size = get_process_size();
        poplar::plain_bonsai_map<int> map;

        Stopwatch sw;
        std::sort(keys.begin(), keys.end());
        insert_by_centroid_path_order(keys.begin(), keys.end(), 0, map);
        auto ram_size = get_process_size() - begin_size;
        // auto time = sw.get_milli_sec();
        auto time = sw.get_sec();

        // bool test_check = true;
        // map.reset_cnt_hash();
        // for(int i=0; i < input_num_keys; i++) {
        //     // std::cout << "*** key" << i << " : " << keys[i] << std::endl;
        //     const int *ptr = map.find(keys[i]);
        //     if(not (ptr != nullptr and *ptr == 1)) {
        //         std::cout << "search_failed : " << keys[i] << std::endl;
        //         //std::cout << "ptr : " << *ptr << std::endl;
        //         test_check = false;
        //         // break;
        //         return 0;
        //     }
        // }
        // std::cout << (test_check ? "ok." : "failed.") << std::endl;

        // std::cout << "Build time(m): " << time / 1000.0 << std::endl;
        std::cout << "Build time(m): " << time << std::endl;
        std::cout << "ram_size : " << ram_size << std::endl;
        std::cout << "capa_size : " << map.capa_size() << std::endl;
        
        map.reset_cnt_hash();
        auto search_time = AllDatasetSearchSpeed(map, test_keys);
        // auto search_time = CalcSearchSpeed(map, keys, input_num_keys);
        // auto search_time = CalcRandomFileSearchSpeed(map, random_test_keys);
        map.show_cnt_hash();
        std::cout << "time_search : " << search_time << std::endl;
    } else if(input_name == "dr") { // 動的に辞書を入れ替える際に使用
        poplar::plain_bonsai_map_dr<int> map;
        bench(map, keys, test_keys, random_test_keys);
        multi_CP_swap(map, keys, test_keys);
    }else {
        std::cout << "そのようなデータセットは存在しません" << std::endl;
    }


    return 0;
}
