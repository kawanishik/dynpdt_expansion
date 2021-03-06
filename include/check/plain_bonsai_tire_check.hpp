#ifndef POPLAR_TRIE_PLAIN_BONSAI_TRIE_CHECK_HPP
#define POPLAR_TRIE_PLAIN_BONSAI_TRIE_CHECK_HPP

#include "../poplar/bit_tools.hpp"
#include "../poplar/bit_vector.hpp"
#include "../poplar/compact_vector.hpp"
#include "../poplar/hash.hpp"

#include <queue>
#include <algorithm>

namespace poplar {

static double node_transition_search_time = 0.0;

template <uint32_t MaxFactor = 90, typename Hasher = hash::vigna_hasher>
class plain_bonsai_trie_check {
  private:
    static_assert(0 < MaxFactor and MaxFactor < 100);

  public:
    static constexpr uint64_t nil_id = UINT64_MAX;
    static constexpr uint32_t min_capa_bits = 16;

    static constexpr auto trie_type_id = trie_type_ids::BONSAI_TRIE;

  public:
    plain_bonsai_trie_check() = default;

    plain_bonsai_trie_check(uint32_t capa_bits, uint32_t symb_bits) {
        capa_size_ = size_p2{std::max(min_capa_bits, capa_bits)};
        symb_size_ = size_p2{symb_bits};
        max_size_ = static_cast<uint64_t>(capa_size_.size() * MaxFactor / 100.0);
        table_ = compact_vector{capa_size_.size(), capa_size_.bits() + symb_size_.bits()};
    }

    ~plain_bonsai_trie_check() = default;

    uint64_t get_root() const {
        assert(size_ != 0);
        return 1;
    }

    void add_root() {
        assert(size_ == 0);
        size_ = 1;
    }

    uint64_t find_child(uint64_t node_id, uint64_t symb) const {
        assert(node_id < capa_size_.size());
        assert(symb < symb_size_.size());

        // Stopwatch sw;

        if (size_ == 0) {
            return nil_id;
        }

        uint64_t key = make_key_(node_id, symb);
        assert(key != 0);

        // int cnt = 0;
        for (uint64_t i = Hasher::hash(key) & capa_size_.mask();; i = right_(i)) {
            // cnt += 1;
            if (i == 0) {
                // table_[0] is always empty so that table_[i] = 0 indicates to be empty.
                continue;
            }
            if (i == get_root()) {
                continue;
            }
            if (table_[i] == 0) {
                // encounter an empty slot
                // node_transition_search_time += sw.get_micro_sec();
                return nil_id;
            }
            if (table_[i] == key) {
                // if(symb == 255) { // ?????????????????????????????????????????????????????????
                    // cnt_linear_proving[cnt] += 1;
                    // cnt_linear_proving_all += 1;
                // }
                // node_transition_search_time += sw.get_micro_sec();
                return i;
            }
        }
    }

    bool add_child(uint64_t& node_id, uint64_t symb) {
        assert(node_id < capa_size_.size());
        assert(symb < symb_size_.size());

        uint64_t key = make_key_(node_id, symb);
        assert(key != 0);

        for (uint64_t i = Hasher::hash(key) & capa_size_.mask();; i = right_(i)) {
            if (i == 0) {
                // table_[0] is always empty so that any table_[i] = 0 indicates to be empty.
                continue;
            }

            if (i == get_root()) {
                continue;
            }

            if (table_[i] == 0) {
                // this slot is empty
                if (size_ == max_size_) {
                    return false;  // needs to expand
                }

                table_.set(i, key);

                ++size_;
                node_id = i;

                return true;
            }

            if (table_[i] == key) {
                node_id = i;
                return false;  // already stored
            }
        }
    }

    // ????????????????????????????????????????????????
    bool add_child_new_table(uint64_t& node_id, uint64_t symb) {
        // std::cout << "symb : " << symb << std::endl;
        // std::cout << "node_id : " << node_id << std::endl;
        assert(node_id < capa_size_.size());
        assert(symb < symb_size_.size());

        uint64_t key = make_key_(node_id, symb); // node???symb????????????????????????node_id?????????????????????ID
        assert(key != 0);

        // int cnt = 0;
        for (uint64_t i = Hasher::hash(key) & capa_size_.mask();; i = right_(i)) { // i?????????????????????????????????????????????i=(i+1) mod mask???hash?????????????????????
            // cnt += 1;
            if (i == 0) {
                // table_[0] is always empty so that any table_[i] = 0 indicates to be empty.
                continue;
            }

            if (i == get_root()) {
                continue;
            }

            if (new_table_[i] == 0) { // ?????????????????????
                // this slot is empty
                // if (size_ == max_size_) {
                //     //cnt_compare.clear();
                //     return false;  // needs to expand
                // }

                new_table_.set(i, key); // ????????????

                ++size_;
                node_id = i;

                //count(cnt);
                return true;
            }

            if (new_table_[i] == key) { // ??????????????????
                node_id = i;
                return false;  // already stored
            }
        }
    }

    // ?????????????????????(lower_bound)
    std::pair<bool, uint64_t> BinarySearch(const std::vector<info_fp>& data, uint64_t match) {
        uint64_t size = data.size();
        if(size == 0) return {false, 0};
        int64_t left = 0;
        int64_t right = size - 1;
        int64_t mid;
        while(left <= right) {
            mid = (left + right) / 2;
            if(data[mid].match == match) {
                return {true, mid};
            } else if(data[mid].match > match) {
                right = mid - 1;
            } else {
                left = mid + 1;
            }
        }
        return {false, uint64_t(left)};
    }

    // CP??????????????????????????????????????????????????????
    std::tuple<std::vector<std::vector<info_fp>>, std::vector<uint64_t>, std::vector<bool>> return_partial_CP_info(const std::array<uint8_t, 256>& restore_codes_) {
        // std::cout << "--- return_partial_CP_info ---" << std::endl;
        uint64_t table_size = table_.size();
        std::vector<parent_info> parent(table_size);
        std::vector<uint64_t> partial_num(table_size, 0); // ???????????????????????????????????????(?????????????????????)
        
        std::vector<std::vector<info_fp>> fork_pos(table_size);// std::vector<std::vector<std::pair<uint64_t, uint64_t>>> fork_pos(table_size); // ???????????????????????????????????????????????????????????????(????????????)
        std::vector<uint64_t> cnt_leaf(table_size, 0); // ?????????????????????????????????????????????????????????????????????????????????
        std::vector<uint64_t> all_branch(table_size, 0); // ??????????????????????????????????????????????????????????????????

        // O(n)??????????????????????????????(??????????????????)????????????
        // ??????????????????????????????????????????????????????????????????????????????
        for(uint64_t i=0; i < table_size; i++) {
            if(table_[i] != 0) {
                auto [p, label] = get_parent_and_symb(i); // ???????????????????????????
                auto [c, match] = std::pair{uint8_t(restore_codes_[label % 256]), label/256}; // ????????????????????????????????????
                partial_num[i] += 1; // ?????????????????????????????????
                partial_num[p] += 1; // ???????????????????????????????????????
                // parent[i].first = p;
                // parent[i].second = match;
                parent[i].parent = p;
                parent[i].match = match;
                parent[i].c = c;
            }
        }

        // queue????????????????????????????????????????????????????????????(CP????????????)
        std::queue<uint64_t> que;
        std::vector<bool> check_bottom(table_size, false);
        // ????????????????????????????????????
        for(uint64_t i=0; i < table_size; i++) {
            if(partial_num[i] == 1) {
                que.push(i);
                check_bottom[i] = true;
            }
        }
        check_bottom[get_root()] = true; // ????????????????????????????????????????????????????????????
        while(!que.empty()) { // ???????????????????????????
            uint64_t q = que.front();
            que.pop();
            auto [node_id, match, c] = parent[q];
            if(c != 0) cnt_leaf[q] += 1; // ??????????????????????????????????????????
            uint64_t p = node_id;
            cnt_leaf[p] += cnt_leaf[q];
            if(match != 0) all_branch[p] += cnt_leaf[q];    // match=0????????????????????????????????????????????????????????????????????????
            // fork_pos???????????????????????????????????????????????????????????????????????????????????????
            // ??????????????????????????????????????????????????????????????????????????????
            auto [flag, pos] = BinarySearch(fork_pos[p], match);    // match?????????????????????????????????????????????????????????
            if(flag) {
                fork_pos[p][pos].cnt += cnt_leaf[q];
                fork_pos[p][pos].children.push_back({q, cnt_leaf[q]});
            } else {
                fork_pos[p].insert(fork_pos[p].begin()+pos, info_fp{match, cnt_leaf[q]});
                fork_pos[p][pos].children.push_back({q, cnt_leaf[q]});
            }

            partial_num[p]--;
            if(partial_num[p] == 1) { // ????????????????????????????????????????????????
                que.push(p);
            }
        }

        return std::tuple{fork_pos, all_branch, check_bottom};
    }

    std::pair<uint64_t, uint64_t> get_parent_and_symb(uint64_t node_id) const {
        assert(node_id < capa_size_.size());

        uint64_t key = table_[node_id];
        if (key == 0) {
            // root or not exist
            return {nil_id, 0};
        }
        // Returns pair (parent, label)
        return std::make_pair(key >> symb_size_.bits(), key & symb_size_.mask());
    };

    class node_map {
      public:
        node_map() = default;

        node_map(compact_vector&& map, bit_vector&& done_flags)
            : map_{std::move(map)}, done_flags_{std::move(done_flags)} {}

        ~node_map() = default;

        uint64_t operator[](uint64_t i) const {
            return done_flags_[i] ? map_[i] : UINT64_MAX;
        }

        uint64_t size() const {
            return map_.size();
        }

        node_map(const node_map&) = delete;
        node_map& operator=(const node_map&) = delete;

        node_map(node_map&& rhs) noexcept = default;
        node_map& operator=(node_map&& rhs) noexcept = default;

      private:
        compact_vector map_;
        bit_vector done_flags_;
    };

    bool needs_to_expand() const {
        return max_size() <= size();
    }

    node_map expand() {
        plain_bonsai_trie_check new_ht{capa_bits() + 1, symb_size_.bits()};
        new_ht.add_root();

#ifdef POPLAR_EXTRA_STATS
        new_ht.num_resize_ = num_resize_ + 1;
#endif

        bit_vector done_flags(capa_size());
        done_flags.set(get_root());

        table_.set(get_root(), new_ht.get_root());

        std::vector<std::pair<uint64_t, uint64_t>> path;
        path.reserve(256);

        // 0 is empty, 1 is root
        for (uint64_t i = 2; i < table_.size(); ++i) {
            if (done_flags[i] || table_[i] == 0) {
                // skip already processed or empty elements
                continue;
            }

            path.clear();
            uint64_t node_id = i;

            do {
                auto [parent, label] = get_parent_and_symb(node_id);
                assert(parent != nil_id);
                path.emplace_back(std::make_pair(node_id, label));
                node_id = parent;
            } while (!done_flags[node_id]);

            uint64_t new_node_id = table_[node_id];

            for (auto rit = std::rbegin(path); rit != std::rend(path); ++rit) {
                new_ht.add_child(new_node_id, rit->second);
                table_.set(rit->first, new_node_id);
                done_flags.set(rit->first);
            }
        }

        node_map node_map{std::move(table_), std::move(done_flags)};
        std::swap(*this, new_ht);

        return node_map;
    }

    bool is_use_table(uint64_t i) {
        if(table_[i] == 0) return false;
        return true;
    }

    void reset_data_() {
        size_ = 0;
    }

    void expand_tmp_table() {
        reset_data_();
        // capa_size_ = size_p2{capa_bits()+1}; // ???????????????????????????????????????????????????expand?????????????????????????????????
        capa_size_ = size_p2{capa_bits()};
        symb_size_ = size_p2{symb_size_.bits()};
        max_size_ = static_cast<uint64_t>(capa_size_.size() * MaxFactor / 100.0); // ?????????????????????????????????
        new_table_ = compact_vector{capa_size_.size(), capa_size_.bits() + symb_size_.bits()};
    }

    void expand_restore_string() {
        reset_data_();
        plain_bonsai_trie_check new_ht{capa_bits(), symb_size_.bits()};
        std::swap(*this, new_ht);
    }

    void move_table() {
        std::swap(table_, new_table_);
        new_table_ = {};
    }

    void set_first_insert(bool flag) {
        new_table_first_insert = flag;
    }

    bool checK_first_insert() {
        if(new_table_first_insert) return true;
        return false;
    }

    void reset_cnt_compare() {
        std::cout << "--- reset_cnt_linear_proving ---" << std::endl;
        cnt_linear_proving.clear();
        cnt_linear_proving_all = 0;
    }

    double get_node_transition_search_time() {
        return node_transition_search_time;
    }

    void show_cnt_compare() {
        std::cout << "--- cnt_linear_proving ---" << std::endl;
        std::cout << "all : " << cnt_linear_proving_all << std::endl;
        uint64_t all = 0;
        for(auto p : cnt_linear_proving) {
            // std::cout << p.first << " : " << (long double)(p.second) / (long double)(cnt_linear_proving_all) << std::endl;
            // std::cout << p.first << " : " << p.second << std::endl;
            std::cout << p.first << std::endl;
            all += p.second;
        }
        std::cout << "---" << std::endl;
        for(auto p : cnt_linear_proving) {
            std::cout << p.second << std::endl;
        }
        std::cout << "all : " << all << std::endl;
    }

    // # of registerd nodes
    uint64_t size() const {
        return size_;
    }
    uint64_t max_size() const {
        return max_size_;
    }
    uint64_t capa_size() const {
        return capa_size_.size();
    }
    uint32_t capa_bits() const {
        return capa_size_.bits();
    }
    uint64_t symb_size() const {
        return symb_size_.size();
    }
    uint32_t symb_bits() const {
        return symb_size_.bits();
    }
#ifdef POPLAR_EXTRA_STATS
    uint64_t num_resize() const {
        return num_resize_;
    }
#endif
    uint64_t alloc_bytes() const {
        return table_.alloc_bytes();
    }

    void show_stats(std::ostream& os, int n = 0) const {
        auto indent = get_indent(n);
        show_stat(os, indent, "name", "plain_bonsai_trie_check");
        show_stat(os, indent, "factor", double(size()) / capa_size() * 100);
        show_stat(os, indent, "max_factor", MaxFactor);
        show_stat(os, indent, "size", size());
        show_stat(os, indent, "alloc_bytes", alloc_bytes());
        show_stat(os, indent, "capa_bits", capa_bits());
        show_stat(os, indent, "symb_bits", symb_bits());
#ifdef POPLAR_EXTRA_STATS
        show_stat(os, indent, "num_resize", num_resize_);
#endif
    }

    plain_bonsai_trie_check(const plain_bonsai_trie_check&) = delete;
    plain_bonsai_trie_check& operator=(const plain_bonsai_trie_check&) = delete;

    plain_bonsai_trie_check(plain_bonsai_trie_check&&) noexcept = default;
    plain_bonsai_trie_check& operator=(plain_bonsai_trie_check&&) noexcept = default;

  private:
    compact_vector table_;
    compact_vector new_table_;
    bool new_table_first_insert = false;
    uint64_t size_ = 0;  // # of registered nodes
    uint64_t max_size_ = 0;  // MaxFactor% of the capacity
    size_p2 capa_size_;
    size_p2 symb_size_;
#ifdef POPLAR_EXTRA_STATS
    uint64_t num_resize_ = 0;
#endif

    uint64_t make_key_(uint64_t node_id, uint64_t symb) const {
        return (node_id << symb_size_.bits()) | symb;
    }
    uint64_t right_(uint64_t slot_id) const {
        return (slot_id + 1) & capa_size_.mask();
    }
};

}  // namespace poplar

#endif  // POPLAR_TRIE_PLAIN_BONSAI_TRIE_CHECK_HPP
