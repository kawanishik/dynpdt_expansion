#ifndef POPLAR_TRIE_PLAIN_BONSAI_NLM_SNR_HPP
#define POPLAR_TRIE_PLAIN_BONSAI_NLM_SNR_HPP

#include <memory>
#include <vector>

#include "../poplar/basics.hpp"
#include "../poplar/compact_vector.hpp"

namespace poplar {

template <typename Value>
class plain_bonsai_nlm_SNR {
  public:
    using value_type = Value;

    static constexpr auto trie_type_id = trie_type_ids::BONSAI_TRIE;

  public:
    plain_bonsai_nlm_SNR() = default;

    explicit plain_bonsai_nlm_SNR(uint32_t capa_bits) : ptrs_(1ULL << capa_bits) {}

    ~plain_bonsai_nlm_SNR() = default;

    std::pair<const value_type*, uint64_t> compare(uint64_t pos, const char_range& key) const {
        assert(pos < ptrs_.size());
        assert(ptrs_[pos]);

        const uint8_t* ptr = ptrs_[pos].get();

        if (key.empty()) {
            return {reinterpret_cast<const value_type*>(ptr), 0};
        }

        for (uint64_t i = 0; i < key.length(); ++i) {
            if (key[i] != ptr[i]) {
                return {nullptr, i};
            }
        }

        return {reinterpret_cast<const value_type*>(ptr + key.length()), key.length()};
    }

    value_type* insert(uint64_t pos, const char_range& key) {
        assert(!ptrs_[pos]);

        ++size_;

        uint64_t length = key.length();
        ptrs_[pos] = std::make_unique<uint8_t[]>(length + sizeof(value_type));
        auto ptr = ptrs_[pos].get();
        copy_bytes(ptr, key.begin, length);

        label_bytes_ += length + sizeof(value_type);

#ifdef POPLAR_EXTRA_STATS
        max_length_ = std::max(max_length_, length);
        sum_length_ += length;
#endif

        auto ret = reinterpret_cast<value_type*>(ptr + length);
        *ret = static_cast<value_type>(0);

        return ret;
    }

    template <typename T>
    void expand(const T& pos_map) {
        std::vector<std::unique_ptr<uint8_t[]>> new_ptrs(ptrs_.size() * 2);
        for (uint64_t i = 0; i < pos_map.size(); ++i) {
            if (pos_map[i] != UINT64_MAX) {
                new_ptrs[pos_map[i]] = std::move(ptrs_[i]);
            }
        }
        ptrs_ = std::move(new_ptrs);
    }

    // ??????ptrs_??????new_map.ptrs????????????????????????
    void move_original_to_new_ptrs_(std::vector<std::unique_ptr<uint8_t[]>>& original_ptrs_, std::vector<std::pair<uint64_t, uint64_t>>& move_pair) {
        // std::cout << "--- move_original_to_new_ptrs ---" << std::endl;
        // std::cout << "pair_size : " << move_pair.size() << std::endl;
        for(auto [pre_node_id, new_node_id] : move_pair) {
            ptrs_[new_node_id] = std::move(original_ptrs_[pre_node_id]);
        }
    }

    std::vector<std::unique_ptr<uint8_t[]>> &return_ptrs_() {
        return ptrs_;
    }

    uint8_t* return_string_pointer(uint64_t pos) const {
        return ptrs_[pos].get();
    }

    uint64_t size() const {
        return size_;
    }
    uint64_t num_ptrs() const {
        return ptrs_.size();
    }
    uint64_t alloc_bytes() const {
        uint64_t bytes = 0;
        bytes += ptrs_.capacity() * sizeof(std::unique_ptr<uint8_t[]>);
        bytes += label_bytes_;
        return bytes;
    }

    void show_stats(std::ostream& os, int n = 0) const {
        auto indent = get_indent(n);
        show_stat(os, indent, "name", "plain_bonsai_nlm_SNR");
        show_stat(os, indent, "size", size());
        show_stat(os, indent, "num_ptrs", num_ptrs());
        show_stat(os, indent, "alloc_bytes", alloc_bytes());
#ifdef POPLAR_EXTRA_STATS
        show_stat(os, indent, "max_length", max_length_);
        show_stat(os, indent, "ave_length", double(sum_length_) / size());
#endif
    }

    plain_bonsai_nlm_SNR(const plain_bonsai_nlm_SNR&) = delete;
    plain_bonsai_nlm_SNR& operator=(const plain_bonsai_nlm_SNR&) = delete;

    plain_bonsai_nlm_SNR(plain_bonsai_nlm_SNR&&) noexcept = default;
    plain_bonsai_nlm_SNR& operator=(plain_bonsai_nlm_SNR&&) noexcept = default;

  private:
    std::vector<std::unique_ptr<uint8_t[]>> ptrs_;
    uint64_t size_ = 0;
    uint64_t label_bytes_ = 0;
#ifdef POPLAR_EXTRA_STATS
    uint64_t max_length_ = 0;
    uint64_t sum_length_ = 0;
#endif
};

}  // namespace poplar

#endif  // POPLAR_TRIE_PLAIN_BONSAI_NLM_SNR_HPP
