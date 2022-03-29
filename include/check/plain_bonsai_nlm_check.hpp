#ifndef POPLAR_TRIE_PLAIN_BONSAI_NLM_CHECK_HPP
#define POPLAR_TRIE_PLAIN_BONSAI_NLM_CHECK_HPP

#include <memory>
#include <vector>

#include "../poplar/basics.hpp"
#include "../poplar/compact_vector.hpp"

namespace poplar {

template <typename Value>
class plain_bonsai_nlm_check {
  public:
    using value_type = Value;

    static constexpr auto trie_type_id = trie_type_ids::BONSAI_TRIE;

  public:
    plain_bonsai_nlm_check() = default;

    explicit plain_bonsai_nlm_check(uint32_t capa_bits) : ptrs_(1ULL << capa_bits) {}

    ~plain_bonsai_nlm_check() = default;

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

    // new_ptrsに対する文字列比較
    std::pair<const value_type*, uint64_t> compare_new_ptrs(uint64_t pos, const char_range& key, const uint64_t start) const {
        // std::cout << " --- compare_new_ptrs ---" << std::endl;
        assert(pos < new_ptrs_.size());
        assert(new_ptrs_[pos]);

        const uint8_t* ptr = new_ptrs_[pos].get();

        if (key.empty()) {
            return {reinterpret_cast<const value_type*>(ptr), 0};
        }

        for (uint64_t i = 0; i < key.length(); ++i) {
            //std::cout << key[i] << ", " << ptr[i] << ", " << i << std::endl;
            if (key[i] != ptr[start+i]) {
                //std::cout << key[i] << ", " << ptr[i] << std::endl;
                return {nullptr, start+i}; // iは間違えた箇所
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

    value_type* insert_new_table(uint64_t pos, const char_range& key) {
        //std::cout << "--- insert ---" << std::endl;
        assert(!new_ptrs_[pos]);
        ++size_;

        uint64_t length = key.length();
        // std::cout << "length : " << length << std::endl;
        // std::cout << "pos : " << pos << std::endl;
        new_ptrs_[pos] = std::make_unique<uint8_t[]>(length + sizeof(value_type));
        auto ptr = new_ptrs_[pos].get();    // unique_ptr::get() : 保持しているポインタを返す
        copy_bytes(ptr, key.begin, length); // ptrの先頭からにkeyを長さlength分，コピーする

// #ifdef POPLAR_EXTRA_STATS
//         max_length_ = std::max(max_length_, length);
//         sum_length_ += length;
// #endif

        auto ret = reinterpret_cast<value_type*>(ptr + length);
        *ret = static_cast<value_type>(0);

        return ret;
    }

    void reset_data_() {
        size_ = 0;
    }

    uint8_t* return_string(uint64_t pos) const {
        return ptrs_[pos].get();
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

    void expand_tmp_ptrs() {
        new_ptrs_.resize(ptrs_.size());
    }

    void move_ptrs() {
        std::swap(ptrs_, new_ptrs_);
        new_ptrs_ = decltype(new_ptrs_)();
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
        show_stat(os, indent, "name", "plain_bonsai_nlm_check");
        show_stat(os, indent, "size", size());
        show_stat(os, indent, "num_ptrs", num_ptrs());
        show_stat(os, indent, "alloc_bytes", alloc_bytes());
#ifdef POPLAR_EXTRA_STATS
        show_stat(os, indent, "max_length", max_length_);
        show_stat(os, indent, "ave_length", double(sum_length_) / size());
#endif
    }

    plain_bonsai_nlm_check(const plain_bonsai_nlm_check&) = delete;
    plain_bonsai_nlm_check& operator=(const plain_bonsai_nlm_check&) = delete;

    plain_bonsai_nlm_check(plain_bonsai_nlm_check&&) noexcept = default;
    plain_bonsai_nlm_check& operator=(plain_bonsai_nlm_check&&) noexcept = default;

  private:
    std::vector<std::unique_ptr<uint8_t[]>> ptrs_;
    std::vector<std::unique_ptr<uint8_t[]>> new_ptrs_;
    uint64_t size_ = 0;
    uint64_t label_bytes_ = 0;
#ifdef POPLAR_EXTRA_STATS
    uint64_t max_length_ = 0;
    uint64_t sum_length_ = 0;
#endif
};

}  // namespace poplar

#endif  // POPLAR_TRIE_PLAIN_BONSAI_NLM_CHECK_HPP
