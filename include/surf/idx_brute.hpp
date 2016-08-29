#pragma once

#include <algorithm>
#include <limits>
#include <queue>
#include <set>

namespace surf {

/*! Trivial implementation of top-k retrieval using brute force search, for
 * testing.
 */
template <typename t_alphabet_tag>
class idx_brute {
public:
    using alphabet_category = t_alphabet_tag;
    using index_category = int;

    void load(sdsl::cache_config& cc) {
        // TODO implement
    }

    using t_result_set = std::vector<std::pair<uint64_t, double>>;

    class top_k_iterator {
    private:
        size_t m_index = 0;
        const t_result_set& m_results;
    public:
        top_k_iterator() = delete;
        top_k_iterator(const top_k_iterator&) = default;
        top_k_iterator(top_k_iterator&&) = default;
        top_k_iterator& operator=(const top_k_iterator&) = default;
        top_k_iterator& operator=(top_k_iterator&&) = default;

        explicit top_k_iterator(const t_result_set& results) : m_results(results) {}

        t_result_set::value_type operator*() const {
            return m_results[m_index];
        }

        bool done() const {
            return m_index >= m_results.size();
        }

        top_k_iterator& operator++() {
            m_index++;
            return *this;
        }

        auto extract_snippet(const size_t k) const {
            return string("foobar");
        }
    };

    template<typename t_pat_iter>
    top_k_iterator topk(t_pat_iter begin, t_pat_iter end,
                   bool multi_occ = false,
                   bool only_match = false) const {
        t_result_set results;
        return top_k_iterator(results);
    }

    void mem_info() const { }

    uint64_t doc_cnt() const{
        // TODO implement
        return 0;
    }

    uint64_t word_cnt() const {
        // TODO implement
        return 0;
    }
};

template <typename t_alphabet_type>
size_t size_in_bytes(const surf::idx_brute<t_alphabet_type>& idx) { return 0; }

template <typename t_alphabet_type>
void construct(idx_brute<t_alphabet_type>& idx,
            const std::string&,
            sdsl::cache_config& cc,
            uint8_t num_bytes)
{
    // TODO implement
}

}
