#pragma once

#include <algorithm>
#include <limits>
#include <queue>
#include <set>

#include "sdsl/int_vector.hpp"
#include "sdsl/sdsl_concepts.hpp"
#include "sdsl/suffix_trees.hpp"
#include "surf/util.hpp"

namespace surf {

/*! Trivial implementation of top-k retrieval using brute force search, for
 * testing.
 */
template <typename t_csa>
class idx_brute {
public:
    using alphabet_category = typename t_csa::alphabet_category;
    using csa_type = t_csa;
    using size_type = sdsl::int_vector<>::size_type;
    using result_set = std::vector<std::pair<uint64_t, double>>;

private:
    csa_type m_csa;
    sdsl::rrr_vector<> m_doc_splitters;
    sdsl::rrr_vector<>::rank_1_type m_doc_splitters_rank;
    result_set m_results;

public:
    void load(sdsl::cache_config& cc) {
        load_from_cache(m_csa, surf::KEY_CSA, cc, true);
        load_from_cache(m_doc_splitters, surf::KEY_DOCBORDER, cc, true);
        load_from_cache(m_doc_splitters_rank, surf::KEY_DOCBORDER_RANK, cc, true);
        m_doc_splitters_rank.set_vector(&m_doc_splitters);
    }

    class top_k_iterator {
    private:
        size_t m_index = 0;
        const result_set& m_results;
    public:
        top_k_iterator() = delete;
        top_k_iterator(const top_k_iterator&) = default;
        top_k_iterator(top_k_iterator&&) = default;
        top_k_iterator& operator=(const top_k_iterator&) = default;
        top_k_iterator& operator=(top_k_iterator&&) = default;

        explicit top_k_iterator(const result_set& results) : m_results(results) {}

        result_set::value_type operator*() const {
            return m_results[m_index];
        }

        bool done() const {
            return m_index >= m_results.size();
        }

        top_k_iterator& operator++() {
            m_index++;
            return *this;
        }

        auto extract_snippet(const size_type k) const {
            return string("foobar");
        }
    };

    template<typename t_pat_iter>
    top_k_iterator topk(t_pat_iter begin, t_pat_iter end,
                        bool multi_occ = false,
                        bool only_match = false) {
        assert(!only_match);
        auto occs = locate(m_csa, begin, end);
        std::sort(occs.begin(), occs.end());

        std::map<uint64_t, double> occs_by_doc;
        for (auto pos : occs) {
            auto doc = m_doc_splitters_rank(pos);
            occs_by_doc[doc] += 1;
        }

        m_results = result_set(occs_by_doc.begin(), occs_by_doc.end());
        std::sort(m_results.begin(), m_results.end(), [&](
                    const std::pair<uint64_t, double>& a,
                    const std::pair<uint64_t, double>& b) {
            return std::make_pair(a.second, a.first) >
                std::make_pair(b.second, b.first);
        });
        return top_k_iterator(m_results);
    }

    void mem_info() const { }

    uint64_t doc_cnt() const{
        return m_doc_splitters_rank(m_csa.size());
    }

    uint64_t word_cnt() const {
        return m_csa.size() - doc_cnt();
    }

    size_type serialize(std::ostream& out,
                        structure_tree_node* v = nullptr,
                        std::string name = "") const {
        sdsl::structure_tree_node* child =
            sdsl::structure_tree::add_child(v, name, util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += m_csa.serialize(out, child, "CSA");
        written_bytes += m_doc_splitters.serialize(out, child, "DOCBORDER");
        written_bytes += m_doc_splitters_rank.serialize(out, child,
                "DOCBORDER_RANK");
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

};

namespace internal {
template <typename t_vec>
void build_doc_splitters(sdsl::cache_config& cc, const t_vec& vec) {
    sdsl::bit_vector bv(vec.size());
    for (size_t i = 0; i < vec.size(); ++i)
        if (vec[i] == 1)
            bv[i] = 1;
    sdsl::rrr_vector<> bvr(bv);
    sdsl::rrr_vector<>::rank_1_type rank(&bvr);
    store_to_cache(bvr, surf::KEY_DOCBORDER, cc, true);
    store_to_cache(rank, surf::KEY_DOCBORDER_RANK, cc, true);
}
}

template <typename t_csa>
void construct(idx_brute<t_csa>& idx, const std::string&, sdsl::cache_config& cc,
               uint8_t num_bytes) {
    using t_idx = idx_brute<t_csa>;
    std::cout << "Need CSA" << std::endl;
    if (!cache_file_exists<t_csa>(surf::KEY_CSA, cc)) {
        std::cout << "  building..." << std::endl;
        t_csa csa;
        construct(csa, "", cc, 0);
        store_to_cache(csa, surf::KEY_CSA, cc, true);
    }
    std::cout << "Need doc borders" << std::endl;
    if (!cache_file_exists<sdsl::rrr_vector<>>(surf::KEY_DOCBORDER, cc)) {
        std::cout << "  building..." << std::endl;
        if (std::is_same<typename t_idx::alphabet_category,
                         sdsl::int_alphabet_tag>::value) {
            sdsl::int_vector<> collection;
            load_from_cache(collection, sdsl::conf::KEY_TEXT_INT, cc, false);
            internal::build_doc_splitters(cc, collection);
        } else {
            sdsl::int_vector<8> collection;
            load_from_cache(collection, sdsl::conf::KEY_TEXT, cc, false);
            internal::build_doc_splitters(cc, collection);
        }
    }
    // TODO implement
}

}
