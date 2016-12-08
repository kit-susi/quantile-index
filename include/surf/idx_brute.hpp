#pragma once

#include <algorithm>
#include <limits>
#include <queue>
#include <set>

#include "sdsl/int_vector.hpp"
#include "sdsl/sdsl_concepts.hpp"
#include "sdsl/suffix_trees.hpp"
#include "surf/util.hpp"
#include "surf/topk_interface.hpp"

namespace surf {

/*! Trivial implementation of top-k retrieval using brute force search, for
 * testing.
 */
template <typename t_csa>
class idx_brute
    : public topk_index_by_alphabet<typename t_csa::alphabet_category>::type {
public:
    using alphabet_category = typename t_csa::alphabet_category;
    using topk_interface = typename topk_index_by_alphabet<alphabet_category>::type;
    using csa_type = t_csa;
    using size_type = sdsl::int_vector<>::size_type;

private:
    using token_type = typename topk_interface::token_type;

    csa_type m_csa;
    sdsl::rrr_vector<> m_doc_splitters;
    sdsl::rrr_vector<>::rank_1_type m_doc_splitters_rank;
    topk_result_set m_results;

public:
    void load(sdsl::cache_config& cc) {
        load_from_cache(m_csa, surf::KEY_CSA, cc, true);
        load_from_cache(m_doc_splitters, surf::KEY_DOCBORDER, cc, true);
        load_from_cache(m_doc_splitters_rank, surf::KEY_DOCBORDER_RANK, cc, true);
        m_doc_splitters_rank.set_vector(&m_doc_splitters);
    }

    std::unique_ptr<typename topk_interface::iter> topk(
        size_t k, const token_type* begin, const token_type* end,
        bool multi_occ = false, bool only_match = false) override {
        auto occs = locate(m_csa, begin, end);

        std::map<uint64_t, double> occs_by_doc;
        for (auto pos : occs) {
            //std::cerr << "occurrence " << pos << std::endl;
            auto doc = m_doc_splitters_rank(pos);
            occs_by_doc[doc] += 1;
        }

        m_results.clear();
        for (auto it : occs_by_doc)
            if (!multi_occ || it.second > 1)
                m_results.emplace_back(it.first, it.second);
        return sort_topk_results<token_type>(&m_results);
    }

    std::unique_ptr<typename topk_interface::iter> topk_intersect(
        size_t k, const typename topk_interface::intersect_query& query,
        bool multi_occ = false, bool only_match = false) override {
        std::map<uint64_t, double> by_doc;

        bool first = true;
        for (const auto& q : query) {
            auto occs = locate(m_csa, q.first, q.second);
            std::map<uint64_t, double> by_doc_new;
            for (auto pos : occs) {
                auto doc = m_doc_splitters_rank(pos);
                if (first || by_doc.count(doc))
                    by_doc_new[doc] += 1;
            }
            std::vector<uint64_t> erase_docs;
            for (auto& it : by_doc_new) {
                if (multi_occ && it.second <= 1)
                    erase_docs.push_back(it.first);
                else
                    it.second += by_doc[it.first];
            }
            for (uint64_t doc : erase_docs)
                by_doc_new.erase(doc);

            std::swap(by_doc_new, by_doc);
            first = false;
        }

        m_results = topk_result_set(by_doc.begin(), by_doc.end());
        return sort_topk_results<token_type>(&m_results);
    }

    void mem_info() const { }

    uint64_t doc_cnt() const {
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
void construct(idx_brute<t_csa>& idx,
               const std::string&, sdsl::cache_config& cc,
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

}  // namespace surf
