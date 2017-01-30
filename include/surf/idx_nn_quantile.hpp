#pragma once

#include <algorithm>
#include <limits>
#include <memory>
#include <queue>
#include <set>
#include <unordered_set>

#include "btree/safe_btree_set.h"
#include "sdsl/rrr_vector.hpp"
#include "sdsl/suffix_trees.hpp"
#include "sdsl/k2_treap.hpp"
#include "surf/construct_col_len.hpp"
#include "surf/df_sada.hpp"
#include "surf/rank_functions.hpp"
#include "surf/topk_interface.hpp"

namespace surf {

/*! Class idx_nn_quantile consists of a
 *   - CSA over the collection concatenation
 *   - H
 */
template<typename t_csa,
         typename t_k2treap,
         int      quantile = 1, // max query time slow down by 32.
         int max_query_length = 0,
         typename t_border = sdsl::sd_vector<>,
         typename t_border_rank = typename t_border::rank_1_type,
         typename t_border_select = typename t_border::select_1_type,
         typename t_h = sdsl::rrr_vector<63>,
         typename t_h_select_0 = typename t_h::select_0_type,
         typename t_h_select_1 = typename t_h::select_1_type,
         bool     offset_encoding = false,
         typename t_doc_offset = sdsl::hyb_sd_vector<>
         >
class idx_nn_quantile
    : public topk_index_by_alphabet<typename t_csa::alphabet_category>::type {
public:
    using size_type = sdsl::int_vector<>::size_type;
    typedef t_csa                                      csa_type;
    typedef t_border                                   border_type;
    typedef t_border_rank                              border_rank_type;
    typedef t_border_select                            border_select_type;
    typedef t_h                                        h_type;
    typedef t_h_select_0                               h_select_0_type;
    typedef t_h_select_1                               h_select_1_type;
    typedef t_k2treap                                  k2treap_type;
    typedef k2_treap_ns::top_k_iterator<k2treap_type>  k2treap_iterator;
    typedef typename t_csa::alphabet_category          alphabet_category;
    typedef t_doc_offset                               doc_offset_type;
    typedef typename t_doc_offset::select_1_type       doc_offset_select_type;
    typedef map_to_dup_type<h_select_1_type>           map_to_h_type;

    using qfilter_type = rrr_vector<>;

    using topk_interface = typename topk_index_by_alphabet<alphabet_category>::type;

    static constexpr string QUANTILE_SUFFIX() {
        return string("_q") + std::to_string(quantile);
    };

private:
    csa_type           m_csa;
    border_type        m_border;
    border_rank_type   m_border_rank;
    border_select_type m_border_select;
    h_type             m_h;
    h_select_0_type    m_h_select_0;
    h_select_1_type    m_h_select_1;
    typename h_type::rank_1_type m_h_rank;
    doc_offset_type    m_doc_offset; // offset representation of documents in node list
    doc_offset_select_type m_doc_offset_select;
    int_vector<>       m_doc; // documents in node lists
    k2treap_type       m_k2treap;
    map_to_h_type      m_map_to_h;
    topk_result_set    m_results;

    qfilter_type m_quantile_filter;
    qfilter_type::rank_1_type m_quantile_filter_rank;
    qfilter_type::select_1_type m_quantile_filter_select;

    // Using k2treap.
    void getTopK(k2treap_iterator k2_iter, uint64_t k) {
        while (k2_iter && k != 0) {
            auto xy_w = *k2_iter;
            uint64_t arrow_id = real(xy_w.first);
            uint64_t doc_id;
            if (offset_encoding) {
                uint64_t h_id = m_quantile_filter_select(arrow_id+1);
                if (m_h[h_id]) { // Singleton.
                    doc_id = sa_to_doc(m_h_rank(h_id));
                } else{
                    uint64_t ones = m_h_rank(h_id); // Offset id.
                    uint64_t zeros = h_id - ones;
                    uint64_t zeros_start = m_h_select_1(ones) - ones;
                    assert(zeros_start < zeros);
                    uint64_t p = ones + 
                        m_doc_offset_select(zeros+2) - m_doc_offset_select(zeros_start+2);
                    doc_id = sa_to_doc(p-1);
                }
            } else {
                //std::cerr << "arrow_id = " << arrow_id << " #m_doc=" << m_doc.size() << endl;
                doc_id = m_doc[arrow_id];
                //std::cerr << "decode: " << doc_id << " " << arrow_id << std::endl;
            }
            m_results.push_back(topk_result(doc_id, xy_w.second));
            ++k2_iter;
            k--;
        }
    }
    // Naive fallback.
    void getTopK(uint64_t s, uint64_t e) {
        std::unordered_map<uint64_t, uint64_t> counts;
        //std::cerr << s << "---" << e << std::endl;
        for (size_t i = s; i <= e; ++i) {
            uint64_t doc_id = sa_to_doc(i);
            counts[doc_id]++;
        }
        // TODO only take top k.
        for (const auto res : counts)
            m_results.push_back(topk_result(res.first, res.second));
    }


public:

    std::unique_ptr<typename topk_interface::iter> topk(
        size_t k,
        const typename topk_interface::token_type* begin,
        const typename topk_interface::token_type* end,
        bool multi_occ, bool only_match) override {
            using std::get;
            m_results.clear();
            uint64_t sp, ep;
            bool valid = backward_search(m_csa, 0, m_csa.size() - 1,
                                      begin, end, sp, ep) > 0;
            //std::cerr << "sp ep = " << sp << " " << ep << std::endl;
            valid &= !only_match;
            uint64_t interval_size = 0;

            if (valid) {
                interval_size = ep - sp + 1;
                //std::cerr<< "interval size: " << interval_size << " " << quantile << " " << k << std::endl;
                // interval_size > 1 handles the special case interval_size = quantile = k = 1
                if (interval_size >= k*quantile && interval_size > 1) { // Use grid.
                    //std::cerr << "using grid" << std::endl;
                    uint64_t depth = end - begin;

                    // round up to succeeding sample
                    uint64_t from = m_quantile_filter_rank(m_h_select_1(sp+1));
                    // round down to preceding sample (`to` is exclusive!)
                    uint64_t to = m_quantile_filter_rank(m_h_select_1(ep+1)+1);

                    //std::cerr
                        //<< "x range = " << get<0>(h_range) << " " << get<1>(h_range)
                        //<< " q range = " << from << "-" << to
                        //<< " depth range = 0-" << depth-1 << std::endl;

                    if (from < to) {
                        --to;
                        if (from <= to) {
                            getTopK(k2_treap_ns::top_k(m_k2treap,
                                        {from, 0},
                                        {to, depth - 1}),
                                    k);
                        }
                    }
                } else { // Naive fallback.
                    //std::cerr << "fallback" << std::endl;
                    getTopK(sp, ep);
                }
            }

            if (this->get_debug_stream())
                (*this->get_debug_stream()) << "INTERVAL_SIZE;" << interval_size << "\n";
            return sort_topk_results<typename topk_interface::token_type>(&m_results);
    }

    // Decode m_doc value at postion index by using offset encoding.
    uint64_t get_doc(const uint64_t index_zero, const uint64_t index_one) const {
        // All sa offsets are relative to sa_base_pos (rightmost leaf of left subtree).
        //uint64_t sa_base_pos = m_h_select_0(index + 1) - index + 1;
        uint64_t sa_base_pos = index_one + 1;
        uint64_t base_index = // Index of first dup entry in the node.
            m_h_select_1(sa_base_pos - 1) + 2 - sa_base_pos;
        uint64_t sa_delta = // Extract delta from base index to index.
            (base_index == 0) ?
            m_doc_offset_select(index_zero + 1) :
            m_doc_offset_select(index_zero + 1) - m_doc_offset_select(base_index);
        --sa_delta; // Because zero deltas can't be encoded otherwise.
        return sa_to_doc(sa_base_pos + sa_delta);
    }

    auto doc(uint64_t doc_id) -> decltype(extract(m_csa, 0, 0)) {
        size_type doc_begin = 0;
        if (doc_id) {
            doc_begin = m_border_select(doc_id) + 1;
        }
        size_type doc_end = m_border_select(doc_id + 1) - 1;
        auto res = extract(m_csa, doc_begin, doc_end);
        return res;
    }

    uint64_t sa_to_doc(const uint64_t sa_pos) const {
        return m_border_rank(m_csa[sa_pos]);
    }

    uint64_t doc_cnt() const {
        return m_border_rank(m_csa.size());
    }

    uint64_t word_cnt() const {
        return m_csa.size() - doc_cnt();
    }


    void load(sdsl::cache_config& cc) {
        load_from_cache(m_csa, surf::KEY_CSA, cc, true);
        if (offset_encoding) {
            // TODO(niklasb) don't forget to use QUANTILE_SUFFIX if we ever fix this
            load_from_cache(m_doc_offset, surf::KEY_DOC_OFFSET + QUANTILE_SUFFIX(), cc, true);
            load_from_cache(m_doc_offset_select, surf::KEY_DOC_OFFSET_SELECT + QUANTILE_SUFFIX(), cc, true);
            m_doc_offset_select.set_vector(&m_doc_offset);
        } else {
            load_from_cache(m_doc, surf::KEY_DUP_G + QUANTILE_SUFFIX(), cc);
        }

        load_from_cache(m_quantile_filter,
                surf::KEY_FILTERED_QUANTILE_FILTER + QUANTILE_SUFFIX(), cc, true);
        load_from_cache(m_quantile_filter_rank,
                surf::KEY_FILTERED_QUANTILE_FILTER_RANK + QUANTILE_SUFFIX(), cc, true);
        load_from_cache(m_quantile_filter_select,
                surf::KEY_FILTERED_QUANTILE_FILTER_SELECT + QUANTILE_SUFFIX(), cc, true);
        m_quantile_filter_rank.set_vector(&m_quantile_filter);
        m_quantile_filter_select.set_vector(&m_quantile_filter);

        load_from_cache(m_border, surf::KEY_DOCBORDER, cc, true);
        load_from_cache(m_border_rank, surf::KEY_DOCBORDER_RANK, cc, true);
        m_border_rank.set_vector(&m_border);
        load_from_cache(m_border_select, surf::KEY_DOCBORDER_SELECT, cc, true);
        m_border_select.set_vector(&m_border);

        load_from_cache(m_h, surf::KEY_FILTERED_H + QUANTILE_SUFFIX(), cc, true);
        load_from_cache(m_h_select_0, surf::KEY_FILTERED_H_SELECT_0 + QUANTILE_SUFFIX(),
                cc, true);
        load_from_cache(m_h_select_1, surf::KEY_FILTERED_H_SELECT_1 + QUANTILE_SUFFIX(),
                cc, true);
        load_from_cache(m_h_rank, surf::KEY_FILTERED_H_RANK + QUANTILE_SUFFIX(),
                cc, true);
        m_h_select_0.set_vector(&m_h);
        m_h_select_1.set_vector(&m_h);
        m_h_rank.set_vector(&m_h);

        m_map_to_h = map_to_h_type(&m_h_select_1);

        auto key_w_and_p = surf::KEY_W_AND_P_G;
        load_from_cache(m_k2treap, key_w_and_p + QUANTILE_SUFFIX(), cc, true);
    }

    size_type serialize(std::ostream& out, structure_tree_node* v = nullptr,
                        std::string name = "")const {
        structure_tree_node* child = structure_tree::add_child(v, name,
                                     util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += m_csa.serialize(out, child, "CSA");
        if (offset_encoding) {
            written_bytes += m_doc_offset.serialize(out, child, "DOC_OFFSET");
            written_bytes += m_doc_offset_select.serialize(out, child,
                             "DOC_OFFSET_SELECT");
        } else {
            written_bytes += m_doc.serialize(out, child, "DOC");
        }
        written_bytes += m_border.serialize(out, child, "BORDER");
        written_bytes += m_border_rank.serialize(out, child, "BORDER_RANK");
        written_bytes += m_h.serialize(out, child, "H");
        written_bytes += m_h_select_0.serialize(out, child, "H_SELECT_0");
        written_bytes += m_h_select_1.serialize(out, child, "H_SELECT_1");
        written_bytes += m_h_rank.serialize(out, child, "H_RANK");
        written_bytes += m_k2treap.serialize(out, child, "W_AND_P");
        written_bytes += m_quantile_filter.serialize(out, child, "QUANTILE_FILTER");
        written_bytes += m_quantile_filter_rank.serialize(out, child, "QUANTILE_FILTER_RANK");
        written_bytes += m_quantile_filter_select.serialize(out, child, "QUANTILE_FILTER_SELECT");
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void mem_info()const {
        std::cout << "Dupsize " << m_doc.size() << std::endl;
        std::cout << sdsl::size_in_bytes(m_csa) +
                  sdsl::size_in_bytes(m_border) +
                  sdsl::size_in_bytes(m_border_rank) << ";"; // CSA
        if (offset_encoding) {
            std::cout << sdsl::size_in_bytes(m_doc_offset)
                      + sdsl::size_in_bytes(m_doc_offset_select) << ";"; // DOC
        } else {
            std::cout << sdsl::size_in_bytes(m_doc) << ";"; // DOC
        }
        std::cout << sdsl::size_in_bytes(m_h)
                  + sdsl::size_in_bytes(m_h_select_0)
                  + sdsl::size_in_bytes(m_h_select_1)
                  + sdsl::size_in_bytes(m_h_rank)
                  + sdsl::size_in_bytes(m_quantile_filter)
                  + sdsl::size_in_bytes(m_quantile_filter_rank)
                  + sdsl::size_in_bytes(m_quantile_filter_select)<< ";"; // H
        std::cout << sdsl::size_in_bytes(m_k2treap) << ";";  // k2treap
        std::cout << m_quantile_filter_rank(m_quantile_filter.size()) << std::endl; // |G_q|.
    }
};

using arrow = std::pair<uint64_t, uint64_t>;
struct arrow_cmp {
    bool operator()(const arrow& a, const arrow& b) const {
        return a.first > b.first || (a.first==b.first && a.second > b.second);
    }
};

template<typename t_csa,
         typename t_k2treap,
         int quantile,
         int max_query_length,
         typename t_border,
         typename t_border_rank,
         typename t_border_select,
         typename t_h,
         typename t_h_select_0,
         typename t_h_select_1,
         bool     offset_encoding,
         typename t_doc_offset
         >
void construct(idx_nn_quantile<t_csa, t_k2treap, quantile, max_query_length, t_border, t_border_rank,
               t_border_select, t_h, t_h_select_0, t_h_select_1, offset_encoding,
               t_doc_offset>& idx, const std::string&, sdsl::cache_config& cc,
               uint8_t num_bytes) {
    using namespace sdsl;
    using namespace std;
    using t_df = DF_TYPE;
    using cst_type = typename t_df::cst_type;
    using t_wtd = WTD_TYPE;
    using idx_type = idx_nn_quantile<t_csa, t_k2treap, quantile, max_query_length, t_border, t_border_rank,
          t_border_select, t_h, t_h_select_0, t_h_select_1, offset_encoding, t_doc_offset>;
    using doc_offset_type = typename idx_type::doc_offset_type;
    using timer = chrono::high_resolution_clock;
    using qfilter_type = rrr_vector<>;

    construct_col_len<t_df::alphabet_category::WIDTH>(cc);

    auto key_w_and_p = surf::KEY_W_AND_P_G + idx_type::QUANTILE_SUFFIX();
    const auto key_p = surf::KEY_P_QUANTILE_G + idx_type::QUANTILE_SUFFIX();
    const auto key_dup = surf::KEY_DUP_G;
    const auto key_weights = surf::KEY_WEIGHTS_G;
    const auto key_df = surf::KEY_SADADF_G;

    cout << "...CSA" << endl; // CSA to get the lex. range
    if (!cache_file_exists<t_csa>(surf::KEY_CSA, cc))
    {
        t_csa csa;
        construct(csa, "", cc, 0);
        store_to_cache(csa, surf::KEY_CSA, cc, true);
    }
    cout << "...WTD" << endl; // Document array and wavelet tree of it
    if (!cache_file_exists<t_wtd>(surf::KEY_WTD, cc)) {
        construct_darray<t_csa::alphabet_type::int_width>(cc, false);
        t_wtd wtd;
        construct(wtd, cache_file_name(surf::KEY_DARRAY, cc), cc);
        cout << "wtd.size() = " << wtd.size() << endl;
        cout << "wtd.sigma = " << wtd.sigma << endl;
        store_to_cache(wtd, surf::KEY_WTD, cc, true);
    }
    cout << "...DF" << endl; //
    if (!cache_file_exists<t_df>(key_df, cc))
    {
        t_df df;
        construct(df, "", cc, 0);
        store_to_cache(df, key_df, cc, true);
        bit_vector h;
        load_from_cache(h, surf::KEY_H_LEFT, cc);
        t_h hrrr = h;
        store_to_cache(hrrr, surf::KEY_H_LEFT, cc, true);
        t_h_select_0 h_select_0(&hrrr);
        t_h_select_1 h_select_1(&hrrr);
        store_to_cache(h_select_0, surf::KEY_H_LEFT_SELECT_0, cc, true);
        store_to_cache(h_select_1, surf::KEY_H_LEFT_SELECT_1, cc, true);
    }
    cout << "...DOC_BORDER" << endl;
    if (!cache_file_exists<t_border>(surf::KEY_DOCBORDER, cc) or
            !cache_file_exists<t_border_rank>(surf::KEY_DOCBORDER_RANK, cc) or
            !cache_file_exists<t_border_select>(surf::KEY_DOCBORDER_SELECT, cc))
    {
        bit_vector doc_border;
        load_from_cache(doc_border, surf::KEY_DOCBORDER, cc);
        t_border sd_doc_border(doc_border);
        store_to_cache(sd_doc_border, surf::KEY_DOCBORDER, cc, true);
        t_border_rank doc_border_rank(&sd_doc_border);
        store_to_cache(doc_border_rank, surf::KEY_DOCBORDER_RANK, cc, true);
        t_border_select doc_border_select(&sd_doc_border);
        store_to_cache(doc_border_select, surf::KEY_DOCBORDER_SELECT, cc, true);
    }
    cout << "...WTD" << endl;
    if (!cache_file_exists<t_wtd>(surf::KEY_WTD, cc)) {
        construct_darray<t_csa::alphabet_type::int_width>(cc, false);
        t_wtd wtd;
        construct(wtd, cache_file_name(surf::KEY_DARRAY, cc), cc);
        cout << "wtd.size() = " << wtd.size() << endl;
        cout << "wtd.sigma = " << wtd.sigma << endl;
        store_to_cache(wtd, surf::KEY_WTD, cc, true);
    }
// P corresponds to up-pointers
    cout << "...P" << endl;
    if (!cache_file_exists(key_p, cc))
    {
        uint64_t max_depth = 0;
        load_from_cache(max_depth, surf::KEY_MAXCSTDEPTH, cc);

        int_vector<> dup;
        load_from_cache(dup, key_dup, cc);
        cout << "dup.size()=" << dup.size() << endl;
        if (dup.size() < 20) {
            cout << "dup=" << dup << endl;
        }

        std::string P_file = cache_file_name(key_p, cc);
        int_vector_buffer<> P_buf(P_file, std::ios::out, 1 << 20,
                                  sdsl::bits::hi(max_depth) + 1);

        t_wtd wtd;
        load_from_cache(wtd, surf::KEY_WTD, cc, true);

        t_h hrrr;
        load_from_cache(hrrr, surf::KEY_H_LEFT, cc, true);
        t_h_select_1 h_select_1;
        t_h_select_0 h_select_0;
        load_from_cache(h_select_1, surf::KEY_H_LEFT_SELECT_1, cc, true);
        load_from_cache(h_select_0, surf::KEY_H_LEFT_SELECT_0, cc, true);
        h_select_1.set_vector(&hrrr);
        h_select_0.set_vector(&hrrr);
        cst_type cst;
        load_from_file(cst, cache_file_name<cst_type>(surf::KEY_TMPCST, cc));
        map_node_to_dup_type<cst_type, t_h_select_1> map_node_to_dup(&h_select_1, &cst);

        uint64_t doc_cnt = 1;
        load_from_cache(doc_cnt, KEY_DOCCNT, cc);
        typedef stack<uint32_t, vector<uint32_t>> t_stack;
        //  HELPER to build the pointer structure
        vector<t_stack> depths(doc_cnt, t_stack(vector<uint32_t>(1, 0))); // doc_cnt stack for last depth
        uint64_t depth = 0;

        auto start = timer::now();

        // empty string, first index in H mapping
        P_buf[0] = 0;

        // DFS traversal of CST
        for (auto it = cst.begin(); it != cst.end(); ++it) {
            auto v = *it; // get the node by dereferencing the iterator
            //cout << "node " << v.i << "-" << v.j << " visit=" << (int)it.visit() << endl;
            if (!cst.is_leaf(v)) {
                auto left_rb = cst.rb(cst.select_child(v, 1));
                //cout << "  left_rb=" << left_rb  << endl;
                auto x1 = h_select_1(left_rb + 1);
                auto left1 = x1 - left_rb;
                auto x2 = h_select_1(left_rb + 2);
                auto left2 = x2 - left_rb - 1;
                //cout << "  " << x1 << " " << x2 << " " << left1 << " "<< left2 << endl;
                range_type r = { left1, left2 - 1 };
                //range_type r = map_node_to_dup(v);
                //cout << "  r=" << get<0>(r) << "-" << get<1>(r) << endl;

                if (it.visit() == 1) {
                    // node visited the first time
                    depth = cst.depth(v);

                    if (!empty(r)) {
                        //cout << "  non-empty" << endl;
                        for (size_t i = get<0>(r); i <= get<1>(r); ++i) {
                            depths[dup[i]].push(depth);
                        }
                    }
                } else { // node visited the second time
                    //range_type r = map_node_to_dup(v);
                    if (!empty(r)) {
                        for (size_t i = get<0>(r); i <= get<1>(r); ++i) {
                            depths[dup[i]].pop();
                            uint64_t idx = h_select_0(i+1);
                            //cout << "  " << i << " " << idx << " "
                                //<< depths[dup[i]].top() << endl;
                            P_buf[idx] = depths[dup[i]].top();
                        }
                    }
                }
            } else if (v.i > 0) {
                uint64_t sa_pos = v.i;
                uint64_t d = wtd[sa_pos];
                uint64_t idx = h_select_1(sa_pos+1);
                //cout << "  singleton " << sa_pos << " " << idx << " " << depths[d].top() << endl;
                if (d < depths.size())
                    P_buf[idx] = depths[d].top();
            }
        }

        uint64_t msecs =
            chrono::duration_cast<chrono::microseconds>(timer::now() - start).count();
        cout << "Computing P took " << setprecision(2) << fixed
            << 1.*msecs/1e6 << " seconds" << endl;

        P_buf.close();
    }

    cout << "...quantile filter" << endl;
    if (!cache_file_exists<bit_vector>(
                surf::KEY_QUANTILE_FILTER + idx_type::QUANTILE_SUFFIX(), cc)) {
        // Compute quantile filter. 1 -> include arrow, 0 -> remove arrow.
        t_h hrrr;
        load_from_cache(hrrr, KEY_H_LEFT, cc, true);
        if (hrrr.size() < 30)
            cout << "H = " << hrrr << endl;
        t_h_select_1 h_select_1;
        load_from_cache(h_select_1, KEY_H_LEFT_SELECT_1, cc, true);
        h_select_1.set_vector(&hrrr);

        // TODO(niklasb) why is hrrr one too large?
        const uint64_t bits =  hrrr.size() - 1;
        bit_vector quantile_filter(bits, 0);
        int_vector<> weights;
        load_from_file(weights, cache_file_name(key_weights, cc));
        if (weights.size() < 20)
            cout << "weights = " << weights << endl;

        cst_type cst;
        load_from_file(cst, cache_file_name<cst_type>(KEY_TMPCST, cc));
        map_node_to_dup_type<cst_type, t_h_select_1> map_node_to_dup(&h_select_1, &cst);
        std::cout << hrrr.size() << " " << weights.size() << std::endl;

        int_vector<> P;
        load_from_cache(P, key_p, cc);
        std::cout << "P.size()=" << P.size() << endl;
        if (P.size() <30) cout << "P=" << P << endl;

        auto start = timer::now();

        using arrow_set = btree::safe_btree_set<arrow, arrow_cmp>;
        map<uint64_t, arrow_set*> arrows;

        vector<cst_type::node_type> path;

        int depth = 0;

        // DFS traversal of CST
        for (auto it = cst.begin(); it != cst.end(); ++it) {
            auto v = *it; // get the node by dereferencing the iterator

            bool leaf = cst.is_leaf(v);

            if (!leaf)
                depth += (it.visit() == 1) ? 1 : -1;

            //cout << "node " << v.i << "-" << v.j
                    //<< " visit=" << (int)it.visit() << " depth=" << depth << endl;
            if (!leaf && depth > 0 && it.visit() == 2) {
                // node visited the second time
                auto depth = cst.depth(v);

                // Locate children in the arrows set (ordered by left node border)
                auto first_child = arrows.find(v.i);
                assert(first_child != arrows.end());

                // Find child with maximum arrows. We are going to use it as
                // the set of arrows for the current node.
                arrow_set* cur = nullptr;
                for (auto it = first_child; it != arrows.end(); ++it) {
                    //cout << "  child=" << it->second.first.i <<"-" << it->second.first.j << endl;
                    auto* a = it->second;
                    if (!cur || (a && a->size() > cur->size()))
                        cur = a;
                }
                assert(cur); // because we're not at a leaf
                // Merge all other children.
                for (auto it = first_child; it != arrows.end(); ++it) {
                    arrow_set* a = it->second;
                    if (a == cur) continue;

                    for (auto arrow : (*a))
                        if (P[arrow.second] < depth)
                            cur->insert(arrow);
                    delete a;
                }
                arrows.erase(next(first_child), arrows.end());
                first_child->second = cur;

                // Insert repetitions associated with current node
                auto left_rb = cst.rb(cst.select_child(v, 1));
                auto x = h_select_1(left_rb+1);
                auto weight_idx = x - left_rb;
                ++x;
                while (x < bits && !hrrr[x]) {
                    auto weight = weights[weight_idx];
                    cur->insert(arrow(weight, x));
                    ++weight_idx;
                    ++x;
                }

                uint64_t interval_size = v.j - v.i + 1;
                uint64_t k = interval_size / quantile;

                // TODO(niklasb) instead of explicitly walking the tree to mark the
                // quantiles, can we instead use order statistics and check if
                // an arrow is in some top quantile as soon as we see it?
                // The problem then would be that we couldn't do the deletions
                // lazily, like we do now.
                for (auto it = cur->begin(); k && it != cur->end();) {
                    if (P[it->second] < depth) {
                        quantile_filter[it->second] = 1;
                        ++it;
                        --k;
                    } else {
                        // Erase arrows fully contained in current subtree
                        cur->erase(it++);
                    }
                }

                // We can delete arrows[v.i] here, if v is
                // a child of the root node.
                if (depth == 1) {
                    delete cur;
                    arrows.erase(v.i);
                }
            } else if (leaf && it.visit() == 1) {
                auto x = h_select_1(v.i+1);
                auto* cur = new arrow_set();
                cur->insert(arrow(0, x));
                arrows[v.i] = cur;
            }
        }
        for (auto it : arrows) delete it.second;
        arrows.clear();

        uint64_t msecs = chrono::duration_cast<chrono::microseconds>(timer::now() - start).count();
        cout << "quantile filtering took " << setprecision(2) << fixed
            << 1.*msecs/1e6 << " seconds" << endl;

        size_t cnt_needed = 0;
        for (uint64_t i = 0; i < quantile_filter.size(); ++i) {
            if (quantile_filter[i]) {
                cnt_needed++;
            }
        }
        cout << "total arrows: " << bits << " after filter: " << cnt_needed << endl;

        store_to_cache(quantile_filter,
                surf::KEY_QUANTILE_FILTER + idx_type::QUANTILE_SUFFIX(), cc, true);
    }

    cout << "...filtering H vector" << endl;
    if (!cache_file_exists<qfilter_type>(surf::KEY_FILTERED_QUANTILE_FILTER + idx_type::QUANTILE_SUFFIX(), cc)) {
        t_h hrrr;
        load_from_cache(hrrr, KEY_H_LEFT, cc, true);
        auto bits = hrrr.size() - 1;

        bit_vector qfilter;
        load_from_cache(qfilter,
                KEY_QUANTILE_FILTER  + idx_type::QUANTILE_SUFFIX(),
                cc, true);

        bit_vector filtered_h(bits, 0);
        bit_vector filtered_qfilter(bits, 0);
        size_t j = 0;
        for (size_t i = 0; i < bits; ++i) {
            if (hrrr[i]) {
                filtered_h[j] = 1;
                filtered_qfilter[j] = qfilter[i];
                ++j;
            } else if (qfilter[i]) {
                filtered_qfilter[j] = 1;
                ++j;
            }
        }
        filtered_h.resize(j);
        filtered_qfilter.resize(j);

        qfilter_type qfilter_filtered(filtered_qfilter);
        qfilter_type::rank_1_type qfilter_filtered_rank(&qfilter_filtered);
        qfilter_type::select_1_type qfilter_filtered_select(&qfilter_filtered);

        t_h h_filtered(filtered_h);
        t_h_select_1 h_filtered_select_1(&h_filtered);
        t_h_select_0 h_filtered_select_0(&h_filtered);
        typename t_h::rank_1_type h_filtered_rank(&h_filtered);

        store_to_cache(h_filtered,
                KEY_FILTERED_H + idx_type::QUANTILE_SUFFIX(),
                cc, true);
        store_to_cache(h_filtered_select_1,
                KEY_FILTERED_H_SELECT_1 + idx_type::QUANTILE_SUFFIX(),
                cc, true);
        store_to_cache(h_filtered_select_0,
                KEY_FILTERED_H_SELECT_0 + idx_type::QUANTILE_SUFFIX(),
                cc, true);
        store_to_cache(h_filtered_rank,
                KEY_FILTERED_H_RANK + idx_type::QUANTILE_SUFFIX(),
                cc, true);

        store_to_cache(qfilter_filtered_rank,
                KEY_FILTERED_QUANTILE_FILTER_RANK + idx_type::QUANTILE_SUFFIX(),
                cc, true);
        store_to_cache(qfilter_filtered_select,
                KEY_FILTERED_QUANTILE_FILTER_SELECT + idx_type::QUANTILE_SUFFIX(),
                cc, true);
        store_to_cache(qfilter_filtered,
                KEY_FILTERED_QUANTILE_FILTER + idx_type::QUANTILE_SUFFIX(),
                cc, true);
    }
   cout << "...W_AND_P" << endl;
    if (!cache_file_exists<t_k2treap>(key_w_and_p, cc) ||
            (offset_encoding && !cache_file_exists<doc_offset_type>(surf::KEY_DOC_OFFSET + idx_type::QUANTILE_SUFFIX(), cc)) ||
            (!offset_encoding && !cache_file_exists(key_dup + idx_type::QUANTILE_SUFFIX(), cc)))
    {
        int_vector_buffer<> P_buf(cache_file_name(key_p, cc));
        std::string W_and_P_file = cache_file_name(key_w_and_p, cc);
        t_h hrrr;
        load_from_cache(hrrr, surf::KEY_H_LEFT, cc, true);
        cout << "P_buf.size()=" << P_buf.size() << endl;

        bit_vector quantile_filter;
        load_from_cache(quantile_filter,
                surf::KEY_QUANTILE_FILTER + idx_type::QUANTILE_SUFFIX(), cc, true);

        {
            size_t keep = 0;
            for (size_t i = 0; i < quantile_filter.size(); ++i)
                keep += quantile_filter[i];
            int_vector<> id_v(keep, 0, bits::hi(keep) + 1);
            util::set_to_id(id_v);
            cout << "x size = " << id_v.size() << endl;
            store_to_file(id_v, W_and_P_file + ".x");
        }
        {
            int_vector<> P;
            load_from_cache(P, key_p, cc);
            filter(P, quantile_filter);
            cout << "y size = " << P.size() << endl;
            if (P.size() < 30) cout << "y=" << P << endl;
            //cout << "y(191)=" << P[191] << endl;
            store_to_file(P, W_and_P_file + ".y");
        }
        {
            int_vector<> darray;
            load_from_cache(darray, surf::KEY_DARRAY, cc);

            int_vector<> dup_nosingletons;
            load_from_cache(dup_nosingletons, key_dup, cc);

            int_vector<> dup(dup_nosingletons);
            dup.resize(quantile_filter.size());

            int_vector<> W_nosingletons;
            load_from_cache(W_nosingletons, key_weights, cc);

            int_vector<> W(W_nosingletons);
            W.resize(quantile_filter.size());

            // Add singletons and filter by quantiles.
            uint64_t dup_idx = 0;
            uint64_t idx = 0;
            uint64_t sa_id = 0;
            uint64_t last_sa_id = -1;
            uint64_t last_val = 0;
            vector<uint64_t> sa_offset;
            for (size_t i = 0; i < quantile_filter.size(); ++i) {
                if (quantile_filter[i]) {
                    if (offset_encoding) { 
                        if (!hrrr[i]) {
                            assert(i == dup_idx + sa_id);
                            uint64_t d = dup_nosingletons[dup_idx];
                            uint64_t j = last_sa_id == sa_id ? last_val : sa_id;
                            while (darray[j] != d) {
                                j++;
                            }
                            if (last_sa_id == sa_id) {
                                //assert(sa_offset[sa_offset.size() -1] < j - sa_id);
                                assert(j > last_val);
                                sa_offset.push_back(j - last_val);
                            } else {
                                // + 1 because 0 cannot be encoded.
                                sa_offset.push_back(j - sa_id + 1);
                            }
                            last_val = j;
                            last_sa_id = sa_id;
                        }
                    } else {
                        dup[idx] = hrrr[i] ? darray[sa_id] : dup_nosingletons[dup_idx];
                    }
                    W[idx] = hrrr[i] ? 1 : (W_nosingletons[dup_idx]+1);
                    ++idx;
                }
                if(hrrr[i])  {
                    sa_id++;
                } else dup_idx++;
            }
            W.resize(idx);
            if (offset_encoding) {
                cout << "offset encoding..." << endl;
                // Compute number of bits.
                uint64_t sd_n = 1; // add 1 at the beginning.
                for (const auto delta : sa_offset) {
                   sd_n += delta; 
                }
                bit_vector plain_bv(sd_n+1, 0);
                plain_bv[0] = 1;
                uint64_t pos = 1;
                for (const auto delta: sa_offset) {
                    pos += delta;
                    plain_bv[pos] = 1;
                }
                doc_offset_type doc_offset(plain_bv);
                typename doc_offset_type::select_1_type doc_offset_select(&doc_offset);
                store_to_cache(doc_offset, surf::KEY_DOC_OFFSET + idx_type::QUANTILE_SUFFIX(), cc, true);
                store_to_cache(doc_offset_select, surf::KEY_DOC_OFFSET_SELECT + idx_type::QUANTILE_SUFFIX(), cc, true);
            } else {
                dup.resize(idx);
                store_to_cache(dup, key_dup + idx_type::QUANTILE_SUFFIX(), cc);
            }
            cout << "w size = " << W.size() << endl;
            if (W.size() < 30) cout << "w=" << W << endl;
            //cout << "w(191)=" << W[191] << endl;
            store_to_file(W, W_and_P_file + ".w");
        }
        cout << "build k2treap" << endl;
        t_k2treap k2treap;
        construct(k2treap, cache_file_name(key_w_and_p, cc));
        store_to_cache(k2treap, key_w_and_p, cc, true);
        sdsl::remove(W_and_P_file + ".x");
        sdsl::remove(W_and_P_file + ".y");
        sdsl::remove(W_and_P_file + ".w");
    }
}

} // end namespace surf
