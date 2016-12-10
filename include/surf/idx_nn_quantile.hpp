#pragma once

#include <algorithm>
#include <limits>
#include <memory>
#include <queue>
#include <set>
#include <unordered_set>

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
    doc_offset_type    m_doc_offset; // offset representation of documents in node list
    doc_offset_select_type m_doc_offset_select;
    int_vector<>       m_doc; // documents in node lists
    k2treap_type       m_k2treap;
    map_to_h_type      m_map_to_h;
    topk_result_set    m_results;

    rrr_vector<>       m_quantile_filter;
    rrr_vector<>::rank_1_type m_quantile_filter_rank;

    // Using k2treap.
    void getTopK(k2treap_iterator k2_iter, uint64_t k) {
        while (k2_iter && k != 0) {
            auto xy_w = *k2_iter;
            uint64_t arrow_id = real(xy_w.first);
            uint64_t doc_id;
            if (offset_encoding) {
                std::cerr << "Not implemented." << std::endl;
                abort();
                // TODO implement and test.
                //uint64_t ones = m_h_rank(arrow_id);
                //uint64_t zeros = arrow_id - ones;
                //if (m_h[arrow_id]) { // singelton.
                //   doc_id = sa_to_doc(ones);
                //} else {
                //    doc_id = get_doc(zeros, ones);
                //}
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
            if (counts.count(doc_id) == 0)
                counts[doc_id] = 0;
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
        bool multi_occ = false, bool only_match = false) override {
            using std::get;
            m_results.clear();
            uint64_t sp,ep;
            bool valid = backward_search(m_csa, 0, m_csa.size() - 1,
                                      begin, end, sp, ep) > 0;
            //std::cerr << "sp ep = " << sp << " " << ep << std::endl;
            valid &= !only_match;
            if (valid) {
                //TODO(niklasb) we should check which one is correct
                //and use that in map_to_dup_type
                range_type h_range {
                    m_h_select_1(sp+1),
                    m_h_select_1(ep+1) };

                //std::cerr << get<0>(h_range) << " " << get<1>(h_range) << std::endl;
                if (!empty(h_range)) {
                    uint64_t interval_size = get<1>(h_range) + 1 - get<0>(h_range);
                    //std::cerr<< "interval size: " << interval_size << " " << quantile << " " << k << std::endl;

                    // interval_size > 1 handles the special case interval_size = quantile = k = 1
                    if (interval_size >= k*quantile && interval_size > 1) { // Use grid.
                        //std::cerr << "using grid" << std::endl;
                        uint64_t depth = end - begin;

                        // round up to succeeding sample
                        uint64_t from = m_quantile_filter_rank(get<0>(h_range));
                        // round down to preceding sample (border is exclusive!)
                        uint64_t to = m_quantile_filter_rank(get<1>(h_range)+1);

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
            }
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
            load_from_cache(m_doc_offset, surf::KEY_DOC_OFFSET, cc, true);
            load_from_cache(m_doc_offset_select, surf::KEY_DOC_OFFSET_SELECT, cc, true);
            m_doc_offset_select.set_vector(&m_doc_offset);
        } else {
            load_from_cache(m_doc, surf::KEY_DUP + QUANTILE_SUFFIX(), cc);
        }

        load_from_cache(m_quantile_filter,
                surf::KEY_QUANTILE_FILTER + QUANTILE_SUFFIX(), cc);
        load_from_cache(m_quantile_filter_rank,
                surf::KEY_QUANTILE_FILTER_RANK + QUANTILE_SUFFIX(), cc);
        m_quantile_filter_rank.set_vector(&m_quantile_filter);

        load_from_cache(m_border, surf::KEY_DOCBORDER, cc, true);
        load_from_cache(m_border_rank, surf::KEY_DOCBORDER_RANK, cc, true);
        m_border_rank.set_vector(&m_border);
        load_from_cache(m_border_select, surf::KEY_DOCBORDER_SELECT, cc, true);
        m_border_select.set_vector(&m_border);
        load_from_cache(m_h, surf::KEY_H_LEFT, cc, true);
        load_from_cache(m_h_select_0, surf::KEY_H_LEFT_SELECT_0, cc, true);
        load_from_cache(m_h_select_1, surf::KEY_H_LEFT_SELECT_1, cc, true);
        m_h_select_0.set_vector(&m_h);
        m_h_select_1.set_vector(&m_h);
        m_map_to_h = map_to_h_type(&m_h_select_1);

        auto key_w_and_p =
            (offset_encoding ? surf::KEY_W_AND_P_G : surf::KEY_W_AND_P) +
            std::to_string(max_query_length);
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
        written_bytes += m_k2treap.serialize(out, child, "W_AND_P");
        written_bytes += m_quantile_filter.serialize(out, child, "QUANTILE_FILTER");
        written_bytes += m_quantile_filter_rank.serialize(out, child, "QUANTILE_FILTER_RANK");
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
                  + sdsl::size_in_bytes(m_h_select_1) << ";"; // H
        std::cout << sdsl::size_in_bytes(m_k2treap) << std::endl;  // k2treap
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

    assert(!offset_encoding);

    construct_col_len<t_df::alphabet_category::WIDTH>(cc);

    auto key_w_and_p = (offset_encoding ?
                              surf::KEY_W_AND_P_G : surf::KEY_W_AND_P) +
                             std::to_string(max_query_length);
    key_w_and_p += idx_type::QUANTILE_SUFFIX();

    const auto key_p = offset_encoding ?
                       surf::KEY_P_G : surf::KEY_P;
    const auto key_dup = offset_encoding ?
                         surf::KEY_DUP_G : surf::KEY_DUP;
    const auto key_weights = offset_encoding ?
                             surf::KEY_WEIGHTS_G : surf::KEY_WEIGHTS;
    const auto key_df = offset_encoding ?
                        surf::KEY_SADADF_G : surf::KEY_SADADF;

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
                P_buf[idx] = depths[d].top();
            }
        }
        P_buf.close();
    }

    cout << "...quantile filter" << endl;
    if (!cache_file_exists(surf::KEY_QUANTILE_FILTER + idx_type::QUANTILE_SUFFIX(), cc)) {
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
        if (P.size() <30) cout << "P=" << P << endl;

        std::vector<std::tuple<uint64_t, uint64_t>> cur_weights;

        // DFS traversal of CST
        for (auto it = cst.begin(); it != cst.end(); ++it) {
            auto v = *it; // get the node by dereferencing the iterator
            if (!cst.is_leaf(v)) {
                if (it.visit() == 1) {
                    // node visited the first time
                    uint64_t start = h_select_1(v.i+1);
                    uint64_t end = h_select_1(v.j+1)+1;

                    // TODO(niklasb) is this still correct?
                    uint64_t weight_idx = start - v.i;

                    uint64_t interval_size = end - start;
                    uint64_t k = interval_size / quantile;

                    if (k > 0) {
                        auto depth = cst.depth(v);
                        cur_weights.clear();
                        for (uint64_t i = start; i < end; ++i) {
                            auto valid = P[i] < depth;
                            if (hrrr[i] == 1) { // singleton.
                                if (valid) cur_weights.emplace_back(0, i);
                            } else { // no singleton.
                                if (valid) cur_weights.emplace_back(weights[weight_idx], i);
                                ++weight_idx;
                            }
                        }

                        k = min(k, cur_weights.size());
                        std::nth_element(cur_weights.begin(), cur_weights.begin()+k, cur_weights.end(),
                                std::greater<std::tuple<uint64_t, uint64_t>>());
                        for (size_t i = 0; i < k; ++i)
                            quantile_filter[get<1>(cur_weights[i])] = 1;
                    }
                }
            }
        }
        size_t cnt_needed = 0;
        for (uint64_t i = 0; i < quantile_filter.size(); ++i) {
            if (quantile_filter[i]) {
                cnt_needed++;
            }
        }
        cout << "total arrows: " << bits << " after filter: " << cnt_needed << endl;

        rrr_vector<> qfilter(quantile_filter);
        rrr_vector<>::rank_1_type qfilter_rank(&qfilter);

        store_to_cache(qfilter,
                surf::KEY_QUANTILE_FILTER + idx_type::QUANTILE_SUFFIX(), cc);
        store_to_cache(qfilter_rank,
                surf::KEY_QUANTILE_FILTER_RANK + idx_type::QUANTILE_SUFFIX(), cc);
    }
    if (offset_encoding) {
        cout << "...DOC_OFFSET" << endl;
        if (!cache_file_exists<doc_offset_type>(surf::KEY_DOC_OFFSET, cc)) {
            int_vector<> darray, dup;
            load_from_cache(darray, surf::KEY_DARRAY, cc);
            load_from_cache(dup, key_dup, cc);
            t_h hrrr;
            load_from_cache(hrrr, surf::KEY_H_LEFT, cc, true);
            t_h_select_1 h_select_1;
            load_from_cache(h_select_1, surf::KEY_H_LEFT_SELECT_1, cc, true);
            h_select_1.set_vector(&hrrr);

            // Iterate through all nodes.
            uint64_t start = 0;
            uint64_t end;
            // For each dup value offset o such that darray[nodeIndex+o+1] = dup.
            vector<uint64_t> sa_offset;
            uint64_t sd_n = 1;
            for (uint64_t i = 1; i <= darray.size(); ++i) {
                end = h_select_1(i) + 1 - i;
                if (start < end) { // Dup lens
                    vector<uint64_t> dup_set(dup.begin() + start, dup.begin() + end);
                    uint64_t sa_pos = i;
                    uint64_t j = 0;
                    while (j < dup_set.size()) {
                        if (dup_set[j] == darray[sa_pos]) {
                            // Store offset.
                            sa_offset.push_back(sa_pos - i);
                            ++j;
                        }
                        if (sa_pos >= darray.size()) {
                            cout << "ERROR: sa_pos is out of bounds." << endl;
                            return;
                        }
                        ++sa_pos;
                    }
                    // sd_n computation.
                    // encode first value + 1.
                    sd_n += sa_offset[start] + 1; // +1 because zero deltas can't be encoded.
                    for (size_t j = start + 1; j < end; ++j)
                        sd_n += sa_offset[j] - sa_offset[j - 1]; // encode deltas.
                }
                start = end;
            }
            sdsl::bit_vector plain_bv(sd_n);
            start = 0;
            uint64_t cur_pos = 0;
            for (uint64_t i = 1; i < darray.size(); ++i) {
                end = h_select_1(i) + 1 - i;
                if (start < end) { // Dup lens
                    cur_pos += sa_offset[start] + 1;
                    plain_bv[cur_pos] = 1;
                    for (size_t j = start + 1; j < end; ++j) {
                        cur_pos += sa_offset[j] - sa_offset[j - 1];
                        plain_bv[cur_pos] = 1;
                    }
                }
                start = end;
            }
            doc_offset_type doc_offset(plain_bv);
            // Build select.
            typename doc_offset_type::select_1_type doc_offset_select(&doc_offset);
            store_to_cache(doc_offset, surf::KEY_DOC_OFFSET, cc, true);
            store_to_cache(doc_offset_select, surf::KEY_DOC_OFFSET_SELECT, cc, true);
        }
    }
    cout << "...W_AND_P" << endl;
    if (!cache_file_exists<t_k2treap>(key_w_and_p, cc))
    {
        int_vector_buffer<> P_buf(cache_file_name(key_p, cc));
        std::string W_and_P_file = cache_file_name(key_w_and_p, cc);
        t_h hrrr;
        load_from_cache(hrrr, surf::KEY_H_LEFT, cc, true);
        cout << "P_buf.size()=" << P_buf.size() << endl;

        rrr_vector<> quantile_filter;
        rrr_vector<>::rank_1_type qfilter_rank(&quantile_filter);
        load_from_cache(quantile_filter,
                surf::KEY_QUANTILE_FILTER + idx_type::QUANTILE_SUFFIX(), cc);
        load_from_cache(qfilter_rank,
                surf::KEY_QUANTILE_FILTER_RANK + idx_type::QUANTILE_SUFFIX(), cc);
        qfilter_rank.set_vector(&quantile_filter);

        {
            auto keep = qfilter_rank(quantile_filter.size());
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
            for (size_t i = 0; i < quantile_filter.size(); ++i) {
                if (quantile_filter[i]) {
                    dup[idx] = hrrr[i] ? darray[sa_id] : dup_nosingletons[dup_idx];
                    W[idx] = hrrr[i] ? 1 : (W_nosingletons[dup_idx]+1);
                    ++idx;
                }
                if(hrrr[i]) sa_id++;
                else dup_idx++;
            }
            dup.resize(idx);
            W.resize(idx);
            store_to_cache(dup, key_dup + idx_type::QUANTILE_SUFFIX(), cc);
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
