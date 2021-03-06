#pragma once

#include <algorithm>
#include <limits>
#include <memory>
#include <queue>
#include <set>
#include <unordered_set>

#include "sdsl/suffix_trees.hpp"
#include "sdsl/k2_treap.hpp"
#include "surf/construct_col_len.hpp"
#include "surf/df_sada.hpp"
#include "surf/rank_functions.hpp"
#include "surf/topk_interface.hpp"

namespace surf {

using range_type = sdsl::range_type;


template<typename t_select>
struct map_to_dup_type {
    const t_select* m_sel;

    map_to_dup_type(const t_select* select = nullptr) :
        m_sel(select)
    {}

    range_type
    operator()(size_t sp, size_t ep, const bool SINGLETONS = false)const {
        if (SINGLETONS) {
            return {(*m_sel)(sp) + 1, (*m_sel)(ep) - 1};
        }
        uint64_t y = (*m_sel)(ep);
        if (y == 0)
            return {0, (int_vector<>::size_type) - 1};
        uint64_t ep_ = (y + 1) - ep - 1; // # of 0 left to y - 1
        uint64_t sp_ = 0;          // # of 0 left to x
        if (0 == sp) {
        } else {
            uint64_t x = (*m_sel)(sp);
            sp_ = (x + 1) - sp;
        }
        return {sp_, ep_};
    }
};

/*! Class idx_nn consists of a
 *   - CSA over the collection concatenation
 *   - H
 */
template<typename t_csa,
         typename t_k2treap,
         int max_query_length = 0,
         typename t_rmq = sdsl::rmq_succinct_sct<>,
         typename t_border = sdsl::sd_vector<>,
         typename t_border_rank = typename t_border::rank_1_type,
         typename t_border_select = typename t_border::select_1_type,
         typename t_h = sdsl::rrr_vector<63>,
         typename t_h_select_0 = typename t_h::select_0_type,
         typename t_h_select_1 = typename t_h::select_1_type,
         bool     offset_encoding = true,
         typename t_doc_offset = sdsl::hyb_sd_vector<>
         >
class idx_nn
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
    typedef t_rmq                                      rmqc_type;
    typedef t_k2treap                                  k2treap_type;
    typedef k2_treap_ns::top_k_iterator<k2treap_type>  k2treap_iterator;
    typedef typename t_csa::alphabet_category          alphabet_category;
    typedef t_doc_offset                               doc_offset_type;
    typedef typename t_doc_offset::select_1_type       doc_offset_select_type;
    typedef map_to_dup_type<h_select_1_type>           map_to_h_type;
    using topk_interface = typename topk_index_by_alphabet<alphabet_category>::type;

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
    rmqc_type          m_rmqc;
    k2treap_type       m_k2treap;
    map_to_h_type      m_map_to_h;
    topk_result_set    m_results;

    class top_k_iterator : public topk_interface::iter {
    public:
        using k2treap_iterator = k2_treap_ns::top_k_iterator<t_k2treap>;
        typedef std::pair<uint64_t, double> t_doc_val;
        typedef std::stack<std::array<uint64_t, 2>> t_stack_array;
    private:
        const idx_nn*      m_idx;
        uint64_t           m_sp;  // start point of lex interval
        uint64_t           m_ep;  // end point of lex interval
        t_doc_val          m_doc_val;  // stores the current result
        bool               m_valid = false;
        k2treap_iterator   m_k2_iter;
        std::set<uint64_t> m_reported;
        std::set<uint64_t> m_singletons;
        t_stack_array      m_states;
        bool               m_multi_occ = false; // true, if document has to occur more than once
    public:
        top_k_iterator() = delete;
        top_k_iterator(const idx_nn* idx,
                       const typename topk_interface::token_type* begin,
                       const typename topk_interface::token_type* end,
                       bool multi_occ, bool only_match) :
            m_idx(idx), m_multi_occ(multi_occ) {
            m_valid = backward_search(m_idx->m_csa, 0, m_idx->m_csa.size() - 1,
                                      begin, end, m_sp, m_ep) > 0;
            m_valid &= !only_match;
            if (m_valid) {
                auto h_range = m_idx->m_map_to_h(m_sp, m_ep);
                if (!empty(h_range)) {
                    uint64_t depth = end - begin;
                    m_k2_iter = k2_treap_ns::top_k(m_idx->m_k2treap,
                    {std::get<0>(h_range), 0},
                    {std::get<1>(h_range), depth - 1});
                }
                m_states.push({m_sp, m_ep});
                this->next();
            }
        }

        typename topk_interface::snippet_type extract_snippet(const size_t k)
        const override {
            size_type s = (m_doc_val.first == 0)
                          ?  0
                          : (m_idx->m_border_select(m_doc_val.first) + 1);
            size_type e = std::min(s + k,
                                   m_idx->m_border_select(m_doc_val.first + 1) - 1);
            auto res = extract(m_idx->m_csa, s, e);
            return {res.begin(), res.end()};
        }

        void next() override {
            if (m_valid) {
                m_valid = false;
                if (m_k2_iter) {   // multiple occurrence result exists
                    auto xy_w = *m_k2_iter;
                    uint64_t doc_id = offset_encoding
                                      ? m_idx->get_doc(real(xy_w.first))
                                      : m_idx->m_doc[real(xy_w.first)];
                    m_doc_val = t_doc_val(doc_id, xy_w.second + 1);
                    m_reported.insert(doc_id);
                    m_valid = true;
                    ++m_k2_iter;
                } else { // search for singleton results
                    while (!m_multi_occ and !m_states.empty()) {
                        auto state = m_states.top();
                        m_states.pop();
                        uint64_t min_idx = m_idx->m_rmqc(state[0], state[1]);
                        uint64_t doc_id  = m_idx->m_border_rank(m_idx->m_csa[min_idx]);
                        if (m_singletons.find(doc_id) == m_singletons.end()) {
                            m_singletons.insert(doc_id);
                            if (min_idx + 1 <= state[1])
                                m_states.push({min_idx + 1, state[1]});
                            if (state[0] + 1 <= min_idx)
                                m_states.push({state[0], min_idx - 1});
                            if (m_reported.find(doc_id) == m_reported.end()) {
                                m_doc_val = t_doc_val(doc_id, 1);
                                m_reported.insert(doc_id);
                                m_valid = true;
                                break;
                            }
                        }
                    }
                }
            }
        }

        t_doc_val get() const override {
            return m_doc_val;
        }

        bool done() const override {
            return !m_valid;
        }
    };

public:

    std::unique_ptr<typename topk_interface::iter> topk(
        size_t k,
        const typename topk_interface::token_type* begin,
        const typename topk_interface::token_type* end,
        bool multi_occ = false, bool only_match = false) override {
        return std::make_unique<top_k_iterator>(
                   this, begin, end, multi_occ, only_match);
    }

    // Decode m_doc value at postion index by using offset encoding.
    uint64_t get_doc(const uint64_t index) const {
        // All sa offsets are relative to sa_base_pos (rightmost leaf of left subtree).
        uint64_t sa_base_pos = m_h_select_0(index + 1) - index + 1;
        uint64_t base_index = // Index of first dup entry in the node.
            m_h_select_1(sa_base_pos - 1) + 2 - sa_base_pos;
        uint64_t sa_delta = // Extract delta from base index to index.
            (base_index == 0) ?
            m_doc_offset_select(index + 1) :
            m_doc_offset_select(index + 1) - m_doc_offset_select(base_index);
        --sa_delta; // Because zero deltas can't be encoded otherwise.
        uint64_t sa_pos = sa_base_pos + sa_delta;
        uint64_t text_pos = m_csa[sa_pos];
        return m_border_rank(text_pos);
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

    uint64_t doc_cnt() const {
        return m_border_rank(m_csa.size());
    }

    uint64_t word_cnt() const {
        return m_csa.size() - doc_cnt();
    }


    void load(sdsl::cache_config& cc) {
        load_from_cache(m_csa, surf::KEY_CSA, cc, true);
        if (offset_encoding) {
            load_from_cache(m_doc_offset, surf::KEY_DOC_OFFSET, cc, true);
            load_from_cache(m_doc_offset_select, surf::KEY_DOC_OFFSET_SELECT, cc, true);
            m_doc_offset_select.set_vector(&m_doc_offset);
        } else {
            load_from_cache(m_doc, surf::KEY_DUP, cc);
        }
        load_from_cache(m_border, surf::KEY_DOCBORDER, cc, true);
        load_from_cache(m_border_rank, surf::KEY_DOCBORDER_RANK, cc, true);
        m_border_rank.set_vector(&m_border);
        load_from_cache(m_border_select, surf::KEY_DOCBORDER_SELECT, cc, true);
        m_border_select.set_vector(&m_border);
        load_from_cache(m_h, surf::KEY_H, cc, true);
        load_from_cache(m_h_select_0, surf::KEY_H_SELECT_0, cc, true);
        load_from_cache(m_h_select_1, surf::KEY_H_SELECT_1, cc, true);
        m_h_select_0.set_vector(&m_h);
        m_h_select_1.set_vector(&m_h);
        m_map_to_h = map_to_h_type(&m_h_select_1);
        load_from_cache(m_rmqc, surf::KEY_RMQC, cc, true);
        const auto key_w_and_p = (offset_encoding ?
                                  surf::KEY_W_AND_P_G : surf::KEY_W_AND_P) +
                                 std::to_string(max_query_length);
        load_from_cache(m_k2treap, key_w_and_p, cc, true);
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
        written_bytes += m_rmqc.serialize(out, child, "RMQ_C");
        written_bytes += m_k2treap.serialize(out, child, "W_AND_P");
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
        std::cout << sdsl::size_in_bytes(m_rmqc)
                  + sdsl::size_in_bytes(m_k2treap) << std::endl;  // rmq + k2treap
    }
};

template<typename t_cst,
         typename t_select>
struct map_node_to_dup_type {
    typedef typename t_cst::node_type t_node;
    const map_to_dup_type<t_select> m_map;
    const t_cst* m_cst;

    map_node_to_dup_type(const t_select* select, const t_cst* cst):
        m_map(select), m_cst(cst)
    { }

    range_type
    operator()(const t_node& v, const bool SINGLETONS = false)const {
        auto left    = 1;
        auto left_rb = m_cst->rb(m_cst->select_child(v, left));
        return m_map(left_rb, left_rb + 1, SINGLETONS);
    }
    // id \in [1..n-1]
    uint64_t id(const t_node& v)const {
        auto left    = 1;
        return m_cst->rb(m_cst->select_child(v, left)) + 1;
    }
};

template<typename t_csa,
         typename t_k2treap,
         int max_query_length,
         typename t_rmq,
         typename t_border,
         typename t_border_rank,
         typename t_border_select,
         typename t_h,
         typename t_h_select_0,
         typename t_h_select_1,
         bool     offset_encoding,
         typename t_doc_offset
         >
void construct(idx_nn<t_csa, t_k2treap, max_query_length, t_rmq, t_border, t_border_rank,
               t_border_select, t_h, t_h_select_0, t_h_select_1, offset_encoding,
               t_doc_offset>& idx, const std::string&, sdsl::cache_config& cc,
               uint8_t num_bytes) {
    using namespace sdsl;
    using namespace std;
    using t_df = DF_TYPE;
    using cst_type = typename t_df::cst_type;
    using t_wtd = WTD_TYPE;
    using idx_type = idx_nn<t_csa, t_k2treap, max_query_length, t_rmq, t_border, t_border_rank,
          t_border_select, t_h, t_h_select_0, t_h_select_1, offset_encoding, t_doc_offset>;
    using doc_offset_type = typename idx_type::doc_offset_type;

    construct_col_len<t_df::alphabet_category::WIDTH>(cc);

    const auto key_w_and_p = (offset_encoding ?
                              surf::KEY_W_AND_P_G : surf::KEY_W_AND_P) +
                             std::to_string(max_query_length);
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
        load_from_cache(h, surf::KEY_H, cc);
        t_h hrrr = h;
        store_to_cache(hrrr, surf::KEY_H, cc, true);
        t_h_select_0 h_select_0(&hrrr);
        t_h_select_1 h_select_1(&hrrr);
        store_to_cache(h_select_0, surf::KEY_H_SELECT_0, cc, true);
        store_to_cache(h_select_1, surf::KEY_H_SELECT_1, cc, true);
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
            cout << dup << endl;
        }

        std::string P_file = cache_file_name(key_p, cc);
        int_vector_buffer<> P_buf(P_file, std::ios::out, 1 << 20,
                                  sdsl::bits::hi(max_depth) + 1);

        t_wtd wtd;
        load_from_cache(wtd, surf::KEY_WTD, cc, true);

        t_h hrrr;
        load_from_cache(hrrr, surf::KEY_H, cc, true);
        t_h_select_1 h_select_1;
        load_from_cache(h_select_1, surf::KEY_H_SELECT_1, cc, true);
        h_select_1.set_vector(&hrrr);
        cst_type cst;
        load_from_file(cst, cache_file_name<cst_type>(surf::KEY_TMPCST, cc));
        map_node_to_dup_type<cst_type, t_h_select_1> map_node_to_dup(&h_select_1, &cst);

        uint64_t doc_cnt = 1;
        load_from_cache(doc_cnt, KEY_DOCCNT, cc);
        typedef stack<uint32_t, vector<uint32_t>> t_stack;
        //  HELPER to build the pointer structure
        vector<t_stack> depths(doc_cnt, t_stack(vector<uint32_t>(1, 0))); // doc_cnt stack for last depth
        uint64_t depth = 0;

        // DFS traversal of CST
        for (auto it = cst.begin(); it != cst.end(); ++it) {
            auto v = *it; // get the node by dereferencing the iterator
            if (!cst.is_leaf(v)) {
                if (it.visit() == 1) {
                    // node visited the first time
                    depth = cst.depth(v);
                    range_type r = map_node_to_dup(v);
                    if (!empty(r)) {
                        for (size_t i = get<0>(r); i <= get<1>(r); ++i) {
                            depths[dup[i]].push(depth);
                        }
                    }
                } else { // node visited the second time
                    range_type r = map_node_to_dup(v);
                    if (!empty(r)) {
                        for (size_t i = get<0>(r); i <= get<1>(r); ++i) {
                            depths[dup[i]].pop();
                            P_buf[i] = depths[dup[i]].top();
                        }
                    }
                }
            }
        }
        P_buf.close();
    }
    if (offset_encoding) {
        cout << "...DOC_OFFSET" << endl;
        if (!cache_file_exists<doc_offset_type>(surf::KEY_DOC_OFFSET, cc)) {
            int_vector<> darray, dup;
            load_from_cache(darray, surf::KEY_DARRAY, cc);
            load_from_cache(dup, surf::KEY_DUP_G, cc);
            t_h hrrr;
            load_from_cache(hrrr, surf::KEY_H, cc, true);
            t_h_select_1 h_select_1;
            load_from_cache(h_select_1, surf::KEY_H_SELECT_1, cc, true);
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

    cout << "...RMQ_C" << endl;
    if (!cache_file_exists<t_rmq>(surf::KEY_RMQC, cc))
    {
        int_vector<> C;
        load_from_cache(C, surf::KEY_C, cc);
        t_rmq rmq_c(&C);
        store_to_cache(rmq_c, surf::KEY_RMQC, cc, true);
    }
    cout << "...W_AND_P" << endl;
    if (!cache_file_exists<t_k2treap>(key_w_and_p, cc))
    {
        int_vector_buffer<> P_buf(cache_file_name(key_p, cc));
        std::string W_and_P_file = cache_file_name(key_w_and_p, cc);
        cout << "P_buf.size()=" << P_buf.size() << endl;
        // Build filter bitvector.
        bit_vector add_to_grid_bv(P_buf.size(), 1);
        if (max_query_length > 0) {
            int_vector<> P;
            load_from_cache(P, key_p, cc);
            uint64_t removed_count = 0;
            for (size_t i = 0; i < P.size(); ++i)
                if (P[i] > max_query_length) {
                    add_to_grid_bv[i] = 0;
                    removed_count++;
                }
            cout << "Removed " << removed_count << " from " << P.size() << " grid points\n";
        }
        {
            int_vector<> id_v(P_buf.size(), 0, bits::hi(P_buf.size()) + 1);
            util::set_to_id(id_v);
            filter(id_v, add_to_grid_bv);
            if (id_v.size() < 30) cout << "x=" << id_v << endl;
            store_to_file(id_v, W_and_P_file + ".x");
        }
        {
            int_vector<> P;
            load_from_cache(P, key_p, cc);
            filter(P, add_to_grid_bv);
            if (P.size() < 30) cout << "y=" << P << endl;
            store_to_file(P, W_and_P_file + ".y");
        }
        {
            int_vector<> W;
            load_from_cache(W, key_weights, cc);
            filter(W, add_to_grid_bv);
            if (W.size() < 30) cout << "w=" << W << endl;
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
