#pragma once
#include <stack>
#include <vector>
const bool SINGLETONS = true;

#include "sdsl/suffix_trees.hpp"
#include "surf/topk_interface.hpp"
#include "sdsl/rmq_succinct_sct.hpp"

namespace surf {

//template<typename t_cst,
//         typename t_select>
//struct map_node_to_dup_type {
//    typedef typename t_cst::node_type t_node;
//    const map_to_dup_type<t_select> m_map;
//    const t_cst* m_cst;
//
//    map_node_to_dup_type(const t_select* select, const t_cst* cst):
//        m_map(select), m_cst(cst)
//    { }
//
//    range_type
//    operator()(const t_node& v)const {
//        auto left    = 1;
//        auto left_rb = m_cst->rb(m_cst->select_child(v, left));
//        return m_map(left_rb, left_rb + 1);
//    }
//    // id \in [1..n-1]
//    uint64_t id(const t_node& v)const {
//        auto left    = 1;
//        return m_cst->rb(m_cst->select_child(v, left)) + 1;
//    }
//};


using rmq_type = rmq_succinct_sct<false>;
template<typename t_csa,
         typename t_border = sdsl::sd_vector<>,
         typename t_border_rank = typename t_border::rank_1_type,
         typename t_border_select = typename t_border::select_1_type,
         int LEVELS = 10,
         typename t_tails = sdsl::wt_huff<>>
class idx_top_down
    : public topk_index_by_alphabet<typename t_csa::alphabet_category>::type {

public:
    typedef t_csa                                      csa_type;
    typedef t_border                                   border_type;
    typedef t_border_rank                              border_rank_type;
    typedef t_border_select                            border_select_type;
    typedef t_tails                                    tails_type;
    using size_type = sdsl::int_vector<>::size_type;
    typedef typename t_csa::alphabet_category          alphabet_category;
    using topk_interface = typename topk_index_by_alphabet<alphabet_category>::type;

    typedef sdsl::rrr_vector<63>                       h_type;
    typedef h_type::rank_1_type                        h_rank_type;
    typedef h_type::select_1_type                      h_select_1_type;
    typedef map_to_dup_type<h_select_1_type>           map_to_h_type;
    typedef sdsl::hyb_sd_vector<>                      doc_offset_type;
    typedef typename doc_offset_type::select_1_type    doc_offset_select_type;

private:
    using token_type = typename topk_interface::token_type;

    csa_type           m_csa;
    border_type        m_border;
    border_rank_type   m_border_rank;
    border_select_type m_border_select;
    tails_type         m_tails;
    dac_vector<>       m_weights;       // Only nonsingletons are stored in H order!
    rmq_type           m_weights_rmq;   // Tails order including singletons.
    doc_offset_type    m_doc_offset;
    doc_offset_select_type m_doc_offset_select;
    h_type             m_h;
    h_rank_type        m_h_rank;
    h_select_1_type    m_h_select_1;
    map_to_h_type      m_map_to_h;

public:

    template <typename t_token>
    class top_down_topk_iterator : public topk_iterator<t_token> {
    private:
        const idx_top_down* m_idx;
        uint64_t           m_sp;  // start point of lex interval
        uint64_t           m_ep;  // end point of lex interval
        bool               m_valid = false;
        bool               m_multi_occ = false; // true, if document has to occur more than once
        struct top_down_result {
            double weight;
            uint64_t document;
            uint64_t tails_id;
            // Interval in [in tails ids].
            uint64_t start;
            uint64_t end; // exclusive.
            // Only for navigation purposes.
            uint64_t offset; // Start of level in tails id.
            uint64_t depth; // Tree depth of interval.
            bool operator<(const top_down_result& obj) const {
                return weight < obj.weight;
            }
        };
        std::priority_queue<top_down_result> m_intervals;
        std::unordered_set<uint64_t> m_reported_docs;

        void init_interval(top_down_result& interval) {
            interval.tails_id = m_idx->m_weights_rmq(interval.start, interval.end - 1);
            assert(interval.tails_id >= interval.start && interval.tails_id < interval.end);
            uint64_t h_pos = m_idx->m_tails.select(
                                 interval.tails_id - interval.offset + 1, interval.depth);
            // Zeros in H are duplicates, and ones are singletons.
            uint64_t singletons = m_idx->m_h_rank(h_pos);
            uint64_t nonsingletons = h_pos - singletons;
            if (m_idx->m_h[h_pos] == 1) { // is singleton.
                interval.weight = 1;
                interval.document = m_idx->sa_to_doc(singletons + 1);
            } else {
                interval.weight = m_idx->m_weights[nonsingletons] + 1;
                interval.document = m_idx->get_doc(nonsingletons, singletons);
            }
        }

    public:
        top_down_topk_iterator() = delete;
        top_down_topk_iterator(const idx_top_down* idx,
                               const typename topk_interface::token_type* begin,
                               const typename topk_interface::token_type* end,
                               bool multi_occ, bool only_match) :
            m_idx(idx), m_multi_occ(multi_occ) {
            m_valid = backward_search(m_idx->m_csa, 0, m_idx->m_csa.size() - 1,
                                      begin, end, m_sp, m_ep) > 0;
            m_valid &= !only_match;
            if (m_valid) {
                auto h_range = m_idx->m_map_to_h(m_sp, m_ep, SINGLETONS);
                std::get<0>(h_range)--;
                std::get<1>(h_range)++;

                if (!empty(h_range)) {
                    uint64_t offset = 0;
                    // <= here instead of < ??? TODO
                    for (uint64_t depth = 0; depth <= (end - begin); ++depth) {
                        top_down_result interval;
                        interval.start = idx->m_tails.rank(std::get<0>(h_range), depth) + offset;
                        interval.end = idx->m_tails.rank(std::get<1>(h_range) + 1, depth) + offset;
                        interval.offset = offset;
                        interval.depth = depth;
                        if (interval.start < interval.end) {
                            init_interval(interval);
                            m_intervals.push(interval);
                        }
                        offset += idx->m_tails.rank(idx->m_tails.size(), depth); // TODO dont do this for last.
                    }
                }
            }
        }

        topk_result get() const {
            if (m_intervals.empty())
                abort();
            return topk_result(m_intervals.top().document,
                               m_intervals.top().weight);
        }

        bool done() const override {
            return m_intervals.empty();
        }

        void next() override {
            if (m_intervals.empty())
                abort();
            auto t = m_intervals.top();
            m_reported_docs.insert(t.document);
            m_intervals.pop();
            top_down_result interval;
            // Left side.
            interval.start = t.start;
            interval.end = t.tails_id;
            interval.offset = t.offset;
            interval.depth = t.depth;
            if (interval.start < interval.end) {
                init_interval(interval);
                m_intervals.push(interval);
            }
            // Right side.
            interval.start = t.tails_id + 1;
            interval.end = t.end;
            interval.offset = t.offset;
            interval.depth = t.depth;
            if (interval.start < interval.end) {
                init_interval(interval);
                m_intervals.push(interval);
            }
            // Remove already reported docs from front.
            while (!m_intervals.empty() &&
                    m_reported_docs.count(m_intervals.top().document) == 1) {
                m_intervals.pop();
            }
        }

        std::vector<t_token> extract_snippet(const size_t k) const override {
            // TODO implement
            return {};
        }
    private:
    };

    std::unique_ptr<typename topk_interface::iter> topk(
        size_t k,
        const token_type* begin,
        const token_type* end,
        bool multi_occ = false, bool only_match = false) override {
        return std::make_unique<top_down_topk_iterator<token_type>>(
                   top_down_topk_iterator<token_type>(this,
                           begin, end, multi_occ, only_match));
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
        load_from_cache(m_border, surf::KEY_DOCBORDER, cc, true);
        load_from_cache(m_border_rank, surf::KEY_DOCBORDER_RANK, cc, true);
        m_border_rank.set_vector(&m_border);
        load_from_cache(m_border_select, surf::KEY_DOCBORDER_SELECT, cc, true);
        m_border_select.set_vector(&m_border);
        load_from_cache(m_tails, surf::KEY_TAILS + std::to_string(LEVELS), cc, true);
        load_from_cache(m_weights, surf::KEY_WEIGHTS_G + std::to_string(LEVELS), cc, true);
        load_from_cache(m_weights_rmq, surf::KEY_WEIGHTS_RMQ + std::to_string(LEVELS), cc, true);
        load_from_cache(m_doc_offset, surf::KEY_DOC_OFFSET + std::to_string(LEVELS), cc, true);
        load_from_cache(m_doc_offset_select, surf::KEY_DOC_OFFSET_SELECT + std::to_string(LEVELS), cc, true);
        m_doc_offset_select.set_vector(&m_doc_offset);

        load_from_cache(m_h, surf::KEY_H, cc, true);
        load_from_cache(m_h_rank, surf::KEY_H_RANK, cc, true);
        load_from_cache(m_h_select_1, surf::KEY_H_SELECT_1, cc, true);
        m_h_rank.set_vector(&m_h);
        m_h_select_1.set_vector(&m_h);
        m_map_to_h = map_to_h_type(&m_h_select_1);
    }

    size_type serialize(std::ostream& out, structure_tree_node* v = nullptr,
                        std::string name = "")const {
        structure_tree_node* child = structure_tree::add_child(v, name,
                                     util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += m_csa.serialize(out, child, "CSA");
        written_bytes += m_border.serialize(out, child, "BORDER");
        written_bytes += m_border_rank.serialize(out, child, "BORDER_RANK");
        written_bytes += m_border_select.serialize(out, child, "BORDER_SELECT");
        written_bytes += m_tails.serialize(out, child, "TAILS");
        written_bytes += m_weights.serialize(out, child, "WEIGHTS");
        written_bytes += m_weights_rmq.serialize(out, child, "WEIGHTSRMQ");
        written_bytes += m_doc_offset.serialize(out, child, "DOC_OFFSET");
        written_bytes += m_doc_offset_select.serialize(out, child, "DOC_OFFSET_SELECT");
        written_bytes += m_h.serialize(out, child, "H");
        written_bytes += m_h_rank.serialize(out, child, "H_RANK");
        written_bytes += m_h_select_1.serialize(out, child, "H_SELECT_1");
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void mem_info()const {
        std::cout << sdsl::size_in_bytes(m_csa) +
                  sdsl::size_in_bytes(m_border) +
                  sdsl::size_in_bytes(m_border_rank) +
                  sdsl::size_in_bytes(m_border_select) << ";"; // CSA

        std::cout << sdsl::size_in_bytes(m_doc_offset)
                  + sdsl::size_in_bytes(m_doc_offset_select) << ";"; // DOC
        std::cout << sdsl::size_in_bytes(m_tails) << ";"; // TAILS
        std::cout << sdsl::size_in_bytes(m_weights) +
                  sdsl::size_in_bytes(m_weights_rmq) << ";"; // WEIGHTS
    }
};

template<typename t_csa,
         typename t_border,
         typename t_border_rank,
         typename t_border_select,
         int LEVELS,
         typename t_tails>
void construct(idx_top_down<t_csa,
               t_border, t_border_rank,
               t_border_select,
               LEVELS,
               t_tails>& idx,
               const std::string&, sdsl::cache_config& cc, uint8_t num_bytes) {
    using t_wtd = WTD_TYPE;
    using t_df = surf::df_sada<sdsl::rrr_vector<63>, sdsl::rrr_vector<63>::select_1_type, sdsl::byte_alphabet_tag>;
    using cst_type = t_df::cst_type;
    using node_type = typename cst_type::node_type;
    using tails_type = t_tails;
    using t_h = sdsl::rrr_vector<63>;
    using t_h_rank = t_h::rank_1_type;
    using t_h_select_1 = t_h::select_1_type;

    const auto key_dup = surf::KEY_DUP_G;
    const auto key_df = surf::KEY_SADADF_G;

    cout << "...CSA" << endl; // CSA to get the lex. range
    if (!cache_file_exists<t_csa>(surf::KEY_CSA, cc))
    {
        t_csa csa;
        construct(csa, "", cc, 0);
        store_to_cache(csa, surf::KEY_CSA, cc, true);
    }
    cout << "...WTD" << endl;
    // Document array and wavelet tree of it
    // Note: This also constructs doc borders.
    if (!cache_file_exists<t_wtd>(surf::KEY_WTD, cc)) {
        construct_darray<t_csa::alphabet_type::int_width>(cc, false);
        t_wtd wtd;
        construct(wtd, cache_file_name(surf::KEY_DARRAY, cc), cc);
        cout << "wtd.size() = " << wtd.size() << endl;
        cout << "wtd.sigma = " << wtd.sigma << endl;
        store_to_cache(wtd, surf::KEY_WTD, cc, true);
    }
    cout << "...DF" << endl; // For h vector and repetition array.
    if (!cache_file_exists<t_df>(key_df, cc))
    {
        t_df df;
        construct(df, "", cc, 0);
        store_to_cache(df, key_df, cc, true);
        bit_vector h;
        load_from_cache(h, surf::KEY_H, cc);
        t_h hrrr(h);
        store_to_cache(hrrr, surf::KEY_H, cc, true);
        t_h_select_1 h_select_1(&hrrr);
        store_to_cache(h_select_1, surf::KEY_H_SELECT_1, cc, true);
        t_h_rank h_rank(&hrrr);
        store_to_cache(h_rank, surf::KEY_H_RANK, cc, true);
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

    cout << "...tails" << endl;
    const auto key_tails = surf::KEY_TAILS + std::to_string(LEVELS);
    if (!cache_file_exists<tails_type>(key_tails, cc))
    {
        uint64_t max_depth = 0;
        load_from_cache(max_depth, surf::KEY_MAXCSTDEPTH, cc);

        int_vector<> dup;
        load_from_cache(dup, key_dup, cc);
        cout << "dup.size()=" << dup.size() << endl;
        if (dup.size() < 20) {
            cout << dup << endl;
        }

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

        int_vector<> old_weights;
        load_from_file(old_weights, cache_file_name(surf::KEY_WEIGHTS_G, cc));
        const uint64_t bits_per_level = SINGLETONS ? hrrr.size() : old_weights.size();
        int_vector<8>  tails_plain(bits_per_level);

        uint64_t doc_cnt = 1;
        load_from_cache(doc_cnt, KEY_DOCCNT, cc);
        typedef std::stack<uint32_t, std::vector<uint32_t>> t_stack;
        //  HELPER to build the pointer structure
        std::vector<t_stack> depths(doc_cnt + 1, t_stack(std::vector<uint32_t>(1, 0))); // doc_cnt stack for last depth
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
                        for (size_t i = std::get<0>(r); i <= std::get<1>(r); ++i) {
                            depths[dup[i]].push(depth);
                        }
                    }
                } else { // node visited the second time
                    range_type r = map_node_to_dup(v);
                    if (!empty(r)) {
                        range_type sr = map_node_to_dup(v, SINGLETONS);
                        assert(std::get<0>(sr) - std::get<0>(r) == std::get<1>(sr) - std::get<1>(r));
                        // Get index of first arrow leaving node.
                        for (size_t i = std::get<0>(r); i <= std::get<1>(r); ++i) {
                            depths[dup[i]].pop();
                            uint64_t idx = i + std::get<0>(sr) - std::get<0>(r);
                            uint64_t depth = depths[dup[i]].top();
                            assert(!tails_plain[idx]);
                            tails_plain[idx] = std::min(depth, (uint64_t)255);
                        }
                    }
                }
            } else if (SINGLETONS && v.i > 0) { // First node is empty string.
                uint64_t sa_pos = v.i;
                uint64_t d = wtd[sa_pos];
                uint64_t depth = depths[d].top();
                uint64_t idx = h_select_1(sa_pos);
                assert(!tails_plain[idx]);
                tails_plain[idx] = depth;
            }
        }
        tails_type tails;
        sdsl::construct_im(tails, tails_plain);
        store_to_cache(tails, key_tails, cc, true);
    }
    cout << "...weights and rmq.";
    const auto key_weights = surf::KEY_WEIGHTS_G + std::to_string(LEVELS);
    const auto key_weights_rmq = surf::KEY_WEIGHTS_RMQ + std::to_string(LEVELS);
    if (!cache_file_exists<rmq_type>(key_weights_rmq, cc)) {
        // Reorder weights.
        cout << "Construct rmq." << endl;
        int_vector<> old_weights;
        int_vector<> dup;
        load_from_file(old_weights, cache_file_name(surf::KEY_WEIGHTS_G, cc));
        load_from_file(dup, cache_file_name(surf::KEY_DUP_G, cc));

        int_vector<> D;
        load_from_file(D, cache_file_name(surf::KEY_DARRAY, cc));
        t_h hrrr;
        load_from_cache(hrrr, surf::KEY_H, cc, true);
        t_h::rank_1_type h_rank(&hrrr);
        tails_type tails;
        load_from_file(tails, cache_file_name(key_tails, cc));
        const uint64_t bits_per_level = SINGLETONS ? hrrr.size() : old_weights.size();
        int_vector<> weights(bits_per_level, 0, old_weights.width());
        cout << bits_per_level << " vs. " << old_weights.size() << endl;
        {
            uint64_t offsets[256];
            offsets[0] = 0;
            // Init offsets.
            for (uint64_t i = 1; i < 256; ++i) {
                offsets[i] = offsets[i - 1] + tails.rank(tails.size(), i - 1);
            }
            for (uint64_t i = 0; i < tails.size(); ++i) {
                uint64_t depth = tails[i];
                if (hrrr[i] == 1) { // Leaf.
                    weights[offsets[depth]++] = 0;
                } else { // Inner node.
                    uint64_t idx =  i - h_rank(i); // Count zeros.
                    weights[offsets[depth]++] = old_weights[idx];
                }
            }
        }
        {
            rmq_type rmq(&weights);
            cout << "rmq_size: " << rmq.size() << endl;
            store_to_cache(rmq, key_weights_rmq, cc, true);
        }
        sdsl::dac_vector<> dac_weights(old_weights);
        store_to_cache(dac_weights, key_weights, cc, true);
    }
    cout << "...DOC_OFFSET" << endl;
    typedef sdsl::hyb_sd_vector<>                      doc_offset_type;
    const auto key_doc_offset = surf::KEY_DOC_OFFSET + std::to_string(LEVELS);
    const auto key_doc_offset_select = surf::KEY_DOC_OFFSET_SELECT + std::to_string(LEVELS);
    if (!cache_file_exists<doc_offset_type>(key_doc_offset, cc)) {
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
        std::vector<uint64_t> sa_offset;
        uint64_t sd_n = 1;
        for (uint64_t i = 1; i <= darray.size(); ++i) {
            end = h_select_1(i) + 1 - i;
            if (start < end) { // Dup lens
                std::vector<uint64_t> dup_set(dup.begin() + start, dup.begin() + end);
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
        store_to_cache(doc_offset, key_doc_offset, cc, true);
        store_to_cache(doc_offset_select, key_doc_offset_select, cc, true);
    }

}
} // namespace surf
