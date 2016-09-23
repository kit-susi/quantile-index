#pragma once

#include "sdsl/suffix_trees.hpp"
#include "surf/topk_interface.hpp"
#include "sdsl/rmq_succinct_sct.hpp"

namespace surf {

using rmq_type = rmq_succinct_sct<false>;
template<typename t_csa,
         typename t_border = sdsl::sd_vector<>,
         typename t_border_rank = typename t_border::rank_1_type,
         typename t_border_select = typename t_border::select_1_type,
         int LEVELS = 10,
         typename t_tails = sdsl::hyb_sd_vector<>,
         typename t_tails_rank = typename t_tails::rank_1_type,
         typename t_tails_select = typename t_tails::select_1_type>
class idx_top_down
    : public topk_index_by_alphabet<typename t_csa::alphabet_category>::type {

public:
    typedef t_csa                                      csa_type;
    typedef t_border                                   border_type;
    typedef t_border_rank                              border_rank_type;
    typedef t_border_select                            border_select_type;
    typedef t_tails                                    tails_type;
    typedef t_tails_rank                               tails_rank_type;
    typedef t_tails_select                             tails_select_type;
    using size_type = sdsl::int_vector<>::size_type;
    typedef typename t_csa::alphabet_category          alphabet_category;
    using topk_interface = typename topk_index_by_alphabet<alphabet_category>::type;

private:
    using token_type = typename topk_interface::token_type;

    csa_type           m_csa;
    border_type        m_border;
    border_rank_type   m_border_rank;
    border_select_type m_border_select;
    tails_type         m_tails;
    tails_rank_type    m_tails_rank;
    tails_select_type  m_tails_select;
    int_vector<>       m_weights;
    rmq_type           m_weights_rmq;

public:

    template <typename t_token>
    class top_down_topk_iterator : public topk_iterator<t_token> {
    private:
        const idx_top_down*m_idx;
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
            bool operator<(const top_down_result& obj) const {
                return weight < obj.weight;
            }
        };
        std::priority_queue<top_down_result> m_intervals;

        void init_interval(top_down_result& interval) {
            interval.tails_id = m_idx->m_weights_rmq(interval.start, interval.end-1);
            assert(interval.tails_id >= interval.start && interval.tails_id < interval.end);
            interval.weight = m_idx->m_weights[interval.tails_id];
            // TODO compute document ids later.
            uint64_t sa_pos = m_idx->m_tails_select(interval.tails_id+1) % m_idx->m_csa.size();
            interval.document = m_idx->m_border_rank(m_idx->m_csa[sa_pos]);
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
                uint64_t offset = idx->m_csa.size();
                std::cerr << "interval size: " << m_ep - m_sp+1<< endl;
                // <= here instead of < ??? TODO
                for (int64_t depth = 0; depth < (end - begin); ++depth) {
                    top_down_result interval;
                    interval.start = idx->m_tails_rank(depth*offset + m_sp);
                    interval.end = idx->m_tails_rank(depth*offset + m_ep+1);
                    std::cerr << interval.start << " " << interval.end << " " << interval.end - interval.start << endl;
                    std::cerr << "with weights: ";
                    for (int i = interval.start; i < interval.end; i++)
                            std::cerr << idx->m_weights[i] << " "; std::cerr << endl;
                    if (interval.start < interval.end) {
                        init_interval(interval);
                        m_intervals.push(interval);     
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
            m_intervals.pop();
            top_down_result interval;
            // Left side.
            interval.start = t.start;
            interval.end = t.tails_id;
            if (interval.start < interval.end) {
                    init_interval(interval);
                    m_intervals.push(interval);
            }
            // Right side.
            interval.start = t.tails_id+1;
            interval.end = t.end;
            if (interval.start < interval.end) {
                    init_interval(interval);
                    m_intervals.push(interval);
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
        load_from_cache(m_tails_rank, surf::KEY_TAILS_RANK + std::to_string(LEVELS), cc, true);
        load_from_cache(m_tails_select, surf::KEY_TAILS_SELECT + std::to_string(LEVELS), cc, true);
        m_tails_rank.set_vector(&m_tails);
        m_tails_select.set_vector(&m_tails);
        load_from_cache(m_weights, surf::KEY_WEIGHTS + std::to_string(LEVELS), cc);
        load_from_cache(m_weights_rmq, surf::KEY_WEIGHTS_RMQ, cc, true);
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
        written_bytes += m_tails_rank.serialize(out, child, "TAILS_RANK");
        written_bytes += m_tails_select.serialize(out, child, "TAILS_SELECT");
        written_bytes += m_weights.serialize(out, child, "WEIGHTS");
        written_bytes += m_weights_rmq.serialize(out, child, "WEIGHTSRMQ");
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void mem_info()const {
        std::cout << sdsl::size_in_bytes(m_csa) +
                  sdsl::size_in_bytes(m_border) +
                  sdsl::size_in_bytes(m_border_rank) +
                  sdsl::size_in_bytes(m_border_select) << ";"; // CSA
        std::cout << sdsl::size_in_bytes(m_tails) +
                  sdsl::size_in_bytes(m_tails_rank) <<
                  sdsl::size_in_bytes(m_tails_select) << ";"; // TAILS
        std::cout << sdsl::size_in_bytes(m_weights) +
                sdsl::size_in_bytes(m_weights_rmq)<< ";"; // WEIGHTS
    }
};

template<typename t_csa,
         typename t_border,
         typename t_border_rank,
         typename t_border_select,
         int LEVELS,
         typename t_tails,
         typename t_tails_rank,
         typename t_tails_select>
void construct(idx_top_down<t_csa,
               t_border, t_border_rank,
               t_border_select,
               LEVELS,
               t_tails,
               t_tails_rank,
               t_tails_select>& idx,
               const std::string&, sdsl::cache_config& cc, uint8_t num_bytes) {
    using t_wtd = WTD_TYPE;
    using cst_type =  sdsl::cst_sct3<t_csa, sdsl::lcp_dac<>, sdsl::bp_support_sada<>, sdsl::bit_vector, sdsl::rank_support_v<>, sdsl::select_support_mcl<>>;
    using node_type = typename cst_type::node_type;
    using tails_type = t_tails;
    using tails_rank_type = t_tails_rank;
    using tails_select_type = t_tails_select;

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
    cout << "...DOC_BORDER" << endl;
    if (!cache_file_exists<t_border>(surf::KEY_DOCBORDER, cc) or
            !cache_file_exists<t_border_rank>(surf::KEY_DOCBORDER_RANK, cc) or
            !cache_file_exists<t_border_select>(surf::KEY_DOCBORDER_SELECT, cc)) {
        bit_vector doc_border;
        load_from_cache(doc_border, surf::KEY_DOCBORDER, cc);
        t_border sd_doc_border(doc_border);
        store_to_cache(sd_doc_border, surf::KEY_DOCBORDER, cc, true);
        t_border_rank doc_border_rank(&sd_doc_border);
        store_to_cache(doc_border_rank, surf::KEY_DOCBORDER_RANK, cc, true);
        t_border_select doc_border_select(&sd_doc_border);
        store_to_cache(doc_border_select, surf::KEY_DOCBORDER_SELECT, cc, true);
    }
    cout << "...TMPCST" << endl;
    {
        if (!cache_file_exists(conf::KEY_SA, cc)) {
            construct_sa<t_csa::alphabet_type::int_width>(cc);
        }
        register_cache_file(conf::KEY_SA, cc);
        if (!cache_file_exists(conf::KEY_LCP, cc)) {
            if (t_csa::alphabet_type::int_width == 8) {
                cout << "byte lcp construct" << endl;
                construct_lcp_semi_extern_PHI(cc);
            } else {
                cout << "int lcp construct" << endl;
                construct_lcp_PHI<t_csa::alphabet_type::int_width>(cc);
            }
        }
        register_cache_file(conf::KEY_LCP, cc);
        if (!cache_file_exists<cst_type>(KEY_TMPCST, cc)) {
            auto event = memory_monitor::event("construct cst");
            cst_type cst = cst_type(cc);
            store_to_file(cst, cache_file_name<cst_type>(surf::KEY_TMPCST, cc));
        }
    }
    cout << "...C" << endl;
    string d_file = cache_file_name(surf::KEY_DARRAY, cc);
    if (!cache_file_exists(surf::KEY_C, cc)) {
        construct_doc_cnt<t_csa::alphabet_type::int_width>(cc);
        uint64_t doc_cnt = 0;
        load_from_cache(doc_cnt, KEY_DOCCNT, cc);
        cout << "doc_cnt = " << doc_cnt << endl;
        int_vector_buffer<> D(d_file);
        cout << "n=" << D.size() << endl;
        int_vector<> C(D.size(), 0, bits::hi(D.size()) + 1);
        int_vector<> last_occ(doc_cnt + 1, D.size(), bits::hi(D.size()) + 1);
        for (size_t i = 0; i < D.size(); ++i) {
            uint64_t d = D[i];
            C[i] = last_occ[d];
            last_occ[d] = i;
        }
        util::bit_compress(C);
        store_to_file(C, cache_file_name(surf::KEY_C, cc));
    }
    cout << "...TAILS" << endl;
    const string key_tails = surf::KEY_TAILS + std::to_string(LEVELS);
    const string key_tails_rank = surf::KEY_TAILS_RANK + std::to_string(LEVELS);
    const string key_tails_select = surf::KEY_TAILS_SELECT + std::to_string(LEVELS);
    cst_type cst;
    if (!cache_file_exists<tails_type>(key_tails, cc)) {
        int_vector<> last_occ;
        // Loading files.
        int_vector_buffer<> D(d_file);
        load_from_file(cst, cache_file_name<cst_type>(surf::KEY_TMPCST, cc));
        load_from_file(last_occ, cache_file_name(surf::KEY_C, cc));
        bit_vector bv(LEVELS * D.size(), 0); // Tails plain bv.

        auto add_lca = [&](uint64_t u, uint64_t v) {
            if (u < D.size() && v < D.size()) { // Check if exists.
                node_type uu = cst.select_leaf(u + 1);
                node_type vv = cst.select_leaf(v + 1);
                uint64_t depth = cst.depth(cst.lca(uu, vv));
                if (depth < LEVELS) {
                    //bv[depth * D.size() + u] = 1;
                    bv[depth * D.size() + v] = 1;
                }
            }
        };
        for (uint64_t pos = 0; pos < D.size(); pos++) {
            add_lca(last_occ[pos], pos);
        }
        {   // Add missing ones in bv.
            cout << "Add missing ones." << endl;
            std::unordered_set<uint64_t> doc_set;
            auto add_missing_ones = [&](node_type v) {
                uint64_t depth = cst.depth(v);
                // Find all doc at node v.
                for (size_t i = v.i; i <= v.j; ++i) {
                    if (bv[depth * D.size() + i])
                        doc_set.insert(D[i]);
                }
                // Mark first occurence of docs.
                size_t i = v.i;
                while (!doc_set.empty()) {
                    if (doc_set.count(D[i]) == 1) {
                        assert(!bv[depth * D.size() + i]);
                        bv[depth * D.size() + i] = 1;
                        doc_set.erase(D[i]);
                    }
                    ++i;
                }
            };
            // DFS.
            std::stack<node_type> s;
            s.push(cst.root());
            while (!s.empty()) {
                auto v = s.top();
                s.pop();
                for (const auto u : cst.children(v)) {
                    if (cst.depth(u) < LEVELS)
                        s.push(u);
                }
                add_missing_ones(v);
            }
        }
        {   // Building and writing tails bv in correct format.
            tails_type tails(bv);
            tails_rank_type tails_rank(&tails);
            tails_select_type tails_select(&tails);
            store_to_cache(tails, key_tails, cc, true);
            store_to_cache(tails_rank, key_tails_rank, cc, true);
            store_to_cache(tails_select, key_tails_select, cc, true);
        }
    }

    const string key_next_occ = surf::KEY_NEXT_OCC + std::to_string(LEVELS);
    cout << "...next_occ" << endl;
    if (!cache_file_exists(key_next_occ, cc)) {
        uint64_t doc_cnt = 0;
        load_from_cache(doc_cnt, KEY_DOCCNT, cc);
        int_vector<> D;
        tails_type tails;
        tails_rank_type tails_rank;
        tails_select_type tails_select;
        load_from_file(D, d_file);
        int_vector_buffer<> next_occ(cache_file_name(key_next_occ, cc), std::ios::out, 1 << 20,
                                     bits::hi(D.size() * LEVELS+ 1) + 1);
        load_from_file(tails, cache_file_name<tails_type>(key_tails, cc));
        load_from_file(tails_rank, cache_file_name<tails_rank_type>(key_tails_rank, cc));
        load_from_file(tails_select, cache_file_name<tails_select_type>(key_tails_select, cc));
        tails_rank.set_vector(&tails);
        tails_select.set_vector(&tails);
        std::vector<uint64_t> cur_next_occ(doc_cnt, D.size()*LEVELS);
        uint64_t d = 0;
        uint64_t num_arrows = tails_rank(LEVELS * D.size());
        cout << "Num arrows:" << num_arrows << " " << D.size() << endl;
        for (uint64_t i = num_arrows; i > 0; --i) {
            d = D[tails_select(i) % D.size()];
            next_occ[i - 1] = cur_next_occ[d];
            cur_next_occ[d] = i; // one based!!
        }
    }
    cout << "...weights" << endl;
    const string key_weights = surf::KEY_WEIGHTS + std::to_string(LEVELS);
    if (!cache_file_exists(key_weights, cc)) {
        using std::tuple;
        using std::make_tuple;
        using std::get;
        uint64_t output_count = 0;
        t_wtd wtd;
        cst_type cst;
        tails_type tails;
        tails_rank_type tails_rank;
        tails_select_type tails_select;
        int_vector<> next_occ;
        load_from_file(wtd, cache_file_name<t_wtd>(surf::KEY_WTD, cc));
        load_from_file(cst, cache_file_name<cst_type>(surf::KEY_TMPCST, cc));
        load_from_file(tails, cache_file_name<tails_type>(key_tails, cc));
        load_from_file(tails_rank, cache_file_name<tails_rank_type>(key_tails_rank, cc));
        load_from_file(tails_select, cache_file_name<tails_select_type>(key_tails_select, cc));
        load_from_file(next_occ, cache_file_name(key_next_occ, cc));
        tails_rank.set_vector(&tails);
        tails_select.set_vector(&tails);

        int_vector_buffer<> weights_buf(cache_file_name(key_weights, cc),
                        std::ios::out, 1 << 20, bits::hi(wtd.size()) + 1);
        auto calc_weights_for_node = [&](node_type v, uint64_t depth) {
            uint64_t offset = depth * wtd.size();
            uint64_t s = tails_rank(v.i + offset);         // inclusive.
            uint64_t e = tails_rank(v.j + 1 + offset); // exclusive.
            for (uint64_t i = s; i < e; ++i) {
                uint64_t sa_start = tails_select(i + 1) - offset; // inclusive.
                uint64_t sa_end = v.j + 1; // exclusive.
                uint64_t d = wtd[sa_start];
                if (next_occ[i] <= e) // next_occ is one based.
                    sa_end = std::min(wtd.size(), tails_select(next_occ[i]) - offset);
                assert(output_count == i);
                uint64_t w = wtd.rank(sa_end, d) - wtd.rank(sa_start, d);
                weights_buf[output_count++] = w;
                //cout << depth << " " << d << " " << w << endl;
            }
        };
        using pq_type = tuple<uint64_t, node_type>;
        std::priority_queue<pq_type, std::vector<pq_type>, std::greater<pq_type>> q;
        q.push(make_tuple(0, cst.root()));
        while (!q.empty()) {
            auto v = get<1>(q.top());
            uint64_t d = get<0>(q.top());
            q.pop();
            for (const auto u : cst.children(v)) {
                uint64_t depth = cst.depth(u);
                if (depth < LEVELS)
                    q.push(make_tuple(depth, u));
            }
            calc_weights_for_node(v, d);
        }
    }
    if (!cache_file_exists<rmq_type>(surf::KEY_WEIGHTS_RMQ, cc)) {
        cout << "Construct rmq." << endl;
        int_vector<> weights;
        load_from_file(weights, cache_file_name(key_weights, cc));
        rmq_type rmq(&weights);
        cout << "rmq_size: " << rmq.size() << endl;
        store_to_cache(rmq, surf::KEY_WEIGHTS_RMQ, cc, true);
    }
}
} // namespace surf
