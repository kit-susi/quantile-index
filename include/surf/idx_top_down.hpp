#pragma once

#include "sdsl/suffix_trees.hpp"
#include "surf/topk_interface.hpp"

namespace surf {

template<typename t_csa,
         typename t_border = sdsl::sd_vector<>,
         typename t_border_rank = typename t_border::rank_1_type,
         typename t_border_select = typename t_border::select_1_type, 
         int LEVELS = 10,
         typename t_tails = sdsl::hyb_sd_vector<>,
         typename t_tails_select = typename t_tails::select_1_type>
class idx_top_down
    : public topk_index_by_alphabet<typename t_csa::alphabet_category>::type {

  public:
    typedef t_csa                                      csa_type;
    typedef t_border                                   border_type;
    typedef t_border_rank                              border_rank_type;
    typedef t_border_select                            border_select_type;
    typedef t_tails                                    tails_type;
    typedef t_tails_select                             tails_select_type;
    using size_type = sdsl::int_vector<>::size_type;
    typedef typename t_csa::alphabet_category          alphabet_category;
    using topk_interface = typename topk_index_by_alphabet<alphabet_category>::type;

  private:
    csa_type           m_csa;
    border_type        m_border;
    border_rank_type   m_border_rank;
    border_select_type m_border_select;
    tails_type         m_tails;
    tails_select_type  m_tails_select;

  public:
    std::unique_ptr<typename topk_interface::iter> topk(
            size_t k,
            const typename topk_interface::token_type* begin,
            const typename topk_interface::token_type* end,
            bool multi_occ = false, bool only_match = false) override {
        std::cerr << "Not implemented yet" << std::endl;
        abort();
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
        load_from_cache(m_tails_select, surf::KEY_TAILS_SELECT + std::to_string(LEVELS), cc, true);
        m_tails_select.set_vector(&m_tails);
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
        written_bytes += m_tails_select.serialize(out, child, "TAILS_SELECT");
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void mem_info()const {
        std::cout << sdsl::size_in_bytes(m_csa) + 
                sdsl::size_in_bytes(m_border) + 
                sdsl::size_in_bytes(m_border_rank) +
                sdsl::size_in_bytes(m_border_select) << ";"; // CSA
        std::cout << sdsl::size_in_bytes(m_tails) + 
                     sdsl::size_in_bytes(m_tails_select) << ";"; // TAILS
    }
};

template<typename t_csa,
         typename t_border,
         typename t_border_rank,
         typename t_border_select,
         int LEVELS,
         typename t_tails,
         typename t_tails_select>
void construct(idx_top_down<t_csa,
               t_border, t_border_rank,
               t_border_select,
               LEVELS,
               t_tails,
               t_tails_select>& idx,
               const std::string&, sdsl::cache_config& cc, uint8_t num_bytes) {
    using t_wtd = WTD_TYPE;
    using cst_type =  sdsl::cst_sct3<t_csa, sdsl::lcp_dac<>, sdsl::bp_support_sada<>, sdsl::bit_vector, sdsl::rank_support_v<>, sdsl::select_support_mcl<>>;
    using node_type = typename cst_type::node_type;
    using tails_type = t_tails;
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
    if (!cache_file_exists(surf::KEY_C, cc)) {
        construct_doc_cnt<t_csa::alphabet_type::int_width>(cc);
        uint64_t doc_cnt = 0;
        load_from_cache(doc_cnt, KEY_DOCCNT, cc);
        cout << "doc_cnt = " << doc_cnt << endl;
        string d_file = cache_file_name(surf::KEY_DARRAY, cc);
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
    {
        cst_type cst;
        t_wtd wtd;
        int_vector<> last_occ;
        load_from_file(cst, cache_file_name<cst_type>(surf::KEY_TMPCST, cc));
        load_from_file(wtd, cache_file_name<t_wtd>(surf::KEY_WTD, cc));
        load_from_file(last_occ, cache_file_name(surf::KEY_C, cc)); 
        bit_vector bv(LEVELS * wtd.size(), 0); // Tails plain bv.
        auto add_lca = [&](uint64_t u, uint64_t v) {
                if (u < wtd.size() && v < wtd.size()) { // Check if exists.
                    node_type uu = cst.select_leaf(u+1);
                    node_type vv = cst.select_leaf(v+1);
                    uint64_t depth = cst.depth(cst.lca(uu,vv));
                    if (depth < LEVELS) {
                        bv[depth*wtd.size() + u] = 1; 
                        bv[depth*wtd.size() + v] = 1; 
                    }
                }
        };
        for (uint64_t pos = 0; pos < wtd.size(); pos++) {
                add_lca(last_occ[pos], pos);
        }
        tails_type tails(bv);
        tails_select_type tails_select(&tails);
        store_to_cache(tails, surf::KEY_TAILS + std::to_string(LEVELS), cc, true); 
        store_to_cache(tails_select, surf::KEY_TAILS_SELECT + std::to_string(LEVELS), cc, true); 
        // Print bv for debug.
        if (bv.size() < 1000) {
            for (size_t i = 0; i < bv.size(); i++) {
                if (i % wtd.size() == 0) cout << endl;
                cout << bv[i];
            } cout << endl;
        }
    }
}
} // namespace surf
