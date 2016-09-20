#pragma once

#include "sdsl/suffix_trees.hpp"
#include "surf/topk_interface.hpp"

namespace surf {

template<typename t_csa,
         typename t_border = sdsl::sd_vector<>,
         typename t_border_rank = typename t_border::rank_1_type,
         typename t_border_select = typename t_border::select_1_type>
class idx_top_down
    : public topk_index_by_alphabet<typename t_csa::alphabet_category>::type {

  public:
    typedef t_csa                                      csa_type;
    typedef t_border                                   border_type;
    typedef t_border_rank                              border_rank_type;
    typedef t_border_select                            border_select_type;
    using size_type = sdsl::int_vector<>::size_type;
    typedef typename t_csa::alphabet_category          alphabet_category;
    using topk_interface = typename topk_index_by_alphabet<alphabet_category>::type;

  private:
    csa_type           m_csa;
    border_type        m_border;
    border_rank_type   m_border_rank;
    border_select_type m_border_select;

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
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void mem_info()const {
        std::cout << sdsl::size_in_bytes(m_csa) + 
                sdsl::size_in_bytes(m_border) + 
                sdsl::size_in_bytes(m_border_rank) +
                sdsl::size_in_bytes(m_border_select)<< ";"; // CSA
    }
};

template<typename t_csa>
void construct(idx_top_down<t_csa>& idx,
               const std::string&, sdsl::cache_config& cc, uint8_t num_bytes) {
        std::cerr << "Not implemented yet" << std::endl;
        abort();
}
} // namespace surf
