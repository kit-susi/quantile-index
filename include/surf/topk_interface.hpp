#pragma once

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

namespace surf {

// (doc id, weight)
using topk_result = std::pair<uint64_t, double>;
using topk_result_set = std::vector<topk_result>;


template <typename t_token>
struct topk_iterator {
    using token_type = t_token;
    using snippet_type = std::vector<token_type>;

    virtual ~topk_iterator() {}
    virtual topk_result get() const = 0;
    virtual bool done() const = 0;
    virtual void next() = 0;
    virtual snippet_type extract_snippet(const size_t k) const = 0;
};

template <typename t_token>
struct topk_index {
    using token_type = t_token;
    using iter = topk_iterator<t_token>;
    using snippet_type = std::vector<token_type>;
    using intersect_query = std::vector<std::pair<const token_type*, const token_type*>>;

    virtual ~topk_index() {}
    virtual std::unique_ptr<iter> topk(
            size_t k, const token_type* begin, const token_type* end,
            bool multi_occ = false, bool match_only = false) = 0;

    virtual std::unique_ptr<iter> topk_intersect(
            size_t k, const intersect_query& query,
            bool multi_occ = false, bool match_only = false) {
        std::cerr << "intersection not implemented" << std::endl;
        abort();
    }

    void set_debug_stream(std::ostream* debug_stream) {
        m_debug_stream = debug_stream;
    }

    std::ostream* get_debug_stream() {
        return m_debug_stream;
    }

private:
    std::ostream* m_debug_stream = nullptr;
};

template <typename t_alphabet_category>
struct topk_index_by_alphabet {
    using type = topk_index<
        typename std::conditional<
            std::is_same<t_alphabet_category, sdsl::int_alphabet_tag>::value,
            uint64_t,
            char>::type>;
};

template <typename t_token>
class vector_topk_iterator : public topk_iterator<t_token> {
public:
    vector_topk_iterator() = delete;
    explicit vector_topk_iterator(const topk_result_set& results) : m_results(results) {}

    topk_result get() const {
        return m_results[m_index];
    }

    bool done() const override {
        return m_index >= m_results.size();
    }

    void next() override {
        m_index++;
    }

    std::vector<t_token> extract_snippet(const size_t k) const override {
        // TODO implement
        return {};
    }

private:
    size_t m_index = 0;
    const topk_result_set& m_results;
};

template <typename t_token>
std::unique_ptr<topk_iterator<t_token>>
sort_topk_results(topk_result_set* results) {
    std::sort(results->begin(), results->end(),
              [&](const topk_result& a, const topk_result& b) {
                  return std::make_pair(-a.second, a.first) <
                    std::make_pair(-b.second, b.first);
              });
    // remove negative weights (we use this as a workaround inside idx_d to
    // implement multi_occ=true)
    while(!results->empty() && results->back().second < 0) results->pop_back();
    return std::make_unique<vector_topk_iterator<t_token>>(*results);
}

}  // namespace surf
