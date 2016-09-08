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
class topk_iterator {
public:
    virtual topk_result get() const = 0;
    virtual bool done() const = 0;
    virtual void next() = 0;
    virtual std::vector<t_token> extract_snippet(const size_t k) const = 0;
};

template <typename t_token>
class topk_index {
public:
    virtual std::unique_ptr<topk_iterator<t_token>> topk(
            size_t k, const t_token* begin, const t_token* end,
            bool multi_occ = false, bool only_match = false);
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
    return std::make_unique<vector_topk_iterator<t_token>>(*results);
}

}  // namespace surf
