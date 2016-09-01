#pragma once

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

namespace surf {

using topk_result = std::pair<uint64_t, double>;
using topk_result_set = std::vector<topk_result>;
using topk_snippet = std::vector<uint64_t>;

class topk_iterator {
public:
    virtual topk_result get() const = 0;
    virtual bool done() const = 0;
    virtual void next() = 0;
    virtual topk_snippet extract_snippet(const size_t k) const = 0;
};

class vector_topk_iterator : public topk_iterator {
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

    topk_snippet extract_snippet(const size_t k) const override {
        // TODO implement
        return topk_snippet();
    }

private:
    size_t m_index = 0;
    const topk_result_set& m_results;
};

std::unique_ptr<topk_iterator>
sort_topk_results(topk_result_set* results) {
    std::sort(results->begin(), results->end(),
              [&](const topk_result& a, const topk_result& b) {
                  return std::make_pair(-a.second, a.first) <
                    std::make_pair(-b.second, b.first);
              });
    return std::make_unique<vector_topk_iterator>(*results);
}

}  // namespace surf
