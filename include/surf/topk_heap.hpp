#pragma once
#include <map>
#include <vector>

namespace surf {

class topk_heap {
private:
    size_t m_k, m_result_size;
    std::vector<std::pair<uint64_t, uint64_t>> m_result;
    std::map<uint64_t, size_t> m_pos;

    const uint64_t inf = std::numeric_limits<uint64_t>::max();

    void swap(size_t i, size_t j) {
        std::swap(m_result[i], m_result[j]);
        //m_pos[m_result[i].second] = i;
        //m_pos[m_result[j].second] = j;
    }

    void sift_down(size_t i) {
        uint64_t w = m_result[i].first;
        uint64_t lw = 2*i + 1 < m_result_size ? m_result[2*i + 1].first : inf;
        uint64_t rw = 2*i + 2 < m_result_size ? m_result[2*i + 2].first : inf;
        if (w <= lw && w <= rw)
            return;
        size_t swap_index = 2*i + 1 + (rw < lw);
        swap(i, swap_index);
        sift_down(swap_index);
    }

    void sift_up(size_t i) {
        if (i == 0)
            return;
        size_t parent = (i - 1) / 2;
        if (m_result[i].first >= m_result[parent].first)
            return;
        swap(i, parent);
        sift_up(parent);
    }

    void check_invariants() {
        return;
        for (size_t i = 0; i < m_result_size; ++i) {
            if (2*i + 1 < m_result_size && m_result[i].first > m_result[2*i + 1].first) {
                std::cout << "heap1 violated" << std::endl;
                std::cout << *(volatile int*)0 << std::endl;
            }
            if (2*i + 2 < m_result_size && m_result[i].first > m_result[2*i + 2].first) {
                std::cout << "heap2 violated" << std::endl;
                std::cout << *(volatile int*)0 << std::endl;
            }
            if (m_pos.count(m_result[i].second) && m_pos[m_result[i].second] != i) {
                std::cout << "m_pos broken" << std::endl;
                std::cout << *(volatile int*)0 << std::endl;
            }
        }
        if (m_pos.size() != m_result_size) {
            std::cout << "m_pos too large: " << m_pos.size() << " " << m_result_size << std::endl;
            std::cout << *(volatile int*)0 << std::endl;
        }
    }

public:
    uint64_t updates = 0, total = 0;

    topk_heap(size_t k)
        : m_k(k), m_result_size(), m_result(k + 1) {}

    uint64_t lower_bound() const {
        return m_result_size == m_k ? m_result[0].first : 0;
    }

    void insert(uint64_t docid, uint64_t weight) {
        if (m_result_size == m_k) {
            // replace min
            //m_pos.erase(m_result[0].second);
            m_result[0] = std::make_pair(weight, docid);
            //m_pos[docid] = 0;
            sift_down(0);
        } else {
            // insert
            m_result[m_result_size++] = std::make_pair(weight, docid);
            //m_pos[m_result[m_result_size - 1].second] = m_result_size - 1;
            sift_up(m_result_size - 1);
        }
        check_invariants();
    }

    void insert_or_update(uint64_t docid, uint64_t weight) {
        if (weight < lower_bound())
            return;
        total++;
        //*
        for (size_t i = 0; i < m_result_size; ++i) {
            if (docid == m_result[i].second) {
                if (weight > m_result[i].first) {
                    m_result[i].first = weight;
                    sift_down(i);
                }
                updates++;
                check_invariants();
                return;
            }
        }
        //*/
        /*
        auto it = m_pos.find(docid);
        if (it != m_pos.end()) {
            size_t i = it->second;
            if (weight > m_result[i].first) {
                m_result[i].first = weight;
                sift_down(i);
            }
            check_invariants();
            return;
        }
        //*/

        insert(docid, weight);
    }

    std::vector<std::pair<uint64_t, uint64_t>> sorted_result() {
        m_result.resize(m_result_size);
        std::sort(m_result.begin(), m_result.end(), [&](const auto& a, const auto& b) {
            return a.first > b.first;
        });
        return m_result;
    }
};


class topk_heap2 {
private:
    size_t m_k, m_result_size;
    // (weight, docid)
    std::vector<std::pair<uint64_t, uint64_t>> m_result;
    const uint64_t inf = std::numeric_limits<uint64_t>::max();

    struct cmp {
        bool operator()(const std::pair<uint64_t, uint64_t>& a,
                        const std::pair<uint64_t, uint64_t>& b) const {
            return a.first > b.first;
        }
    };

public:
    topk_heap2(size_t k) : m_k(k), m_result_size(0), m_result(k + 1) {}

    uint64_t lower_bound() const {
        return m_result_size == m_k ? m_result[0].first : 0;
    }

    void insert(uint64_t docid, uint64_t weight) {
        m_result[m_result_size] = {weight, docid};
        if (m_result_size == m_k)
            std::pop_heap(m_result.data(), m_result.data() + m_result_size + 1, cmp());
        else
            std::push_heap(m_result.data(), m_result.data() + (++m_result_size), cmp());
    }

    std::vector<std::pair<uint64_t, uint64_t>> sorted_result() {
        m_result.resize(m_result_size);
        std::sort(m_result.begin(), m_result.end(), [&](const auto& a, const auto& b) {
            return a.first > b.first;
        });
        return m_result;
    }
};

}  // namespace surf
