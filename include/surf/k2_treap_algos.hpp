#pragma once

#include <queue>

#include "sdsl/k2_treap.hpp"
#include "surf/topk_heap.hpp"

namespace surf {
namespace k2_treap_algos {

template <typename t_k2_treap>
// items of result vector are (weight, docid)
std::vector<std::pair<uint64_t, uint64_t>>
topk_increasing_y(const t_k2_treap& t, size_t k, uint64_t x_lo, uint64_t x_hi) {
    using node = std::tuple<
        uint64_t, // south border
        uint8_t, // level (root = maximal level)
        typename t_k2_treap::node_type, // the node
        bool // parent doc id (if inserted, otherwise -1)
        >;
    struct cmp {
        bool operator()(const node& a, const node& b) {
            if (std::get<0>(a) != std::get<0>(b))
                return std::get<0>(a) > std::get<0>(b);
            return std::get<1>(a) < std::get<1>(b);
        }
    };

    // items are (south border, level), increasingly sorted by (south border, -level)
    // TODO apparently it is better to sort from root to leaf, why?
    // Intuitively it should be better to sort from leaf to root
    std::priority_queue<node, std::vector<node>, cmp> q;

    topk_heap result(k);

    if (t.size() == 0)
        return {};

    auto root = t.root();
    q.emplace(root.south(t), root.t, root, true);
    uint64_t d = 0;
    size_t max_q_size = 0;
    uint64_t dequeued  =0;
    while (!q.empty()) {
        max_q_size = std::max(max_q_size, q.size());
        auto v = std::get<2>(q.top());
        bool consider_inserting_max = std::get<3>(q.top());
        q.pop();
        if (v.north(t) < d || v.max_v <= result.lower_bound()
                || v.east(t) < x_lo || v.west(t) > x_hi)
            continue;
        dequeued++;
        uint64_t x = real(v.max_p), docid = imag(v.max_p);
        bool max_in_range = x_lo <= x && x <= x_hi;
        if (consider_inserting_max && max_in_range && docid >= d) {
            result.insert_or_update(docid, v.max_v);
        }
        if (t.is_leaf(v)) {
            result.insert_or_update(docid, v.max_v);
            d = docid + 1;
        } else {
            for (auto w : t.children(v)) {
                if (w.north(t) < d || w.max_v <= result.lower_bound()
                        || w.east(t) < x_lo || w.west(t) > x_hi)
                    continue;
                auto child_x = real(w.max_p);
                bool child_max_in_range = x_lo <= child_x && child_x <= x_hi;
                q.emplace(w.south(t), w.t, w,
                        docid != imag(w.max_p) || max_in_range != child_max_in_range);
            }
        }
    }
    //std::cerr << "nodes dequeued = " << dequeued << std::endl;
    //std::cerr << result.updates << "/" << result.total << std::endl;
    return result.sorted_result();
}

template <typename t_k2_treap>
// items of result vector are (weight, docid)
std::vector<std::pair<uint64_t, uint64_t>>
topk_increasing_y2(const t_k2_treap& t, size_t k, uint64_t x_lo, uint64_t x_hi) {
    topk_heap result(k);
    if (t.size() == 0)
        return {};

    auto it = top_k(t,
                    {x_lo, 0}, {x_hi, 100},
                    result.lower_bound());
    //std::cerr << "0-100 = " << (!!it) << std::endl;

    uint64_t d = 0;
    uint64_t top = t.root().north(t);
    while (d <= top) {
        //std::cerr << "looking for " << d << std::endl;
        uint64_t hi = 1;
        bool found = 0;
        for (;;) {
            assert(d + hi - 1 <= top);
            //std::cerr << "  hi = " << hi << std::endl;
            auto it = top_k(t,
                            {x_lo, d}, {x_hi, d + hi - 1},
                            result.lower_bound());
            if (it) {
                //std::cerr << "  yep." << std::endl;
                found = true;
                break;
            }
            //std::cerr << "  nope." << std::endl;

            if (d + hi - 1 == top)
                break;
            hi = std::min(hi << 1, top + 1 - d);
        }
        if (!found)
            break;

        uint64_t lo = 1;  // max(1, hi >> 1)
        while (lo < hi) {
            uint64_t mid = (lo + hi)/2;
            auto it = top_k(t,
                            {x_lo, d}, {x_hi, d + mid - 1},
                            result.lower_bound());
            //std::cerr << "mid=" << mid << " "  << (!!it) << std::endl;
            if (it) hi = mid;
            else lo = mid + 1;
        }
        //std::cerr << "lo=" << lo << " hi="  << hi << std::endl;

        assert(lo == hi);
        auto it = top_k(t, {x_lo, d + lo - 1}, {x_hi, d + lo - 1},
                        result.lower_bound());
        assert(it);
        auto point = *it;
        auto docid = imag(point.first);
        //std::cerr << "found docid " << docid << std::endl;
        result.insert(docid, point.second);
        d = docid + 1;
    }

    return result.sorted_result();
}

}  // namespace k2_treap_algos
}  // namespace surf
