#pragma once

#include <queue>

#include "sdsl/k3_treap.hpp"
#include "surf/topk_heap.hpp"

namespace surf {
namespace k3_treap_algos {

// TODO(niklasb) Ideas for optimization:
// * We use mimimal z order here to order quad tree nodes. We should
//   potentially use minimum d actually stored in the node. This can potentially
//   be sampled and only be stored on each i-th level
// * Use min y coordinate of a node for pruning, because most nodes have small y
// * Handle y differently somehow, because it has very special structure

template <typename t_k3_treap>
// items of result vector are (weight, z)
std::vector<std::pair<uint64_t, uint64_t>>
topk_increasing_z(const t_k3_treap& t, size_t k,
                  uint64_t x_lo, uint64_t x_hi,
                  uint64_t y_lo, uint64_t y_hi) {
    using node = std::tuple<
        uint64_t, // lowest possible z coordinate
        uint8_t, // level (root = maximal level)
        typename t_k3_treap::node_type // the node
        >;
    struct cmp {
        bool operator()(const node& a, const node& b) {
            if (std::get<0>(a) != std::get<0>(b))
                return std::get<0>(a) > std::get<0>(b);
            return std::get<1>(a) < std::get<1>(b);
        }
    };

    // items are (zmin, level, node), increasingly sorted by (zmin, -level)
    // TODO apparently it is better to sort from root to leaf, why?
    // Intuitively it should be better to sort from leaf to root
    std::priority_queue<node, std::vector<node>, cmp> q;

    topk_heap2 result(k);

    if (t.size() == 0)
        return {};

    uint64_t d = 0;

    auto is_valid = [&](const typename t_k3_treap::node_type& v) {
        const auto bb = t.bounding_box(v);
        const auto lo = bb.first;
        const auto hi = bb.second;
        return !(hi[2] < d || v.max_v <= result.lower_bound()
                || hi[0] < x_lo || lo[0] > x_hi
                || hi[1] < y_lo || lo[1] > y_hi);
    };

    auto root = t.root();
    q.emplace(t.bounding_box(root).first[0], root.t, root);
    while (!q.empty()) {
        const auto v = std::get<2>(q.top());
        q.pop();
        if (!is_valid(v)) continue;

        uint64_t x = v.max_p[0], y = v.max_p[1], docid = v.max_p[2];

        if (x_lo <= x && x <= x_hi &&
                y_lo <= y && y <= y_hi &&
                docid >= d)
            result.insert(docid, v.max_v);

        if (t.is_leaf(v)) {
            result.insert(docid, v.max_v);
            d = docid + 1;
        } else {
            for (auto w : t.children(v)) {
                if (!is_valid(w)) continue;
                q.emplace(t.bounding_box(w).first[2], w.t, w);
            }
        }
    }
    return result.sorted_result();
}

}  // namespace k3_treap_algos
}  // namespace surf
