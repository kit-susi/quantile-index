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

using xy_point = std::array<uint64_t, 2>;
using xy_range = std::pair<xy_point, xy_point>;

template <typename t_k3_treap>
bool find_successor(const t_k3_treap& t,
                    k3_treap_ns::point_type range_lo,
                    k3_treap_ns::point_type range_hi,
                    uint64_t weight_lower_bound,
                    k3_treap_ns::point_type* result,
                    uint64_t* result_weight) {
    if (t.size() == 0) return false;

    using node = std::tuple<
        bool, // true = the node itself, false = its max
        uint64_t, // lowest possible z coordinate
        uint8_t, // level (root = maximal level)
        typename t_k3_treap::node_type // the node
        >;
    struct cmp {
        bool operator()(const node& a, const node& b) {
            if (std::get<1>(a) != std::get<1>(b))
                return std::get<1>(a) > std::get<1>(b);
            return std::get<2>(a) < std::get<2>(b);
        }
    };

    //bool dbg = range_lo[2] == 185;
    //if (dbg)
        //std::cerr << "DEBUG ON " << std::endl;

    // items are (zmin, level, node), increasingly sorted by (zmin, -level)
    // TODO apparently it is better to sort from root to leaf, why?
    // Intuitively it should be better to sort from leaf to root
    std::priority_queue<node, std::vector<node>, cmp> q;

    auto is_valid = [&](const typename t_k3_treap::node_type& v) {
        const auto bb = t.bounding_box(v);
        const auto lo = bb.first;
        const auto hi = bb.second;
        return !(v.max_v <= weight_lower_bound
                || hi[0] < range_lo[0] || lo[0] > range_hi[0]
                || hi[1] < range_lo[1] || lo[1] > range_hi[1]
                || hi[2] < range_lo[2] || lo[2] > range_hi[2]
                );
    };

    auto in_range = [&](const typename k3_treap_ns::point_type& p) {
        return range_lo[0] <= p[0] && p[0] <= range_hi[0]
            && range_lo[1] <= p[1] && p[1] <= range_hi[1]
            && range_lo[2] <= p[2] && p[2] <= range_hi[2];
    };

    auto root = t.root();
    q.emplace(true, t.bounding_box(root).first[0], root.t, root);
    //std::cerr << "lo = " << range_lo[0] << " " << range_lo[1] << " " << range_lo[2] << std::endl;
    //std::cerr << "hi = " << range_hi[0] << " " << range_hi[1] << " " << range_hi[2] << std::endl;

    while (!q.empty()) {
        bool first_time = std::get<0>(q.top());
        const auto v = std::get<3>(q.top());
        q.pop();

        auto bb = t.bounding_box(v);
        //std::cerr << "============" << std::endl;
        //std::cerr << "  "<<bb.first[0] << " " << bb.first[1] << " " << bb.first[2] << endl;
        //std::cerr << "  "<<bb.second[0] << " " << bb.second[1] << " " << bb.second[2] << endl;
        //std::cerr << "  " << v.max_p[0] << " " << v.max_p[1] << " " << v.max_p[2] << " " << in_range(v.max_p) << endl;
        if (!is_valid(v)) continue;
        //std::cerr << "  " << in_range(v.max_p) << std::endl;

        if ((!first_time || t.is_leaf(v)) && in_range(v.max_p)) {
            *result = v.max_p;
            *result_weight = v.max_v;
            return true;
        }

        // queue children
        for (auto w : t.children(v)) {
            if (!is_valid(w)) continue;
            q.emplace(true, t.bounding_box(w).first[2], w.t, w);
        }

        // re-queue the maximum
        if (first_time && !t.is_leaf(v))
            q.emplace(false, v.max_p[2], 0, v);
    }
    return false;
}

template <typename t_k3_treap>
bool find_successor2(const t_k3_treap& t,
                    k3_treap_ns::point_type range_lo,
                    k3_treap_ns::point_type range_hi,
                    uint64_t weight_lower_bound,
                    k3_treap_ns::point_type* result,
                    uint64_t* result_weight) {
    if (t.size() == 0) return false;

    auto is_valid = [&](const typename t_k3_treap::node_type& v) {
        const auto bb = t.bounding_box(v);
        const auto lo = bb.first;
        const auto hi = bb.second;
        return !(v.max_v <= weight_lower_bound
                || hi[0] < range_lo[0] || lo[0] > range_hi[0]
                || hi[1] < range_lo[1] || lo[1] > range_hi[1]
                || hi[2] < range_lo[2] || lo[2] > range_hi[2]
                );
    };

    auto in_range = [&](const typename k3_treap_ns::point_type& p) {
        return range_lo[0] <= p[0] && p[0] <= range_hi[0]
            && range_lo[1] <= p[1] && p[1] <= range_hi[1]
            && range_lo[2] <= p[2] && p[2] <= range_hi[2];
    };

    std::queue<typename t_k3_treap::node_type> q;
    q.emplace(t.root());

    uint64_t min_z = std::numeric_limits<uint64_t>::max();
    bool found = false;
    while (!q.empty()) {
        const auto v = q.front();
        q.pop();
        if (!is_valid(v)) continue;
        if (v.max_v > weight_lower_bound && in_range(v.max_p) && v.max_p[2] < min_z) {
            min_z = v.max_p[2];
            found = true;
            *result = v.max_p;
            *result_weight = v.max_v;
        }
        for (auto w : t.children(v))
            q.emplace(w);
    }
    return found;
}

template <typename t_k3_treap, typename t_weight, typename t_weight_reduce>
// items of result vector are (weight, z)
std::vector<std::pair<uint64_t, uint64_t>>
topk_intersect(
        const t_k3_treap& t, size_t k,
        std::vector<xy_range>& ranges,
        t_weight weight_init,
        t_weight_reduce weight_reduce) {
    if (!t.size()) return {};
    uint64_t d = 0;

    topk_heap2 result(k);
    const auto inf = std::numeric_limits<uint64_t>::max();

    for (;;) {
        std::cerr << "d  =" << d << std::endl;
        std::vector<k3_treap_ns::point_type> successors(ranges.size());
        uint64_t max_succ = 0;
        bool shared_result = true;
        t_weight combined_weight = weight_init;
        for (size_t i = 0; i < ranges.size(); ++i) {
            uint64_t succ_weight;
            if (!find_successor2(t,
                        {ranges[i].first[0], ranges[i].first[1], d},
                        {ranges[i].second[0], ranges[i].second[1], inf},
                        result.lower_bound(),
                        &successors[i],
                        &succ_weight)) {
                std::cerr << "no successor in range " << i << std::endl;
                // can't find successor in range i, so we are finished
                return result.sorted_result();
            }
            uint64_t succ = successors[i][2];
            std::cerr << "  succ = " << succ << std::endl;
            if (succ != max_succ && i > 0)
                shared_result = false;
            if (succ > max_succ || i == 0)
                max_succ = succ;
            combined_weight = weight_reduce(combined_weight, succ_weight);
        }
        std::cerr << "max_succ = " << max_succ << "  shared = " << shared_result << std::endl;
        if (shared_result) {
            result.insert(max_succ, combined_weight);
            d = max_succ + 1;
        } else {
            d = max_succ;
        }
    }

    return result.sorted_result();
}

template <typename t_k3_treap, typename t_weight, typename t_weight_reduce>
// items of result vector are (weight, z)
std::vector<std::pair<uint64_t, uint64_t>>
topk_intersect2(
        const t_k3_treap& t, size_t k,
        std::vector<xy_range>& ranges,
        t_weight weight_init,
        t_weight_reduce weight_reduce) {
    const auto inf = std::numeric_limits<uint64_t>::max();
    if (t.size() == 0)
        return {};

    topk_heap2 result(k);
    uint64_t d_lo = 0;

    auto max_w = [&](uint64_t d_lo, uint64_t d_hi) {
        uint64_t w = weight_init;
        for (const auto& range: ranges) {
            auto it = k3_treap_ns::top_k(t,
                    {range.first[0], range.first[1], d_lo},
                    {range.second[0], range.second[1], d_hi});
            if (!it) return uint64_t{0};
            w = weight_reduce(w, std::get<1>(*it));
        }
        return w;
    };

    auto eligible_d_hi = [&](uint64_t d_hi) {
        uint64_t w = max_w(d_lo, d_hi);
        return w != -1 && w > result.lower_bound();
    };

    while (eligible_d_hi(inf)) {
        uint64_t window_size = 1;
        while (!eligible_d_hi(d_lo + window_size - 1))
            window_size <<= 1;

        uint64_t lo = d_lo, hi = d_lo + window_size - 1;
        while (lo < hi) {
            uint64_t mid = (lo + hi) / 2;
            if (eligible_d_hi(mid))
                hi = mid;
            else
                lo = mid + 1;
        }

        uint64_t w = max_w(lo, lo);
        if (w != -1)
            result.insert(lo, w);
        d_lo = lo + 1;
    }

    return result.sorted_result();
}


uint64_t total=0;
template <typename t_k3_treap>
struct range {
    const t_k3_treap& m_treap;
    k3_treap_ns::point_type m_lo, m_hi;
    k3_treap_ns::top_k_iterator<t_k3_treap> m_top_iter;

    range(const t_k3_treap& treap,
            const k3_treap_ns::point_type& lo,
            const k3_treap_ns::point_type& hi)
        : m_treap(treap), m_lo(lo), m_hi(hi) {}
    range(const range& other)
        : m_treap(other.m_treap), m_lo(other.m_lo), m_hi(other.m_hi),
          m_top_iter(other.m_top_iter) {}

    range& operator=(const range& other) {
        m_lo = other.m_lo;
        m_hi = other.m_hi;
        m_top_iter = other.m_top_iter;
        return *this;
    }

    std::pair<range, range> split() const {
        uint64_t mid_z = (m_lo[2] + m_hi[2])/2;
        return std::make_pair(
            range<t_k3_treap>(m_treap, m_lo, {m_hi[0], m_hi[1], mid_z}),
            range<t_k3_treap>(m_treap, {m_lo[0], m_lo[1], mid_z + 1}, m_hi)
        );
    }

    uint64_t split_z() const {
        uint64_t mid_z = (m_lo[2] + m_hi[2])/2;
        return mid_z;
    }

    bool find_max(k3_treap_ns::point_type* node, uint64_t* w) {
        if (!m_top_iter.is_initialized()) {
            total++;
            m_top_iter = k3_treap_ns::top_k(m_treap, m_lo, m_hi);
        }
        if (m_top_iter) {
            if (node)
                *node = std::get<0>(*m_top_iter);
            if (w)
                *w = std::get<1>(*m_top_iter);
            return true;
        }
        return false;
    }

    bool is_leaf() {
        if (!m_top_iter.is_initialized()) {
            total++;
            m_top_iter = k3_treap_ns::top_k(m_treap, m_lo, m_hi);
        }

        if (!m_top_iter) abort();
        auto it = m_top_iter;
        ++it;
        return !it;
    }

    bool empty() {
        return !find_max(nullptr, nullptr);
    }

    uint64_t get_max() {
        uint64_t w;
        if (!find_max(nullptr, &w))
            abort();
        return w;
    }

    uint64_t get_max_z() {
        k3_treap_ns::point_type p;
        if (!find_max(&p, nullptr))
            abort();
        return p[2];
    }
};

template <typename t_k3_treap>
// items of result vector are (weight, z)
std::vector<std::pair<uint64_t, uint64_t>>
topk_intersect3(const t_k3_treap& t, size_t k,
                std::vector<xy_range>& ranges,
                uint64_t weight_lower_bound = 0) {
    if (ranges.empty()) abort();
    if (!t.size()) return {};

    const uint64_t inf = std::numeric_limits<uint64_t>::max();
    topk_heap2 result(k);
    auto get_cur_min = [&]() { return result.lower_bound(); };

    std::vector<std::vector<range<t_k3_treap>>> stacks(ranges.size());
    uint64_t d_max = t.bounding_box(t.root()).second[2];
    //uint64_t d_max = 255;

    std::vector<range<t_k3_treap>> vt;
    size_t U = 0;
    size_t count=0;
    for (const auto& r : ranges) {
        //std::cerr << r.first[0] << " " << r.first[1] << " - " << r.second[0] << " " << r.second[1] << std::endl;
        range<t_k3_treap> k3r(t,
                {r.first[0], r.first[1], 0},
                {r.second[0], r.second[1], d_max});
        if (k3r.empty())
            return {};
        count++;
        vt.push_back(k3r);
        U += k3r.get_max() + 1;
    }

    uint64_t d = 0;
    uint64_t L = weight_lower_bound;
    auto change_v = [&](uint64_t t, const range<t_k3_treap>& newv) {
        assert(vt[t].get_max() != inf);
        assert(newv.get_max() != inf);
        U -= vt[t].get_max() + 1;
        vt[t] = newv; // TODO no copy
        U += vt[t].get_max() + 1;
    };
    auto change_d = [&](uint64_t newd) {
        assert(d != inf);
        // TODO(niklas) why is this assert wrong?
        //assert(newd > d)
        d = std::max(d + 1, newd);
        if (d == inf)
            return;
        for (size_t t = 0; t < ranges.size(); ++t) {
            auto& stack2 = stacks[t];
            auto v = vt[t]; // TODO no copy
            while (stack2.size() && d > stack2.back().split_z()) {
                v = stack2.back();
                stack2.pop_back();
            }
            change_v(t, v);
        }
    };

    while (d < inf) {
        //std::cerr << "d=" << d << " L=" << L << " U=" << U << std::endl;
        //std::cerr << "total range queries = " << total << std::endl;

        //std::cerr << "vt = " << std::endl;
        //for (size_t i = 0 ; i < vt.size(); ++i)
            //std::cerr << "   " << vt[i].m_lo[2] << " - " << vt[i].m_hi[2]
                //<< " / " << vt[i].get_max() << std::endl;

        while (U < L) {
            uint64_t newd = inf;
            for (size_t t = 0; t < ranges.size(); ++t) {
                if (stacks[t].size())
                    newd = std::min(newd, stacks[t].back().split_z());
            }
            if (newd == inf)
                break;
            change_d(newd);
        }
        size_t t_shortest = -1;
        //std::cerr << "  max_z: ";
        for (size_t t = 0; t < ranges.size(); ++t) {
            assert(!vt[t].empty());
            //if (vt[t].is_leaf()) std::cerr << t << ":" << vt[t].get_max_z() << " ";
            if (!vt[t].is_leaf() || vt[t].get_max_z() != d) {
                // TODO better selection strategy?
                if (t_shortest == -1 || vt[t].get_max() < vt[t_shortest].get_max())
                    t_shortest = t;
            }
        }
        //std::cerr << std::endl << "  t_shortest = " << t_shortest << std::endl;
        if (t_shortest == -1) {
            //std::cerr << "REPORT " << d << " " << U << std::endl;
            result.insert(d, U);
            L = std::max(weight_lower_bound, get_cur_min() + 1);
            change_d(d + 1);
        } else {
            size_t t = t_shortest;
            auto left_right = vt[t].split();
            auto left = left_right.first;
            auto right = left_right.second;
            //std::cerr << "  changing " << t_shortest << " split_z=" << vt[t].split_z() << endl;
            if (d <= vt[t].split_z()) {
                if (!left.empty()) {
                    stacks[t].push_back(vt[t]);
                    assert(!stacks[t].back().empty());
                    change_v(t, left);
                } else {
                    // TODO(niklas) maybe not +1 here
                    change_d(vt[t].split_z());
                }
            } else {
                if (!right.empty()) {
                    assert(!right.empty());
                    change_v(t, right);
                } else {
                    if (stacks[t].empty())
                        break;
                    change_v(t, stacks[t].back());
                    if (stacks[t].size())
                        stacks[t].pop_back();
                    // TODO(niklas) maybe not +1 here
                    change_d(vt[t].split_z());
                }
            }
        }
    }

    return result.sorted_result();
}

}  // namespace k3_treap_algos
}  // namespace surf
