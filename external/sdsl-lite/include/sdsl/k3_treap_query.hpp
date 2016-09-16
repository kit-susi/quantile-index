//
// Created by Roberto Konow on 10/19/15.
//

#ifndef SURF_K3_TREAP_QUERY_HPP
#define SURF_K3_TREAP_QUERY_HPP
#include <iostream>
#include <sdsl/k3_treap.hpp>
#include <sdsl/k3_treap_algorithm.hpp>
#include <sdsl/k3_treap_helper.hpp>

#include <vector>
#include <queue>
#include <stack>

namespace surf {
namespace k3_treap_intersection {

typedef std::pair<size_t,float>  result_type;


template<typename k3_treap_type>
struct state {
    struct comp {
        bool operator()(result_type &a, result_type &b) const {
            return a.second > b.second;
        }
    };
    typedef std::priority_queue <result_type, std::vector<result_type>, comp>  queue_type;
    using stack_type = std::stack<sdsl::k3_treap_ns::top_k_iterator<k3_treap_type> >;
    using query_vector_type = std::vector<sdsl::k3_treap_ns::top_k_iterator<k3_treap_type>>;

    float                       upper;
    float                       lower;
    size_t                      document;
    size_t                      n_docs;
    std::vector<stack_type>     st;
    query_vector_type           query_it;
    std::vector<result_type>    result;
    queue_type                  current_result;

    state() {
        upper = 0.0;
        lower = 0.0;
        document = 0;
    }

    void dump_result() {
        result.reserve(current_result.size());
        while(!current_result.empty()) {
            auto data = current_result.top();
            result.push_back({data.first, data.second});
            current_result.pop();
        }
    }
};

template <typename k3_treap_type>
float get_score(sdsl::k3_treap_ns::top_k_iterator<k3_treap_type> &it) {
    auto result = *it;
    return std::get<1>(result);
}

template <typename k3_treap_type>
size_t get_document(sdsl::k3_treap_ns::top_k_iterator<k3_treap_type>  &it) {
    auto result = *it;
    auto point = std::get<0>(result);
    return std::get<2>(point);
}

template<typename k3_treap_type>
void changev(state<k3_treap_type>& t, size_t next, sdsl::k3_treap_ns::top_k_iterator<k3_treap_type>& it) {
    t.upper -= get_score(t.query_it[next]);
    t.query_it[next] = it;
    t.upper += get_score(t.query_it[next]);
}

template<typename k3_treap_type>
void changed(state<k3_treap_type>& t,size_t d) {
    t.document = d;
    for (size_t i = 0; i < t.query_it.size(); i++) {
        auto it = t.query_it[i];
        while (t.st[i].size() > 1 and t.document > get_document(t.st[i].top())) {
            it = t.st[i].top();
            t.st[i].pop();
        }
        changev(t,i, it);
    }
}

template<typename k3_treap_type>
void report(state<k3_treap_type>& t, size_t k) {
    t.current_result.emplace(t.document,t.upper);
    if (t.current_result.size() > k) {
        t.current_result.pop();
        t.lower = t.current_result.top().second;
    }
}

template<typename k3_treap_type>
std::vector<result_type> k3_treap_intersection(std::vector<sdsl::k3_treap_ns::top_k_iterator<k3_treap_type>> &q,
                                          size_t k,
                                          size_t n_docs=std::numeric_limits<size_t>::max())
{
    typedef sdsl::k3_treap_ns::top_k_iterator<k3_treap_type> k3_iterator;
    state<k3_treap_type> t;
    t.query_it = q; // <--- check if this is copying or not..
    t.st.resize(q.size());
    size_t next_possible = 0;
    bool search = false;
    for (size_t i = 0; i < t.query_it.size(); i++) {
        t.upper += get_score(t.query_it[i]);
        auto it = k3_iterator();
        t.st[i].push(it); // adding dummy iterator which is not valid;
        t.st[i].push(t.query_it[i]);
    }

    while (t.document < n_docs) {
        while (t.upper < t.lower) {
            size_t min_doc = std::numeric_limits<size_t>::max();
            for (size_t i = 0; i < q.size(); i++) {
                if (t.st[i].size() == 0) {
                    t.dump_result();
                    return t.result;
                }
                auto it = t.st[i].top();
                min_doc = std::min(get_document(it), min_doc);
            }
            if (min_doc == std::numeric_limits<size_t>::max()) {
                t.dump_result();
                return t.result;
            }
            changed(t,min_doc);
        }

        for (size_t i = 0; i < q.size(); i++) {
            if (!t.query_it[i]) {
                t.dump_result();
                return t.result;
            }
            if (get_document(t.query_it[i]) != t.document) {
                next_possible = i;
                search = true;
                break;
            }
        }
        if (!search) {
            report(t, k);
            changed(t,t.document + 1);
        } else {
            search = false;
            auto split = t.query_it[next_possible].split();
            if (t.document < get_document(t.query_it[next_possible])) {
                auto left = split[0];
                auto it = std::get<0>(*left);
                if (left) {
                    t.st[next_possible].push(t.query_it[next_possible]);
                    changev(t,next_possible, left);
                } else {
                    changed(t,get_document(t.query_it[next_possible]));
                }
            } else {
                auto right = split[1];
                auto it = std::get<0>(*right);
                if (right) {
                    changev(t,next_possible, right);
                } else {
                    auto it2 = std::get<0>(*t.st[next_possible].top());
                    changev(t,next_possible, t.st[next_possible].top());
                    t.st[next_possible].pop();
                    changed(t,get_document(t.query_it[next_possible]));
                }
            }
        }
    }
    t.dump_result();
    return t.result;
}
};
};

#endif //SURF_K3_TREAP_QUERY_HPP
