//
// Created by Roberto Konow on 10/19/15.
//
#include <iostream>
#include "sdsl/int_vector_buffer.hpp"
#include "sdsl/k3_treap.hpp"
#include "sdsl/k3_treap_algorithm.hpp"
#include "sdsl/k3_treap_helper.hpp"
#include "surf/k3_treap_query.hpp"
#include "sdsl/rrr_vector.hpp"
using namespace sdsl;
using namespace std;

typedef k3_treap<2,rrr_vector<63>>                      k3_treap_type;
typedef array<k3_treap_ns::point_type,2>                query_element;
typedef vector<query_element>                           query_type;
typedef k3_treap_ns::top_k_iterator<k3_treap_type>      k3_iterator;
typedef vector<k3_iterator>                             query_vector_type;

int main(int argc, char* argv[])
{
    typedef k3_treap<2, rrr_vector<63>> k3_rrr;
    k3_rrr k3t;
    construct_im(k3t,
                 {{1, 0, 1, 3}, {3, 2, 4, 100}, {3, 2, 3, 1001}, {3, 1, 1, 1}, {5, 5, 5, 5}, {4, 4, 5, 1},
                  {6, 9, 5, 10}, {7, 8, 3, 11}, {8, 8, 1, 1003}, {1, 2, 3, 4}});

    k3_treap_ns::point_type p1 = {1, 0, 0};
    k3_treap_ns::point_type p2 = {5, 5, 10};
    k3_treap_ns::point_type pp1 = {6, 6, 0};
    k3_treap_ns::point_type pp2 = {10, 10, 10};
    query_element query1 = {p1, p2};
    query_element query2 = {pp1, pp2};
    query_type query = {query1, query2};

    query_vector_type qv;
    auto it_q1 = top_k(k3t, query1[0], query1[1]);
    auto it_q2 = top_k(k3t, query2[0], query2[1]);

    qv.push_back(it_q1);
    qv.push_back(it_q2);
    vector<pair<size_t, float> > result;
    result = surf::k3_treap_intersection::k3_treap_intersection(qv, 10);
    for (auto x : result) {
        cout << x.first << " , " << x.second << endl;
    }
}
