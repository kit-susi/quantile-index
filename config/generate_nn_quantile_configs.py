
def write_config_lg(s, q):
    file_name = 'IDX_NN_QUANTILE_LG_%d_%d' % (s, q)
    f = open(file_name+".config", "w")
    f.write('NAME=' + file_name + '\n')
    f.write('CSA_TYPE=sdsl::csa_wt<sdsl::wt_huff<sdsl::hyb_vector<>>,%d,%d, sdsl::text_order_sa_sampling<>, sdsl::text_order_isa_sampling_support<>>\n' % (s,s))
    f.write('DF_TYPE=surf::df_sada<CSA_TYPE,sdsl::rrr_vector<63>,sdsl::rrr_vector<63>::select_1_type,true,true>\n')
    f.write('WTD_TYPE=sdsl::wt_int<sdsl::bit_vector, sdsl::rank_support_v5<>, sdsl::select_support_scan<1>, sdsl::select_support_scan<0>>\n')
    f.write('KTWOTREAP_TYPE=sdsl::k2_treap<2,sdsl::rrr_vector<63>>\n')
    f.write('INDEX_TYPE=surf::idx_nn_quantile<CSA_TYPE, KTWOTREAP_TYPE, %d, 0, sdsl::sd_vector<>, sdsl::sd_vector<>::rank_1_type, sdsl::sd_vector<>::select_1_type, sdsl::rrr_vector<63>, sdsl::rrr_vector<63>::select_0_type, sdsl::rrr_vector<63>::select_1_type, true>\n' % q)

def write_config(s, q):
    file_name = 'IDX_NN_QUANTILE_%d_%d' % (s, q)
    f = open(file_name+".config", "w")
    f.write('NAME=' + file_name + '\n')
    f.write('CSA_TYPE=sdsl::csa_wt<sdsl::wt_huff<sdsl::hyb_vector<>>,%d,%d, sdsl::text_order_sa_sampling<>, sdsl::text_order_isa_sampling_support<>>\n' % (s,s))
    f.write('DF_TYPE=surf::df_sada<CSA_TYPE,sdsl::rrr_vector<63>,sdsl::rrr_vector<63>::select_1_type,true,true>\n')
    f.write('WTD_TYPE=sdsl::wt_int<sdsl::bit_vector, sdsl::rank_support_v5<>, sdsl::select_support_scan<1>, sdsl::select_support_scan<0>>\n')
    f.write('KTWOTREAP_TYPE=sdsl::k2_treap<2,sdsl::rrr_vector<63>>\n')
    f.write('INDEX_TYPE=surf::idx_nn_quantile<CSA_TYPE, KTWOTREAP_TYPE, %d, 0>\n' % q)

for s in [4, 8, 16, 32, 64]:
    for q in [8, 16, 32, 64, 128]:
        write_config(s,q)
        write_config_lg(s,q)
