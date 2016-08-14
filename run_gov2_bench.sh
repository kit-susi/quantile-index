./build/surf_query-IDX_NN_INT -c /data01/labeit/collections/gov2 -q benchmark/pattern/GOV2INT.1.int.txt -k 10 -s 0 > gov_2_bench
./build/surf_query-IDX_O_INT -c /data01/labeit/collections/gov2 -q benchmark/pattern/GOV2INT.1.int.txt -k 10 -s 0 >> gov_2_bench
./build/surf_query-IDX_NN_INT -c /data01/labeit/collections/gov2 -q benchmark/pattern/GOV2INT.1.int.txt -k 10 -s 30 >> gov_2_bench
./build/surf_query-IDX_O_INT -c /data01/labeit/collections/gov2 -q benchmark/pattern/GOV2INT.1.int.txt -k 10 -s 30 >> gov_2_bench
