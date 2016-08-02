./build/surf_index-IDX_NN -c /fast-data1/labeit/collections/ENWIKIBIG -m m > mem_big_nn &
./build/surf_index-IDX_O -c /fast-data1/labeit/collections/ENWIKIBIG -m m > mem_big_o &

./build/surf_index-IDX_NN -c /fast-data1/labeit/collections/ENWIKISML -m m > mem_sml_nn &
./build/surf_index-IDX_O -c /fast-data1/labeit/collections/ENWIKISML -m m > mem_sml_o &

./build/surf_index-IDX_NN_INT -c /fast-data1/labeit/collections/ENWIKIBIGINT -m m > mem_big_nn_int &
./build/surf_index-IDX_O_INT -c /fast-data1/labeit/collections/ENWIKIBIGINT -m m > mem_big_o_int &

./build/surf_index-IDX_NN_INT -c /fast-data1/labeit/collections/ENWIKISMLINT -m m > mem_sml_nn_int &
./build/surf_index-IDX_O_INT -c /fast-data1/labeit/collections/ENWIKISMLINT -m m > mem_sml_o_int &
wait
echo "Done"
