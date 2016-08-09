for COL in ENWIKISML ENWIKIBIG REVISIONS REDDIT
do
	echo $COL
	./build/surf_index-IDX_NN -c collections/$COL -m m | tail -2 
	for IDX in O_SD O_HYB O
	do
		./build/surf_index-IDX_$IDX -c collections/$COL -m m | tail -1 
	done
	echo $COL"INT"
	./build/surf_index-IDX_NN_INT -c collections/$COL"INT" -m m | tail -2 
	for IDX in O_SD_INT O_HYB_INT O_INT
	do
		./build/surf_index-IDX_$IDX -c collections/$COL"INT" -m m | tail -1 
	done
done
