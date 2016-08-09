# Experiment 0 and 1
sh run_memory_usage.sh >> result_mem_usage.txt
# Experiment 2
cd benchmark
rm results/*
make experiment2
mv results/experiment2.txt experiment2a.txt
rm results/*
make experiment2
mv results/experiment2.txt experiment2b.txt

