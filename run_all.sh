# Experiment 0 and 1
sh run_memory_usage.sh >> result_mem_usage.txt
# Experiment 2
cd benchmark
# A
rm results/*
echo "10;5;1;0;0" >> experiment2.config
make experiment2
mv results/experiment2.txt experiment2a.txt
# B
rm results/*
echo "10;5;1;140;30" >> experiment2.config
make experiment2
mv results/experiment2.txt experiment2b.txt

