#!/usr/bin/env python2
import os
import argparse
from multiprocessing import Process
import subprocess

sampling = [4, 8, 16, 32, 64]
quantiles = [8, 16, 32, 64, 128]

output_file = open('result_file', "w")
def output(s):
    output_file.write(s + "\n")
    output_file.flush()

def quantile_config(s, q):
    return 'IDX_NN_QUANTILE_%d_%d' % (s, q)
def quantile_lg_config(s, q):
    return 'IDX_NN_QUANTILE_LG_%d_%d' % (s, q)
def nn_config(s):
    return 'IDX_NN_%d' % s
def nn_lg_config(s):
    return 'IDX_NN_LG_%d' % s

def build_executable(configs):
    os.system('./scripts/build_config.sh -d %s' % ' '.join(configs))
    os.system('./scripts/build_config.sh %s' % ' '.join(configs))

def build_index(build_dir, config, c):
    print 'Building index with sampling %d and quantile %d' % (s,q)
    os.system('./%s/surf_index-%s -c %s' % (build_dir, config, c))

def run_sequential(processes):
    for p in processes:
        p.start()
        p.join()

def run_parallel(processes):
    for p in processes:
        p.start()
    for p in processes:
        p.join()

def lastLine(s):
    return s.split('\n')[-2]

def printIndexSizes(configs, collection, build_dir):
    for config in configs:
    	proc = subprocess.Popen(['./%s/surf_index-%s' % (build_dir, config), "-c", collection, "-m m"], stdout=subprocess.PIPE)
    	out = proc.communicate()[0]
    	output(config + ", " + lastLine(out.upper()))

def printIndexSpeed(configs, collection, build_dir):
    command = ['./scripts/compare.py'] + configs + ["-c", collection, "-b", build_dir, '--no_multi_occ', '-q 10000', '-k 10', '-n 4', '-r 1']
    print command
    proc = subprocess.Popen(command, stdout=subprocess.PIPE)
    output(proc.communicate()[0])

def check_results(configs, collection, build_dir):
    print "Checking results"
    command = ['./scripts/compare.py'] + configs + ['-c', collection, '-b', build_dir, '--no_multi_occ', '-q 1', '-k 10', '-n 4', '-r 10']
    proc = subprocess.Popen(command, stdout=subprocess.PIPE)
    print proc.communicate()[0]

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('-c', dest='collection', required=True, metavar='DIRECTORY',
            help='Input collection')
    p.add_argument('-b', dest='build_dir', default='build/release',
            help='Build directory (default: build/release)')
    p.add_argument('--rebuild', default=False, action='store_true',
            help='Rebuilds all executables')
    p.add_argument('--check', default=False, action='store_true',
            help='Only check if results all match.')
    args = p.parse_args()
 
    quantile_configs = [quantile_config(s,q) for s in sampling for q in quantiles]
    quantile_lg_configs = [quantile_lg_config(s,q) for s in sampling for q in quantiles]
    nn_configs = [nn_config(s) for s in sampling]
    nn_lg_configs = [nn_lg_config(s) for s in sampling]
    configs = quantile_configs + quantile_lg_configs + nn_configs + nn_lg_configs
    if args.rebuild:
	build_executable(configs)

    if args.check:
        check_results(configs, args.collection, args.build_dir)
    else:
        print "Done building all executables"
        run_sequential([Process(target=build_index, args=(args.build_dir, config, args.collection)) for config in configs])

        print "Index sizes" 
        printIndexSizes(configs , args.collection, args.build_dir)

        print "Comparing indizes"
        printIndexSpeed(configs, args.collection, args.build_dir)
