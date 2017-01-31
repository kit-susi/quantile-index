#!/usr/bin/env python2
import os
import argparse
import subprocess

quantiles = [8, 16, 32, 64, 128]
sampling = [4, 8, 16, 32, 64]
collections = ['ENWIKIBIG', 'ENWIKISML', 'SOURCES', 'REVISIONS']
#collections = ['TEST_TXT']
n = 5

def quantile_config(s, q):
    return 'IDX_NN_QUANTILE_%d_%d' % (s, q)
def quantile_lg_config(s, q):
    return 'IDX_NN_QUANTILE_LG_%d_%d' % (s, q)
def nn_config(s):
    return 'IDX_NN_%d' % s
def nn_lg_config(s):
    return 'IDX_NN_LG_%d' % s

def printconf(config):
    return config.replace('_', ';')

def queries(collection):
    return '/tmp/queries_%s' % collection
def ex(out, key):
    res = filter(lambda l: key in l, out.split('\n'))
    if len(res) > 0:
        return res[0].split("=")[1]
    else:
        return '0.0'

def getTiming(config, collection, col):
    proc = subprocess.Popen(['./build/release/surf_query-%s' % config, '-c', col, '-q', queries(collection)],
            stdout=subprocess.PIPE)
    out = proc.communicate()[0]
    print out
    return '%s; %s; %s; %s; %s; %s' % (
            ex(out, 'time_per_query_avg'),
            ex(out, 'time_per_query_median'),
            ex(out, 'time_per_query_max'),
            ex(out, 'time_per_query_sigma'),
            ex(out, 'index_size'),
            ex(out, 'input_size'))

def generate_queries(base):
    for collection in collections:
        col = base + '/' + collection
        proc = subprocess.Popen(['./scripts/gen_queries.py', '-q', '50000', '-c', col, '-n', str(n), queries(collection)])
	proc.communicate()

def printSpaceTime(base, build_dir, output):
    output_file = open(output, "w")
    output_file.write('index; collection; s; q; avg; media; max; sigma; indexsize; inputsize\n')
    for collection in collections:
        col = base + '/' + collection
        for s in sampling:
            for q in quantiles:
                for (desc,config) in [("QUANTILE", quantile_config(s, q)), ("QUANTILELG", quantile_lg_config(s,q))]:
    	            output_file.write(desc + "; " + collection + "; " + str(s) + "; " + str(q) + "; " + getTiming(config, collection, col) + "\n")
            q = 1
            for (desc,config) in [("NN", nn_config(s)), ("NNLG", nn_lg_config(s))]:
                output_file.write(desc + "; " + collection + "; " + str(s) + "; " + str(q) + "; " + 
                        getTiming(config, collection, col) + "\n")


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('-c', dest='collection', required=True, metavar='DIRECTORY',
            help='Input collection base dir')
    p.add_argument('-o', dest='output', required=True, metavar='DIRECTORY',
            help='Output csv file')
    p.add_argument('-b', dest='build_dir', default='build/release',
            help='Build directory (default: build/release)')
    args = p.parse_args()
    generate_queries(args.collection)
    printSpaceTime(args.collection, args.build_dir, args.output)
