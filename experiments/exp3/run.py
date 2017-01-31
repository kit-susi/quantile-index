#!/usr/bin/env python2
import os
import argparse
import subprocess

q = 64
s = 16
collections = ['ENWIKIBIG', 'ENWIKISML', 'SOURCES', 'REVISIONS']

def quantile_config(s, q):
    return ('QUANTILE', 'IDX_NN_QUANTILE_%d_%d' % (s, q))
def quantile_lg_config(s, q):
    return ('QUANTILELG', 'IDX_NN_QUANTILE_LG_%d_%d' % (s, q))
def nn_config(s):
    return ('NN', 'IDX_NN_%d' % s)
def nn_lg_config(s):
    return ('NNLG', 'IDX_NN_LG_%d' % s)

def lastLine(s):
    return s.split('\n')[-2]

def printconf(config):
    return config.replace('_', ';')

configs = [quantile_config(s, q), quantile_lg_config(s, q), nn_config(s), nn_lg_config(s)]
def printIndexSizes(base, build_dir, output):
    output_file = open(output, "w")
    output_file.write('index; collection; s; q; csa; doc; bv; grid; G_q\n')
    for collection in collections:
        col = base + '/' + collection
        for (desc, config) in configs:
            proc = subprocess.Popen(['./%s/surf_index-%s' % (build_dir, config), "-c", col, "-m m"],
                    stdout=subprocess.PIPE)
            out = proc.communicate()[0]
            output_file.write(desc + "; " + collection + "; " + str(s) + "; " + str(q) + "; "
                    + lastLine(out.upper()) + "\n")

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('-c', dest='collection', required=True, metavar='DIRECTORY',
            help='Input collection base dir')
    p.add_argument('-o', dest='output', required=True, metavar='DIRECTORY',
            help='Output csv file')
    p.add_argument('-b', dest='build_dir', default='build/release',
            help='Build directory (default: build/release)')
    args = p.parse_args()
    print "Index sizes" 
    printIndexSizes(args.collection, args.build_dir, args.output)
