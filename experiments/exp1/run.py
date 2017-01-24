#!/usr/bin/env python2
import os
import argparse
import subprocess

quantiles = [1, 2, 4, 8, 16, 32, 64, 128, 256]
sampling = 4
collections = ['ENWIKISML', 'REVISION', 'SOURCES']

def quantile_config(s, q):
    return 'IDX_NN_QUANTILE_LG_%d_%d' % (s, q)

def lastLine(s):
    return s.split('\n')[-2]

def printconf(config):
    return config.replace('_', ';')

def printIndexSizes(base, build_dir, output):
    output_file = open(output, "w")
    output_file.write('index; collection; s; q; csa; doc; bv; grid; |G_q|\n')
    for collection in collections:
        for q in quantiles:
            col = base + '/' + collection
            config = quantile_config(sampling, q)
    	    proc = subprocess.Popen(['./%s/surf_index-%s' % (build_dir, config), "-c", col, "-m m"],
                    stdout=subprocess.PIPE)
    	    out = proc.communicate()[0]
    	    output_file.write(config + "; " + collection + "; " + str(sampling) + "; " + str(q) + "; "
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
