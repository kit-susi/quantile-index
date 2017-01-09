#!/usr/bin/env python2
import os
import argparse

sampling = [4, 8, 16, 32, 64]
quantiles = [8, 16, 32, 64, 128]

def config(s,q):
    return 'IDX_NN_QUANTILE_LG_%d_%d' % (s, q)
def build_executable(s, q):
    os.system('./scripts/build_config.sh -d %s' % config(s, q))

def build_index(s, q, c, build_dir):
    print 'Building index with sampling %d and quantile %d' % (s,q)
    os.system('./%s/surf_index-%s -c %s' % (build_dir, config(s, q), c))

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('-c', dest='collection', required=True, metavar='DIRECTORY',
            help='Input collection')
    p.add_argument('-b', dest='build_dir', default='build/debug',
            help='Build directory (default: build/debug)')
    p.add_argument('--rebuild', default=False, action='store_true',
            help='Rebuilds all executables')
    args = p.parse_args()

    if args.rebuild:
        for s in sampling:
            for q in quantiles:
                build_executable(s, q)
    print "Done building all executables"
    for s in sampling:
        for q in quantiles:
            build_index(s, q, args.collection, args.build_dir)
