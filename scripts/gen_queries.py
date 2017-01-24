#!/usr/bin/env python2
import sys, os
sys.path.append(os.path.dirname(__file__) + "/../lib/py")

import argparse
import random
import re
import subprocess
import tempfile
import threading
import time
import traceback
from subprocess import check_output

def exe(cmd):
    try:
        return check_output(cmd)
    except Exception, e:
        print 'Error while running `%s`: %s' % (' '.join(cmd), e)
        raise

def gen_queries(n, args, seed):
    return exe(['%s/gen_patterns' % args.build_dir,
        '-c', args.collection,
        '-m', str(args.n),
        '-s', str(seed),
        '-x', str(n),
        ]
        + (['-o', str(args.min_sampling)] if args.min_sampling else [])
        + (['-i'] if get_collection_type(args.collection) == 'int' else [])
        ).rstrip('\r\n').splitlines()

def get_collection_type(directory):
    if os.path.exists('%s/text_int_SURF.sdsl' % directory):
        return 'int'
    elif os.path.exists('%s/text_SURF.sdsl' % directory):
        return 'text'
    else:
        raise Exception("Not a SUSI collection: %d" % directory)

def make_query_file(args, seed):
    if args.intersection == 1:
        # singleterm
        queries = '\n'.join(gen_queries(args.q, args, seed=seed)) + '\n'
    elif args.intersection > 1:
        # multi-term with intersection
        queries = ''
        for i in range(args.q):
            queries += '\1'.join(gen_queries(args.intersection,
                                             args,
                                             seed=seed + i)) + '\n'
    else:
        assert 0, "Invalid value for -i: %d" % args.intersection
    return queries

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('-c', dest='collection', required=True, metavar='DIRECTORY',
            help='Input collection')
    p.add_argument('-n', default=3, type=int, metavar='INT',
            help='ngram size for queries (default: 3)')
    p.add_argument('-q', default=1000, type=int, metavar='INT',
            help='Number of queries to create (default: 1000)')
    p.add_argument('-o', dest='min_sampling', default=0, type=int, metavar='INT',
            help='Only use ngrams that occur at least once every X samples. Useful for intersection queries')
    p.add_argument('-i', dest='intersection', default=1, type=int, metavar='INT',
            help='Generate intersection queries with the given number of terms')
    p.add_argument('--seed', default=random.randrange(1000000), type=int, metavar='INT',
            help='Random seed')
    p.add_argument('-b', dest='build_dir', default='build/release',
            help='Build directory (default: build/release)')

    p.add_argument('query_file', help='Store queries in the given file')

    args = p.parse_args()

    with open(args.query_file, 'w') as f:
        f.write(make_query_file(args, args.seed))
