#!/usr/bin/env python2
import sys, os
sys.path.append(os.path.dirname(__file__) + "/../lib/py")

from prettytable import PrettyTable
import argparse
import random
import re
import subprocess
import tempfile
import threading
import time
import traceback
from gen_queries import exe, make_query_file, get_collection_type

def is_sorted(lst, eps):
    """ Checks of lst of floats is sorted. """
    for a, b in zip(lst, lst[1:]):
        if b < a - eps:
            return False
    return True

def result_makes_sense(result, eps):
    """ result is list of (docid, score), where docid is int and score is float. """
    return is_sorted([score for _, score in reversed(result)], eps)

def delete_singletons(lst, eps):
    return [(d, w) for d, w in lst if abs(w - 1) > eps]

def results_are_same(a, b, eps, ignore_singletons=False):
    """ a, b are lists of (docid, score), where docid is int and score is float. """
    if ignore_singletons:
        a = delete_singletons(a, eps)
        b = delete_singletons(b, eps)
    if len(a) != len(b):
        return False
    # scores are the same
    for (d1, s1), (d2, s2) in zip(a, b):
        if abs(s1 - s2) > eps:
            return False
    scores1 = dict(a)
    scores2 = dict(b)
    # doc -> score association is the same for common docs
    for doc in set(scores1.keys()) & set(scores2.keys()):
        if abs(scores1[doc] - scores2[doc]) > eps:
            return False
    return True

def print_side_by_side(a, b):
    for i in range(max(len(a), len(b))):
        d1, s1 = a[i] if i < len(a) else ('','')
        d2, s2 = b[i] if i < len(b) else ('','')
        print '%8s%10s   |%8s%10s' % (d1, s1, d2, s2)

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('targets', metavar='CONFIG', nargs='+',
            help='Configs to test')
    p.add_argument('--clear', default=False, action='store_true',
            help='Delete indexes before testing')
    p.add_argument('-c', dest='collection', required=True, metavar='DIRECTORY',
            help='Input collection')
    p.add_argument('-n', default='3', metavar='INT_OR_RANGE',
            help='ngram size for queries (default: 3). Can also be a range: -n 3-10')
    p.add_argument('-r', default=20, type=int, metavar='INT',
            help='Number of testing rounds')
    p.add_argument('-q', default=1, type=int, metavar='INT',
            help='Number of queries per round. For q > 1, no verification is done')
    p.add_argument('-k', default=20, type=int, metavar='INT',
            help='Retrieve top k documents')
    p.add_argument('-e', default=1e-6, type=float, metavar='FLOAT',
            help='Epsilon for score comparisons')
    p.add_argument('-o', dest='min_sampling', default=0, type=int, metavar='INT',
            help='Only use ngrams that occur at least once every X samples. Useful for intersection queries')
    p.add_argument('-i', dest='intersection', default=1, type=int, metavar='INT',
            help='Generate intersection queries with the given number of terms')
    p.add_argument('-m', dest='multi_occ', default=True, action='store_true',
            help='Find only multi-occurences')
    p.add_argument('--no_multi_occ', dest='multi_occ', action='store_false',
            help='Also find singleton occs')
    p.add_argument('--seed', default=random.randrange(1000000), type=int, metavar='INT',
            help='Random seed')
    p.add_argument('--ignore_singletons', default=False, action='store_true',
            help='Ignore single-document occurences (weight 1)')
    p.add_argument('--query_file',
            help='Store queries in the given file')
    p.add_argument('-b', dest='build_dir', default='build/release',
            help='Build directory (default: build/release)')

    # NOTE currently parallel construction does not work, so this can not be turned off
    p.add_argument('--sequential', default=True, action='store_true',
            help='Build collections in sequence, instead of parallel')

    args = p.parse_args()
    collection_type = get_collection_type(args.collection)

    if args.clear:
        print 'Deleting old indexes'
        subprocess.check_call(['rm', '-r', args.collection + '/index'])

    print 'Using build dir = %s' % args.build_dir

    threads = []
    for config in args.targets:
        def build(config):
            print 'Building index for config %s' % config
            cmd = [
                '%s/surf_index-%s' % (args.build_dir, config),
                '-c', args.collection
            ]
            print '    Running command: %s' % ' '.join(cmd)
            try:
                exe(cmd)
            except Exception:
                sys.exit(1)

        t = threading.Thread(target=build, args=(config,))
        t.start()
        if args.sequential:
            t.join()
        else:
            threads.append(t)

    print 'Waiting for %d builder threads' % len(threads)
    for t in threads:
        t.join()

    if args.q != 1:
        print 'WARNING: No correctness will be checked. Use -q 1 if you want to do that'

    seed = args.seed
    for _ in range(args.r):
        print 'Seed = %d' % seed
        random.seed(seed)

        print 'Generating queries...'
        queries = make_query_file(args, seed)

        print 'Starting new round'
        last_result = None
        with tempfile.NamedTemporaryFile() as f:
            if args.query_file:
                print 'Writing queries to %s' % args.query_file
                with open(args.query_file, 'w') as g:
                    g.write(queries)
            f.write(queries)
            f.flush()

            cols = ['IDX','Avg','Median','Max','Sigma']
            res = PrettyTable(cols)
            res.align[cols[0]] = 'l'
            for c in cols[1:]: res.align[c] = 'r'

            result_count = 0
            for config in args.targets:
                print '    Running with %s' % config
                cmd = [
                    '%s/surf_query-%s' % (args.build_dir, config),
                    '-c', args.collection,
                    '-q', f.name,
                    '-k', str(args.k),
                ]
                if args.q == 1:
                    cmd += ['-v']
                if args.multi_occ:
                    cmd += ['-m']
                if args.intersection > 1:
                    cmd += ['-i']

                try:
                    out = exe(cmd)
                except Exception, e:
                    print 'FAIL: Query program crashed'
                    sys.exit(1)

                if args.q == 1:
                    # check correctness only for q = 1
                    lines = out.strip().splitlines()
                    parts = [l.split(';') for l in lines]
                    result = [(int(docid), float(score))
                                for _, _, docid, score in parts]
                    result_count = max(result_count, len(result))

                    if not result_makes_sense(result, args.e):
                        print 'Queries:', repr(queries)
                        print 'FAIL with seed %d: Result is not sorted:' % seed, result
                        assert 0
                    if (last_result is not None
                            and not results_are_same(last_result, result,
                                args.e, args.ignore_singletons)):
                        print 'Queries:', repr(queries)
                        print 'FAIL with seed %d: Different results:' % seed
                        print_side_by_side(last_result, result)
                        sys.exit(1)
                    last_result = result
                else:
                    # otherwise print time
                    avg = float(out.split('time_per_query_avg = ')[1].split()[0])
                    median = int(out.split('time_per_query_median = ')[1].split()[0])
                    maxi = int(out.split('time_per_query_max = ')[1].split()[0])
                    sigma = float(out.split('time_per_query_sigma = ')[1].split()[0])
                    res.add_row([config,avg,median,maxi,sigma])
            print 'Results: %d' % result_count
            if args.q != 1:
                print 'Time per query (sorted by Avg)'
                print res.get_string(sortby='Avg')
        seed = random.randrange(1000000)
