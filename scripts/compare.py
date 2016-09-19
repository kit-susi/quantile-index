#!/usr/bin/env python2
from subprocess import check_output
from prettytable import PrettyTable
import argparse
import os
import random
import re
import subprocess
import tempfile
import threading
import time
import traceback
import sys

def exe(cmd):
    try:
        return check_output(cmd)
    except Exception, e:
        print 'Error while running `%s`: %s' % (' '.join(cmd), e)
        raise

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

def get_collection_type(directory):
    if os.path.exists('%s/text_int_SURF.sdsl' % directory):
        return 'int'
    elif os.path.exists('%s/text_SURF.sdsl' % directory):
        return 'text'
    else:
        raise Exception("Not a susie collection: %d" % directory)

def print_side_by_side(a, b):
    for i in range(max(len(a), len(b))):
        d1, s1 = a[i] if i < len(a) else ('','')
        d2, s2 = b[i] if i < len(b) else ('','')
        print '%8s%10s   |%8s%10s' % (d1, s1, d2, s2)

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


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('targets', metavar='CONFIG', nargs='+',
            help='Configs to test')
    p.add_argument('--clear', default=False, action='store_true',
            help='Delete indexes before testing')
    p.add_argument('-c', dest='collection', required=True, metavar='DIRECTORY',
            help='Input collection')
    p.add_argument('-n', default=3, type=int, metavar='INT',
            help='ngram size for queries')
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

        print 'Starting new round'
        last_result = None
        with tempfile.NamedTemporaryFile() as f:
            if args.query_file:
                print 'Writing queries to %s' % args.query_file
                with open(args.query_file, 'w') as g:
                    g.write(queries)
            f.write(queries)
            f.flush()
            res = PrettyTable(['IDX','Avg','Median','Max','Sigma']);
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
                    sys.exit(1)

                if args.q == 1:
                    # check correctness only for q = 1
                    lines = out.strip().splitlines()
                    parts = [l.split(';') for l in lines]
                    result = [(int(docid), float(score)) for _, _, docid, score in parts]
                    if not result_makes_sense(result, args.e):
                        print 'Queries:', repr(queries)
                        print 'Result is not sorted:', result
                        assert 0
                    if (last_result is not None
                            and not results_are_same(last_result, result,
                                args.e, args.ignore_singletons)):
                        print 'Queries:', repr(queries)
                        print 'Different results:'
                        print_side_by_side(last_result, result)
                        sys.exit(1)
                    last_result = result
                else:
                    # otherwise print time
                    avg = float(out.split('time_per_query_avg = ')[1].split()[0])
                    median = int(out.split('time_per_query_median = ')[1].split()[0])
                    maxi = int(out.split('time_per_query_max = ')[1].split()[0])
                    sigma = float(out.split('time_per_query_sigma = ')[1].split()[0])
                    res.add_row([config,avg,median,maxi,sigma]);
                    #print '      Time per query (avg / median / max / sigma): %d %.2f %d %.2f' % (
                            #avg, median, maxi, sigma)
            if args.q != 1:
            	print 'Time per query (sorted by Avg)'
            	print res.get_string(sortby='Avg')
        seed = random.randrange(1000000)
