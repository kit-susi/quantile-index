#!/usr/bin/env python2
import argparse
import random
import re
import subprocess
import tempfile
import threading
import time

def get_ngram(text, n):
    while True:
        i = random.randrange(len(text) - n)
        ngram = text[i:i+n]
        if '\1' not in ngram and '\n' not in ngram:
            return ngram

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
    p.add_argument('-c', required=True, metavar='DIRECTORY',
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
    p.add_argument('--seed', default=random.randrange(1000000), type=int, metavar='INT',
            help='Random seed')
    # NOTE currently parallel construction does not work
    p.add_argument('--sequential', default=True, action='store_true',
            help='Build collections in sequence, instead of parallel')
    p.add_argument('--ignore_singletons', default=False, action='store_true',
            help='Ignore single-document occurences (weight 1)')
    p.add_argument('--query_file',
            help='Store queries in the given file')

    args = p.parse_args()
    if args.clear:
        print 'Deleting old indexes'
        subprocess.check_call(['rm', '-r', args.c + '/index'])

    threads = []
    for config in args.targets:
        def build(config):
            print 'Building index for config %s' % config
            cmd = [
                'build/surf_index-%s' % config,
                '-c', args.c
            ]
            print '    Running command: %s' % ' '.join(cmd)
            subprocess.check_output(cmd)

        t = threading.Thread(target=build, args=(config,))
        t.start()
        if args.sequential:
            t.join()
        else:
            threads.append(t)

    print 'Waiting for %d builder threads' % len(threads)
    for t in threads:
        t.join()

    with open(args.c + '/text_SURF.sdsl') as f:
        # strip 8-byte header
        text = f.read()[8:]

    print 'Seed = %d' % args.seed
    random.seed(args.seed)

    if args.q != 1:
        print 'WARNING: No correctness will be checked. Use -q 1 if you want to do that'

    for _ in range(args.r):
        queries = []
        for _ in range(args.q):
            queries.append(get_ngram(text, args.n))
        print 'Starting new round'
        last_result = None
        with tempfile.NamedTemporaryFile() as f:
            if args.query_file:
                print 'Writing queries to %s' % args.query_file
                with open(args.query_file, 'w') as g:
                    g.write('\n'.join(queries) + '\n')
            f.write('\n'.join(queries) + '\n')
            f.flush()
            for config in args.targets:
                print '    Running with %s' % config
                cmd = [
                    'build/surf_query-%s' % config,
                    '-c', args.c,
                    '-q', f.name,
                    '-k', str(args.k),
                ]
                if args.q == 1:
                    cmd += ['-v']

                out = subprocess.check_output(cmd)
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
                        assert 0
                    last_result = result
                else:
                    # otherwise print time
                    avg = int(out.split('time_per_query_avg = ')[1].split()[0])
                    median = int(out.split('time_per_query_median = ')[1].split()[0])
                    maxi = int(out.split('time_per_query_max = ')[1].split()[0])
                    print '      Time per query (avg / median / max): %d %d %d' % (avg, median, maxi)
        seed = random.randrange(1000000)
        print "New seed = %d" % seed
        random.seed(seed)
