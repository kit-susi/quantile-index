#!/usr/bin/env python2
import argparse
import random
import subprocess
import tempfile

def get_ngram(text, n):
    while True:
        i = random.randrange(len(text) - n)
        ngram = text[i:i+n]
        if '\1' not in ngram:
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

def results_are_same(a, b, eps):
    """ a, b are lists of (docid, score), where docid is int and score is float. """
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


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('targets', metavar='CONFIG', nargs='+',
            help='Configs to test')
    p.add_argument('--clear', default=False, action='store_true',
            help='Delete indexes before testing')
    p.add_argument('-c', metavar='DIRECTORY', help='Input collection')
    p.add_argument('-n', default=3, type=int, metavar='INT',
            help='ngram size for queries')
    p.add_argument('-s', default=20, type=int, metavar='INT',
            help='Number of sample queries')
    p.add_argument('-k', default=20, type=int, metavar='INT',
            help='Retrieve top k documents')
    p.add_argument('-e', default=1e-6, type=float, metavar='FLOAT',
            help='Epsilon for score comparisons')

    args = p.parse_args()
    with open(args.c + '/text_SURF.sdsl') as f:
        # strip 8-byte header
        text = f.read()[8:]

    if args.clear:
        print 'Deleting old indexes'
        subprocess.check_call(['rm', '-r', args.c + '/index'])

    for config in args.targets:
        print 'Building index for config %s' % config
        cmd = [
            'build/surf_index-%s' % config,
            '-c', args.c
        ]
        print '    Running command: %s' % ' '.join(cmd)
        subprocess.check_output(cmd)

    queries = [get_ngram(text, args.n) for _ in range(args.s)]
    for q in queries:
        print 'Query:', repr(q)
        last_result = None
        with tempfile.NamedTemporaryFile() as f:
            f.write(q + '\n')
            f.flush()
            for config in args.targets:
                print '    Running with %s' % config
                out = subprocess.check_output([
                    'build/surf_query-%s' % config,
                    '-c', args.c,
                    '-q', f.name,
                    '-v'
                ])
                lines = out.strip().splitlines()
                parts = [l.split(';') for l in lines]
                result = [(int(docid), float(score)) for _, _, docid, score in parts]
                if not result_makes_sense(result, args.e):
                    print 'Result is not sorted:', result
                    assert 0
                if (last_result is not None
                        and not results_are_same(last_result, result, args.e)):
                    print last_result, '!=', result
                    assert 0
