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

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('targets', metavar='CONFIG', nargs='+',
            help='Configs to test')
    p.add_argument('--clear', default=False, action='store_true',
            help='Delete indexes before testing')
    p.add_argument('-c', metavar='DIRECTORY', help='Input collection')
    p.add_argument('-n', default=3, type=int, metavar='N',
            help='ngram size for queries')
    p.add_argument('-s', default=20, type=int, metavar='N',
            help='Number of sample queries')

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
        last_weights = None
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
                weights = [int(l.split(';')[3]) for l in lines]
                if weights != last_weights and last_weights is not None:
                    print last_weights, '!=', weights
                    assert 0
