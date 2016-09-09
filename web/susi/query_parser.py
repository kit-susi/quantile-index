from collections import namedtuple
import re

from funcparserlib.lexer import make_tokenizer, Token, LexerError
from funcparserlib.parser import (Parser, some, a, maybe, many, finished, skip,
                                  forward_decl, pure, NoParseError)

ENCODING = u'ascii'
regexps = {
    u'escaped': ur'''
        \\
          ((?P<standard>["\\/'bfnrt])                    # Standard escapes
          |(x(?P<code>[0-9A-Fa-f]{2})))               # Character code escape
        ''',
    u'escaped_regex': ur'\\.',
    u'unescaped_str': ur'''
        [^"\\]                              # Unescaped double-quoted string: avoid ["\\]
        ''',
    u'unescaped_regex': ur'''
        [^/\\]                              # Unescaped regex: avoid [/\\]
        ''',
}
re_esc = re.compile(regexps[u'escaped'], re.VERBOSE)


def tokenize(s):
    """str -> Sequence(Token)"""
    specs = [
        (u'Space', (ur'[ \t\r\n]+',)),
        (u'String', (ur'"(%(unescaped_str)s | %(escaped)s)*"' % regexps, re.VERBOSE)),
        (u'Regex', (ur'/(%(unescaped_regex)s | %(escaped_regex)s)*/[i]*' % regexps, re.VERBOSE)),
        (u'Op', (ur'or|and|not|[\(\)]',)),
        (u'Prefix', (ur'client:|server:|any:',)),
    ]
    useless = [u'Space']
    t = make_tokenizer(specs)
    return [x for x in t(s) if x.type not in useless]

def unescape_str(s):
    std = {
        u'"': u'"', u"'": u"'", u"/": u"/", u'\\': u'\\', u'/': u'/',
        u'b': u'\b', u'f': u'\f', u'n': u'\n', u'r': u'\r', u't': u'\t',
    }

    def sub(m):
            if m.group(u'standard') is not None:
                return std[m.group(u'standard')]
            else:
                return unichr(int(m.group(u'code'), 16))

    return re_esc.sub(sub, s)

# TODO this seems like a good idea, but we are now injectable with arbitrary
# binary characters.
def unescape_regex(s):
    return s

String = namedtuple('String', ('value', 'context'))
Regex = namedtuple('Regex', ('value', 'modifiers', 'context'))
Or = namedtuple('Or', ('x', 'y'))
And = namedtuple('And', ('x', 'y'))
Not = namedtuple('Not', ('x',))
Any = namedtuple('Any', ())

def _parse(seq):
    const = lambda x: lambda _: x
    tokval = lambda x: x.value
    toktype = lambda t: some(lambda x: x.type == t) >> tokval
    op = lambda s: a(Token(u'Op', s)) >> tokval
    op_ = lambda s: skip(op(s))

    def make_string(args):
        context, value = args
        if not context: context = 'any:'
        return String(unescape_str(value[1:-1]), context[:-1])

    def make_regex(args):
        context, value = args
        value, modifiers = value.rsplit('/', 1)
        value = value[1:]
        if not context: context = 'any:'
        return Regex(unescape_regex(value), modifiers, context[:-1])

    def make_or(args):
        return Or(*args)

    def make_and(args):
        return And(*args)

    def make_not(x):
        return Not(x)

    context = maybe(toktype(u'Prefix'))
    string = (context + toktype(u'String')) >> make_string
    regex = (context + toktype(u'Regex')) >> make_regex

    par_term = forward_decl()
    simple_term = forward_decl()
    term = forward_decl()
    not_term = forward_decl()
    and_term = forward_decl()
    or_term = forward_decl()

    par_term.define(op_(u'(') + term + op_(u')'))

    simple_term.define(
        par_term
        | string
        | regex
    )
    not_term.define(
        op_('not') + not_term >> make_not
        | simple_term
    )
    and_term.define(
        not_term + op_('and') + and_term >> make_and
        | not_term
    )
    or_term.define(
        and_term + op_('or') + or_term >> make_or
        | and_term
    )

    term.define(or_term)

    eof = skip(toktype(u'EOF'))
    filter_expr = (term + eof) | (eof >> const(Any()))
    return filter_expr.parse(seq)

def parse_query(s):
    tokens = tokenize(s) + [Token(u'EOF', '')]
    return _parse(tokens)


import unittest
class Test(unittest.TestCase):
    def test_success(self):
        tests = [
            (
                r'not ("a" and client:"b" or server:"f" and client:/g/) or not not any:"c" and server:/d/',
                Or(Not(Or(And(String('a', 'any'), String('b', 'client')), And(String('f', 'server'),
                    Regex('g', '', 'client')))), And(Not(Not(String('c', 'any'))), Regex('d', '', 'server')))
            ),
            (
                r'''"abc\n\"\b/\/\\\'\x12\x41"''',
                String('abc\n"\b//\\\'\x12A', 'any')
            ),
            (
                '',
                Any()
            ),
            (
                r'''/abc\n\"\b\/\\\'\x41\w\s/''',
                Regex(r'''abc\n\"\b\/\\\'\x41\w\s''', '', 'any')
            ),
            (
                r'''/foobar/i''',
                Regex(r'''foobar''', 'i', 'any')
            ),
        ]
        for inp, outp in tests:
            self.assertEqual(parse_query(inp), outp)

    def test_fail(self):
        with self.assertRaises(LexerError):
            parse_query(r'"foo" 123')
        with self.assertRaises(NoParseError):
            parse_query(r'"foo" and')
