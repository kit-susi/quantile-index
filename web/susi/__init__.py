from subprocess import Popen, PIPE
import os
import sys

from query_parser import And, Or, Not, String, Regex, Any, parse_query, LexerError, NoParseError
from repoze.lru import LRUCache
import susi_native


class Collection(object):
    def __init__(self, directory, cache_size=128):
        self.dir = directory
        self.index = susi_native.load(self.dir)
        self._cache = LRUCache(cache_size)

    def lookup(self, pattern, k=10, snippet_size=1000000):
        key = (pattern, k, snippet_size)
        cached_result = self._cache.get(key)
        if cached_result is not None:
            return cached_result

        if len(pattern) < 3:
            return []
        res = susi_native.search(self.index, pattern.encode('latin1'), k, snippet_size)
        self._cache.put(key, res)
        return res
