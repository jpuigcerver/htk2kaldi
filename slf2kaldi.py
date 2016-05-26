#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
# The MIT License (MIT)
#
# Copyright (c) 2016 Joan Puigcerver <joapuipe@prhlt.upv.es>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import argparse
import math
import re
import sys

FLOAT_RE = '(?:[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)'
S_RE='((?:S|START)=(?P<S>\d+))'
E_RE='((?:E|END)=(?P<E>\d+))'
W_RE='((?:W|WORD)=(?P<W>\S+))'
v_RE='((?:var|v)=(?P<v>\d+))'
d_RE='((?:div|d)=(?P<d>\S+))'
a_RE='((?:acoustic|a)=(?P<a>%s))' % FLOAT_RE
n_RE='((?:ngram|n)=(?P<n>%s))' % FLOAT_RE
l_RE='((?:language|l)=(?P<l>%s))' % FLOAT_RE

ARC_RE = re.compile(
    '^J=\d+(?:\s+(?:%s|%s|%s|%s|%s|%s|%s|%s))+$' % (
        S_RE, E_RE, W_RE, v_RE, d_RE, a_RE, n_RE, l_RE))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Convert a SLF lattice into Kaldi\'s format.')
    parser.add_argument('input', type=argparse.FileType('r'),
                        nargs='?', default=sys.stdin,
                        help='Input model in HTK text format')
    parser.add_argument('output', type=argparse.FileType('w'),
                        nargs='?', default=sys.stdout,
                        help='Output model in Kaldi text format')
    args = parser.parse_args()

    def err_msg(m, f, l):
        return '%s at %s:%d' % (m, f.name, l)

    N = set()   # Set of nodes
    O = set()   # Set of nodes with output arcs

    lnum = 0
    for line in args.input:
        lnum += 1
        line = line.strip()
        m = ARC_RE.match(line)
        if not m: continue
        S = m.group('S')
        E = m.group('E')
        W = m.group('W')
        a = 0.0 if not m.group('a') else -float(m.group('a'))
        l = 0.0 if not m.group('l') else -float(m.group('l'))
        n = 0.0 if not m.group('n') else -float(m.group('n'))
        assert S is not None, err_msg('Expected field S', args.input, lnum)
        assert E is not None, err_msg('Expected field E', args.input, lnum)
        assert W is not None, err_msg('Expected field W', args.input, lnum)
        print S, E, W, '%g,%g,%s' % (l + n, a, W)
        N.add(int(S))
        N.add(int(E))
        O.add(int(S))

    for n in N.difference(O):
        print n, '0,0'
