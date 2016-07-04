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

REGEX_FLOAT = r'(?:[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)'

class GMM(object):
    def __init__(self):
        self._mixes = []

    @property
    def dimension(self):
        return (len(self._mixes[0][0]) if len(self._mixes) > 0 else None)

    @property
    def num_mixtures(self):
        return (len(self._mixes) if self._mixes is not None else None)

    def add_mixture(self, m, v, w):
        assert isinstance(w, float), 'Mixture weight must '
        if not isinstance(m, tuple): m = tuple(m)
        if not isinstance(v, tuple): v = tuple(v)
        assert len(m) == len(v), \
            'Dimension of mean and variance vectors must agree'
        self._mixes.append((m, v, w))
        if len(self._mixes) > 1:
            assert len(self._mixes[0][0]) == len(self._mixes[-1][0]), \
                'All mixtures in a GMM must have the same dimension'

    def write(self, f):
        f.write('<DiagGMM>\n')
        f.write('<WEIGHTS> [ %s ]\n' % ' '.join(
            map(lambda x: str(x[-1]), self._mixes)))
        f.write('<MEANS_INVVARS> [\n')
        for (m, v, _) in self._mixes:
            for mi, vi in zip(m, v):
                f.write(' %g' % (mi / vi))
            f.write('\n')
        f.write(']\n')
        f.write('<INV_VARS> [\n')
        for (_, v, _) in self._mixes:
            for vi in v:
                f.write(' %g' % (1.0 / vi))
            f.write('\n')
        f.write(']\n')

        f.write('</DiagGMM>\n')


class HMM(object):
    def __init__(self):
        self._states = []
        self._pinit = None

    @property
    def dimension(self):
        return self._states[0][0].dimension if len(self._states) > 0 \
            else None

    def add_state(self, gmm):
        self._states.append([gmm, []])

    def set_initial(self, pinit):
        self._pinit = pinit

    def set_transition(self, i, ptrans):
        assert i >= 0 and i < len(self._states), 'Unknown state %d' % i
        if not isinstance(ptrans, tuple): ptrans = tuple(ptrans)
        self._states[i][1] = ptrans

    def write_transitions(self, f):
        if self._pinit[1] < 1.0:
            # Initial State
            f.write('<State> 0\n')
            for j, p in enumerate(self._pinit):
                if p <= 0.0: continue
                f.write('<Transition> %d %g\n' % (j, p))
            # Emitting states
            for s, (_, trans) in enumerate(self._states):
                f.write('<State> %d <PdfClass> %d\n' % (s + 1, s))
                for j, p in enumerate(trans):
                    if p <= 0.0: continue
                    f.write('<Transition> %d %g\n' % (j, p))
                f.write('</State>\n')
            # Final state
            f.write('<State> %d\n' % (len(self._states) + 1))
            f.write('</State>\n')
        else:
            # Emitting states
            for s, (_, trans) in enumerate(self._states):
                f.write('<State> %d <PdfClass> %d\n' % (s, s))
                for j, p in enumerate(trans):
                    if p <= 0.0: continue
                    f.write('<Transition> %d %g\n' % (j - 1, p))
                f.write('</State>\n')
            # Final state
            f.write('<State> %d\n' % len(self._states))
            f.write('</State>\n')

    def write_pdfs(self, f):
        for (gmm, _) in self._states:
            gmm.write(f)

    def get_states_pdfs(self):
        if self._pinit[1] < 1.0:
            return map(lambda x: (x + 1, x), xrange(len(self._states)))
        else:
            return map(lambda x: (x, x), xrange(len(self._states)))

    def get_transition_logprobs(self):
        logprobs = []
        if self._pinit[1] < 1.0:
            for p in self._pinit:
                if p <= 0.0: continue
                logprobs.append(math.log(p))
        for (_, trans) in self._states:
            for p in trans:
                if p <= 0.0: continue
                logprobs.append(math.log(p))
        return logprobs


class Model(object):
    def __init__(self):
        self._dim = None
        self._hmm = {}

    @property
    def dimension(self):
        return self._dim

    @dimension.setter
    def dimension(self, value):
        assert isinstance(value, int)
        self._dim = value

    def add_hmm(self, i, hmm):
        assert isinstance(i, int)
        assert self._dim == hmm.dimension
        self._hmm[i] = hmm

    def write(self, f):
        if isinstance(f, str) or isinstance(f, unicode):
            f = open(f, 'w')
        f.write('<TransitionModel>\n')
        f.write('<Topology>\n')
        for hmm_id, hmm in self._hmm.iteritems():
            f.write('<TopologyEntry>\n')
            f.write('<ForPhones> %d </ForPhones>\n' % hmm_id)
            hmm.write_transitions(f)
            f.write('</TopologyEntry>\n')
        f.write('</Topology>\n')
        # Triples
        triples, logprobs = [], [0]
        for hmm_id, hmm in self._hmm.iteritems():
            base_pdfs = len(triples)
            for (s, pdf) in hmm.get_states_pdfs():
                triples.append((hmm_id, s, pdf + base_pdfs))
            logprobs.extend(hmm.get_transition_logprobs())
        f.write('<Triples> %d\n' % len(triples))
        for (h, s, p) in triples:
            f.write('%d %d %d\n' % (h, s, p))
        f.write('</Triples>\n')
        f.write('<LogProbs>\n[ %s ]\n</LogProbs>\n' % \
                ' '.join(map(lambda x: str(x), logprobs)))
        f.write('</TransitionModel>\n')
        num_pdfs = max(map(lambda x: x[-1], triples)) + 1
        f.write('<DIMENSION> %d <NUMPDFS> %d\n' % (self.dimension, num_pdfs))
        for hmm in self._hmm.itervalues():
            hmm.write_pdfs(f)
        f.close()




class HTKParser(object):
    def __init__(self, symbols_table):
        assert isinstance(symbols_table, dict)
        self._n = 0
        self._line = None
        self._data = None
        self._model = None
        self._symbols = symbols_table

    def _next(self):
        if not self._data:
            return False
        self._n += 1
        if self._n < len(self._data):
            self._line = self._data[self._n]
            return True
        else:
            self._line = None
            return False

    def _msg(self, s):
        return '%s at line %d: "%s"' % (s, self._n + 1, self._line)

    def _read_dim(self):
        assert self._line is not None, 'Unexpected EOF'
        m = re.match('^<STREAMINFO> 1 ([0-9]+)$', self._line)
        if not m: return False
        self._model.dimension = int(m.group(1))
        return self._next()

    def _read_hmm(self):
        assert self._line is not None, 'Unexpected EOF'
        m = re.match('^~h "([0-9A-Za-z]+)"$', self._line)
        if not m: return False
        hmm_sym = m.group(1)
        hmm_id = self._symbols.get(hmm_sym, None)
        assert hmm_id is not None, self._msg(
            'HMM symbol \"%s\" is not found in the symbols table' % hmm_sym)
        hmm = HMM()
        # Read <BEGINHMM>
        assert self._next(), 'Unexpected EOF'
        assert re.match('^<BEGINHMM>$', self._line), \
            self._msg('Expected "<BEGINHMM>"')
        # Read <NUMSTATES>
        assert self._next(), 'Unexpected EOF'
        m = re.match('^<NUMSTATES> ([0-9]+)$', self._line)
        assert m, self._msg('Expected "<NUMSTATES> %%d"')
        num_states = int(m.group(1))
        assert num_states > 2, self._msg('Unexpected number of states')
        # Process all <STATE>s
        assert self._next(), 'Unexpected EOF'
        for s in xrange(2, num_states):
            # Check state ID
            m = re.match('^<STATE> %d$' % s, self._line)
            assert m, self._msg('Expected "<STATE> %d"' % s)
            # Read number of mixtures
            assert self._next(), 'Unexpected EOF'
            m = re.match('^<NUMMIXES> ([0-9]+)$', self._line)
            gmm = GMM()
            if not m:
                ## Special case: 1 Gaussian
                num_mixes = 1
                weight = 1.0
                # Check <MEAN>
                m = re.match('^<MEAN> %d$' % self._model.dimension, self._line)
                assert m, self._msg('Expected "<MEAN> %d"' % \
                                    self._model.dimension)
                # Read mean vector
                assert self._next(), 'Unexpected EOF'
                mean = tuple([float(x) for x in self._line.split()])
                assert len(mean) == self._model.dimension, \
                    self._msg('Invalid mean vector')
                # Check <VARIANCE>
                assert self._next(), 'Unexpected EOF'
                m = re.match('^<VARIANCE> %d$' % self._model.dimension,
                             self._line)
                assert m, self._msg('Expected "<VARIANCE> %d"' % \
                                    self._model.dimension)
                # Read variance vector
                assert self._next(), 'Unexpected EOF'
                var = tuple([ float(x) for x in self._line.split() ])
                assert len(var) == self._model.dimension, \
                    self._msg('Invalid variance vector')
                # Read Gconst. It is not used, since it can be obtained from
                # the mean and the variance vectors.
                assert self._next(), 'Unexpected EOF'
                m = re.match('^<GCONST> (%s)$' % REGEX_FLOAT, self._line)
                assert m, self._msg('Expected "<GCONST> %%f"')
                # Add mixture to the gmm
                gmm.add_mixture(mean, var, weight)
                assert self._next(), 'Unexpected EOF'
            else:
                ## General case, Gaussian Mixture Model
                num_mixes = int(m.group(1))
                assert num_mixes > 0, self._msg('Invalid number of mixtures')
                # Process all <MIXTURE>s
                while True:
                    # Read mixture weight
                    assert self._next(), 'Unexpected EOF'
                    m = re.match('^<MIXTURE> \d+ (%s)$' % REGEX_FLOAT,
                                 self._line)
                    if not m: break
                    weight = float(m.group(1))
                    assert weight > 0.0 and weight <= 1.0, \
                        self._msg('Invalid mixture weight')
                    # Check <MEAN>
                    assert self._next(), 'Unexpected EOF'
                    m = re.match('^<MEAN> %d$' % self._model.dimension,
                                 self._line)
                    assert m, self._msg('Expected "<MEAN> %d"' % \
                                        self._model.dimension)
                    # Read mean vector
                    assert self._next(), 'Unexpected EOF'
                    mean = tuple([float(x) for x in self._line.split()])
                    assert len(mean) == self._model.dimension, \
                        self._msg('Invalid mean vector')
                    # Check <VARIANCE>
                    assert self._next(), 'Unexpected EOF'
                    m = re.match('^<VARIANCE> %d$' % self._model.dimension,
                                 self._line)
                    assert m, self._msg('Expected "<VARIANCE> %d"' % \
                                        self._model.dimension)
                    # Read variance vector
                    assert self._next(), 'Unexpected EOF'
                    var = tuple([ float(x) for x in self._line.split() ])
                    assert len(var) == self._model.dimension, \
                        self._msg('Invalid variance vector')
                    # Read Gconst. It is not used, since it can be obtained from
                    # the mean and the variance vectors.
                    assert self._next(), 'Unexpected EOF'
                    m = re.match('^<GCONST> (%s)$' % REGEX_FLOAT, self._line)
                    assert m, self._msg('Expected "<GCONST> %%f"')
                    # Add mixture to the gmm
                    gmm.add_mixture(mean, var, weight)
            # Add state to the HMM
            hmm.add_state(gmm)
        # Process transitions, Read <TRASNP>.
        #assert self._next(), 'Unexpected EOF'
        m = re.match('^<TRANSP> %d$' % num_states, self._line)
        assert m, self._msg('Expected "<TRANSP> %d"' % self._n)
        # Process initial state
        assert self._next(), 'Unexpected EOF'
        t = tuple([float(x) for x in self._line.split()])
        assert len(t) == num_states, \
            self._msg('Unexpected number of transitions')
        hmm.set_initial(t)
        # Process emitting states
        for i in range(num_states - 2):
            assert self._next(), 'Unexpected EOF'
            t = tuple([ float(x) for x in self._line.split() ])
            assert len(t) == num_states, \
                self._msg('Unexpected number of transitions' % self._n)
            hmm.set_transition(i, t)
        # Process final state (just ignore it)
        assert self._next(), 'Unexpected EOF'
        t = tuple([ float(x) for x in self._line.split() ])
        assert len(t) == num_states, \
            'Unexpected number of transitions at line %d' % self._n
        # Check <ENDHMM>
        assert self._next(), 'Unexpected EOF'
        m = re.match('^<ENDHMM>$', self._line)
        self._model.add_hmm(hmm_id, hmm)
        return True

    def parse(self, f):
        if isinstance(f, str) or isinstance(f, unicode):
            f = open(f, 'r')
        self._model = Model()
        self._data = [line.strip() for line in f]
        f.close()
        assert len(self._data) > 0, 'Unexpected EOF'
        self._n = -1
        while self._next():
            if not self._read_dim() and not self._read_hmm():
                sys.stderr.write('IGNORED LINE: \"%s\"\n' % self._line)
        return self._model


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Convert a HMM/GMM model file from HTK to Kaldi format.')
    parser.add_argument('symbols_table', type=argparse.FileType('r'),
                        help='HMM symbols table')
    parser.add_argument('input', type=argparse.FileType('r'),
                        nargs='?', default=sys.stdin,
                        help='Input model in HTK text format')
    parser.add_argument('output', type=argparse.FileType('w'),
                        nargs='?', default=sys.stdout,
                        help='Output model in Kaldi text format')
    args = parser.parse_args()

    # Read symbols table
    symbols_table = {}
    for line in args.symbols_table:
        line = line.split()
        assert len(line) == 2, 'Wrong symbols table format'
        symbols_table[line[0]] = int(line[1])

    htk_parser = HTKParser(symbols_table)
    model = htk_parser.parse(args.input)
    model.write(args.output)
