#!/usr/bin/env python
__author__ = "Apuã Paquola"


from .. import genome as _myself


def test_ttaaaa():
    from pyslavseq.features.ttaaaa import ENSearch
    from pyslavseq.genome import Genome

    mini = Genome(_myself.__path__[0] + '/mini.chromsizes')

    e = ENSearch(mini, _myself.__path__[0] + '/mini.fa', 50, 15)
    r = e.pos_and_score('chr2', 63, 1)
    assert r == ('TTAAAA', -20, 1423)
