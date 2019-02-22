#!/usr/bin/env python
__author__ = "Apuã Paquola"

from .. import genome as _myself


def test_genome():
    from pyslavseq.genome import Genome
    mini = Genome(_myself.__path__[0] + '/mini.chromsizes')
    assert mini.chromsizes['chr1'] == 77
    assert mini.chromsizes['chr2'] == 73

