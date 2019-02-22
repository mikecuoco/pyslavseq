#!/usr/bin/env python
__author__ = "Apuã Paquola"

from .. import tabix as _myself


def test_tabixsamwithpolya():
    import pysam
    from pyslavseq.features.TabixSamWithPolyA import TabixSamWithPolyA

    tabixsam = TabixSamWithPolyA(pysam.Tabixfile(_myself.__path__[0] + "/dummy.bgz"))
    for r in tabixsam.fetch('chr3', 19416000, 19416750):
        if 'r2_poly_a_length' in r:
            assert r['r2_poly_a_length'] == 17
            break
