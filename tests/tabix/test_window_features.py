#!/usr/bin/env python
__author__ = "Apuã Paquola"

from .. import tabix as _myself


def test_windowfeatures():
    import pysam
    from pyslavseq.features.TabixSamWithPolyA import TabixSamWithPolyA
    from pyslavseq.features.WindowFeatures import WindowFeatures

    tabixsam = TabixSamWithPolyA(pysam.Tabixfile(_myself.__path__[0] + "/dummy.bgz"))

    # TODO: add this test with a proper fasta file that corresponds to the bgz file
    # from .ttaaaa import ENSearch
    # e = ENSearch(hg19, 'hg19.fa', 55, 15)

    e = None
    wf = WindowFeatures(tabixsam, 'chr3', 19502755, 19502755+750, 40, 20, 15, 20, '3', e)
    f = wf.features()
    assert f['all_reads.count'] == 12
    assert f['yayg_reads.count'] == 0
    assert f['all_reads.r2_peak_position'] == 19503197
    assert f['all_reads.r2_peak_count'] == 11

