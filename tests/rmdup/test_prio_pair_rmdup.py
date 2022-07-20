#!/usr/bin/env python
__author__ = "Rohini Gadde"

from .. import rmdup as _myself

def test_perloutput():
    from pyslavseq.rmdup.slavseq_rmdup_hts import prio_pair_rmdup

    # Created by running samtools view -s 123.01 on USD07_A1_S87_MDA.bam
    bam_in = _myself.__path__[0] + "/no_dedup.bam"

    prio_pair_rmdup(bam_in, _myself.__path__[0] + "/selected.txt")

    selected = open(_myself.__path__[0] + "/selected.txt", "r")
    ids = [i.strip().split("\t")[0] for i in selected]

    perl_dedup = open(_myself.__path__[0] + "/perl_dedup.sam", "r")
    for line in perl_dedup:
        if line.strip().split("\t")[0] not in ids:
            assert False
    
    os.remove(_myself.__path__[0] + "/selected.txt")
