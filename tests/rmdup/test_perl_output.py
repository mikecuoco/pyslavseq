#!/usr/bin/env python
__author__ = "Rohini Gadde"

from .. import rmdup as _myself

# Compare contents of perl and python scripts
def test_perloutput():
    import filecmp
    from pyslavseq.rmdup.slavseq_rmdup_hts import *

    bam_in = _myself.__path__[0] + "sample.bam"
    bam_out_py = _myself.__path__[0] + "py_dedup.bam"
    bam_out_perl = _myself.__path__[0] + "perl_dedup.bam"

    main(bam_in, bam_out_py)

    # perl = open(_myself.__path__[0] + "perl_dedup.bam", "rb")
    # py = open(_myself.__path__[0] + "py_dedup.bam", "rb")

    check = filecmp.cmp(bam_out_perl, bam_out_py, shallow=False)

    assert check == True