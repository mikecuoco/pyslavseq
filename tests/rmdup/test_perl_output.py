#!/usr/bin/env python
__author__ = "Rohini Gadde"

from .. import rmdup as _myself

# Compare contents of perl and python scripts
def test_perloutput():
    import os
    import subprocess
    import filecmp

    # Created by running samtools view -s 123.01 on USD07_A1_S87_MDA.bam
    bam_in = _myself.__path__[0] + "/no_dedup.bam"

    # Created by running slavseq_rmdup_hts.pl, taking sample.bam as input
    bam_out_perl = _myself.__path__[0] + "/perl_dedup.bam"

    # Created by running slavseq_rmdup_hts.py (see below)
    bam_out_py = _myself.__path__[0] + "/py_dedup.bam"

    os.system("../../pyslavseq/rmdup/slavseq_rmdup_hts.py -b " + bam_in + " -o " + bam_out_py)
    # p1 = subprocess.Popen(["git", "rev-parse", "--show-toplevel"], stdout=subprocess.PIPE, text=True)
    # p2 = subprocess.run(
        # ["$1/pyslavseq/rmdup/slavseq_rmdup_hts.py"), "-b", bam_in, "-o", bam_out_py], 
        # stdin=p1.stdout, check=True, text=True)
    # perl = open(_myself.__path__[0] + "perl_dedup.bam", "rb")
    # py = open(_myself.__path__[0] + "py_dedup.bam", "rb")

    check = filecmp.cmp(bam_out_perl, bam_out_py, shallow=False)

    assert check == True
    # assert 1 == 1