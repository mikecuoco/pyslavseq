#!/usr/bin/env python
__author__ = "Michael Cuoco"

from pathlib import Path
import pyfaidx
import pysam
import pytest

# setup class for testing
from pyslavseq.tags.add_tags import setup_aligner, set_parameters
aligner = setup_aligner()
CONSENSUS, PREFIX_LENGTH, R1_FLANK_LENGTH, R2_FLANK_LENGTH, SOFT_CLIP_LENGTH_THRESHOLD = set_parameters()

ref_file = Path(__file__).parent / "hs37d5_hg19_chr22.fa"
bam_file = Path(__file__).parent / 'test.bam'
bam = pysam.AlignmentFile(str(bam_file), 'rb')
ref = pyfaidx.Fasta(str(ref_file))

def test_calculate_yr(read, SOFT_CLIP_LENGTH_THRESHOLD):
    from pyslavseq.tags.add_tags import calculate_yr
    yr = calculate_yr(read, SOFT_CLIP_LENGTH_THRESHOLD)
    assert int(yr) == read.get_tag('YR')

def test_calculate_ys(read):
    from pyslavseq.tags.add_tags import calculate_ys
    ys = calculate_ys(read)
    assert int(ys) == read.get_tag('YS')

read_id = ""
for count, read in enumerate(bam):  
    # skip secondary and supplementary alignments
    if read.is_secondary or read.is_supplementary: continue
    print(f"Processing read #{count+1} with id: {read.query_name}")
    if read.is_read1:
        r1 = read
    else:
        r2 = read

    # check if read is a mate of the previous read
    if read.query_name != read_id:
        read_id = read.query_name
        continue

    # calculate YR and YS if read has cigarstring
    if r2.cigarstring:
        test_calculate_yr(r2, SOFT_CLIP_LENGTH_THRESHOLD)
        test_calculate_ys(r2)


