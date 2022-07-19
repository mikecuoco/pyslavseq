#!/usr/bin/env python
__author__ = "Michael Cuoco"

from pathlib import Path
import pyfaidx
import pysam
import pytest

# setup class for testing
class TestAddTagsPerl:
    def __init__(self):
        from pyslavseq.tags.add_tags import setup_aligner, set_parameters
        self.aligner = setup_aligner()
        self.parameters = set_parameters()

    def bam(self, bam_file):
        self.bam = pysam.AlignmentFile(str(bam_file), 'rb')

    def ref(self,ref_file):
        self.ref = pyfaidx.Fasta(str(ref_file))

# define fixture 
@pytest.fixture
def test():
    test = TestAddTagsPerl()
    test.bam(bam_file = Path(__file__).parent / 'test.bam')
    test.ref(ref_file = Path(__file__).parent / "hs37d5_hg19_chr22.fa")
    return test

@pytest.fixture
def read(test):
    read_id = ""
    for count, read in enumerate(test.bam):  
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
            yield r2

    test.bam.close()


def test_calculate_yr(read, SOFT_CLIP_LENGTH_THRESHOLD):
    from pyslavseq.tags.add_tags import calculate_yr
    yr = calculate_yr(read, SOFT_CLIP_LENGTH_THRESHOLD)
    assert int(yr) == read.get_tag('YR')

def test_calculate_ys(read):
    from pyslavseq.tags.add_tags import calculate_ys
    yr = calculate_ys(read)
    assert int(yr) == read.get_tag('YR')
