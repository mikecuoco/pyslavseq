#!/usr/bin/env python
__author__ = "Michael Cuoco"

from pathlib import Path
import pyfaidx
import pysam
from Bio import Align
from Bio.Seq import Seq
import pytest
from pyslavseq.tags.add_tags import *

aligner = setup_aligner()
CONSENSUS, PREFIX_LENGTH, R1_FLANK_LENGTH, R2_FLANK_LENGTH, SOFT_CLIP_LENGTH_THRESHOLD = set_parameters()
    
bam_file = Path(__file__).parent / 'test.sam'
bam = pysam.AlignmentFile(str(bam_file), 'rb')
ref_file = Path(__file__).parent / 'hs37d5_hg19_chr22.fa'
ref = pyfaidx.Fasta(str(ref_file))

read_id = ""
yr_inputs, ys_inputs, ya_inputs, y2_inputs, yg_r1_inputs, yg_r2_inputs = [],[],[],[],[],[]
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
        yr_inputs.append(([r2, SOFT_CLIP_LENGTH_THRESHOLD], read.get_tag('YR')))
        ys_inputs.append((r2, read.get_tag('YS')))

    # get what should be LINE-1 sequence
    seq2 = Seq(r2.seq).reverse_complement() if r2.is_reverse else Seq(r2.seq)
    prefix = seq2[:PREFIX_LENGTH]

    # test Y2
    y2_inputs.append((seq2, r1.get_tag('Y2')))


    # calculate YG
    if r2.is_unmapped and r2.is_proper_pair:
        yg_r2_inputs.append(([r2, prefix, PREFIX_LENGTH, R2_FLANK_LENGTH, aligner, ref], read.get_tag('YG')))
    else:
        yg_r1_inputs.append(([r1, prefix, R1_FLANK_LENGTH, aligner, ref], read.get_tag('YG')))

    # calculate YA
    ya_inputs.append(([CONSENSUS, prefix, aligner], read.get_tag('YA')))

bam.close()

@pytest.mark.parametrize('inputs, expected', yr_inputs)
def test_calculate_yr(inputs, expected):
    yr = calculate_yr(*inputs)
    assert int(yr) == expected

@pytest.mark.parametrize('inputs, expected', ys_inputs)
def test_calculate_ys(inputs, expected):
    ys = calculate_ys(inputs)
    assert int(ys) == expected

@pytest.mark.parametrize('inputs, expected', y2_inputs)
def test_y2(inputs, expected):
    assert str(inputs) == str(expected)

@pytest.mark.parametrize('inputs, expected', yg_r2_inputs)
def test_calculate_yg_from_r2(inputs, expected):
    yg = calculate_yg_from_r2(*inputs)
    assert int(yg) == expected

@pytest.mark.parametrize('inputs, expected', yg_r1_inputs)
def test_calculate_yg_from_r1(inputs, expected):
    yg = calculate_yg_from_r1(*inputs)
    assert int(yg) == expected

@pytest.mark.parametrize('inputs, expected', ya_inputs)
def test_calculate_ya(inputs, expected):
    ya = calculate_ya(*inputs)
    assert int(ya) == expected


