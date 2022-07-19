#!/usr/bin/env python
# Author: Mike Cuoco
# Created on: 7/14/22, 9:35 AM
#
# Description: add SLAV-seq tags to a BAM file
# Usage: python add_tags.py {input.bam} {output.bam}

import pysam
import pyfaidx
import argparse
from Bio import Align
from Bio.Seq import Seq
from pathlib import Path

CONSENSUS = Seq('ATGTACCCTAAAACTTAGAGTATAATAAA')
PREFIX_LENGTH = len(CONSENSUS)+2
R1_FLANK_LENGTH=750
R2_FLANK_LENGTH=PREFIX_LENGTH
SOFT_CLIP_LENGTH_THRESHOLD=5

# We add the custom flags below to the bam file:
# YR: a 0 or 1 valued integer, indicating if the pair is a reference read pair or not. 
#   1: the beginning of read 2 contains less than SOFT_CLIP_LENGTH_THRESHOLD soft clip letters, and it is in proper pair with read 1. 
#   0: otherwise.
# YS: size of the clip
# YA: score of semiglobal alignment between the first $prefix_length bp of read2 and $consensus.
# YG: score of prefix vs flanking genome.

def calculate_yr(read, soft_clip_threshold):
    if read.is_unmapped: return 0
    if not read.is_proper_pair: return 0

    is_clipped = False
    for cigar in read.cigartuples:
        is_clipped = True if cigar[0] == 4 and cigar[1] >= soft_clip_threshold else is_clipped

    if not read.is_reverse and is_clipped : return 0
    if read.is_reverse and is_clipped: return 0
    
    return 1

def calculate_ys(read):

    clip_len = 0
    for cigar in read.cigartuples:
        clip_len = cigar[1] if cigar[0] == 4 else clip_len
            
    if not read.is_reverse and clip_len: 
        return clip_len

    end_clip_len = read.cigartuples[-1][0] if read.cigartuples[-1][0] == 4 else 0
    if read.is_reverse and end_clip_len: 
        return end_clip_len

    return 0

def calculate_yg_from_r1(read, prefix, extra_flank):
    if read.is_unmapped: return 0
    if read.is_reverse:
        b = max(0, read.reference_end-extra_flank)
        e = read.reference_end + read.reference_length
        flank = ref.get_seq(str(read.reference_name), b, e)
    else:
        b = max(0, read.reference_end)
        e = read.reference_end + extra_flank
        flank = ref.get_seq(str(read.reference_name), b, e).reverse.complement

    return aligner.score(flank.seq, prefix)

def calculate_yg_from_r2(read, prefix, prefix_length, extra_flank):

    if read.is_reverse:
        e = read.reference_end + read.reference_length - 1 + prefix_length + extra_flank
        b = max(0, e-2*extra_flank-prefix_length+1)
        flank = ref.get_seq(str(read.reference_name), b, e).reverse.complement
    else:
        e = read.reference_end - 1 + extra_flank
        b = max(0,read.reference_end - prefix_length - extra_flank)
        flank = ref.get_seq(str(read.reference_name), b, e)

    return aligner.score(flank.seq, prefix)

def calculate_ya(consensus, prefix):
    return aligner.score(consensus, prefix)

def parse_args():

    parser = argparse.ArgumentParser(description="Add tags to bam file alignments for SLAV-seq random forest training")
    parser.add_argument(
        "-o", "--out", type=Path, required=True, help="the filename to save the tagged BAM file"
    )
    parser.add_argument(
        "-b", "--bam", type=Path, required=True, help="input BAM file to tag"
    )
    parser.add_argument(
        "-r", "--reference", type=Path, required=True, help="reference fasta file"
    )
    args = parser.parse_args()

    return args

def main():
    
    # get arguments
    args = parse_args()

    # sort alignment by read names and then open with pysam
    # TODO: maybe remove this sort step and instead pipe samtools output into this script. 
    sorted_bam = args.bam.with_suffix(".namesorted.sam")
    pysam.sort("-n", str(args.bam), "-o", str(sorted_bam))
    reads = pysam.AlignmentFile(sorted_bam, "rb")
    out = pysam.AlignmentFile(str(args.out), "wb", template=reads)

    # read in the reference
    global ref
    ref = pyfaidx.Fasta(str(args.reference))

    # setup pairwise aligner
    # used scores from from https://github.com/apuapaquola/gapafim/blob/main/Gapafim/sw.h
    # note: algorithm is automatically chosen based on the scoring matrix, need to check if scores match perl implementation
    global aligner
    aligner = Align.PairwiseAligner()
    aligner.alphabet = "ACGTNactgn"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1
    aligner.mode = "local"

    read_id = ""
    for count, read in enumerate(reads):  
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
            yr = calculate_yr(read, SOFT_CLIP_LENGTH_THRESHOLD)
            print(f"yr = {yr}")
            ys = calculate_ys(read)
            print(f"ys = {ys}")

        # get what should be LINE-1 sequence
        seq2 = Seq(r2.seq).reverse_complement() if r2.is_reverse else Seq(r2.seq)
        prefix = seq2[:PREFIX_LENGTH]

        # calculate YG
        if r2.is_unmapped and r2.is_proper_pair:
            yg = calculate_yg_from_r2(r2, prefix, PREFIX_LENGTH, R2_FLANK_LENGTH)
        else:
            yg = calculate_yg_from_r1(r1, prefix, R1_FLANK_LENGTH)
        print(f"yg = {yg}")

        # calculate YA
        ya = calculate_ya(CONSENSUS, prefix)
        print(f"ya = {ya}")
        # save YR, YS, YG, YA to read, write to file
        for r in [r1,r2]:
            if r.is_read1:                
                r.set_tag("Y2",str(seq2),"Z")
            r.set_tag("YR", int(yr), "i")
            r.set_tag("YS", int(ys), "i")
            r.set_tag("YG", int(yg), "i")
            r.set_tag("YA", int(ya), "i")
            out.write(r)    

    reads.close()
    out.close()

if __name__ == "__main__":
    main()