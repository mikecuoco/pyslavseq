#!/usr/bin/env python
__author__ = "Apuã Paquola"

import pysam
import MOODS.tools
import MOODS.scan

from ..genome import Interval
from ..util import rc
import numpy as np
import functools

def scan_best_hits_dna(seq, matrices, target, iterations, MULT, LIMIT_MULT, window_size=7):
    bg = MOODS.tools.bg_from_sequence_dna(seq, 0.01)

    # try to guess a good initial threshold
    p = ((1 + MULT) / 3.0) * target / len(seq)

    if (LIMIT_MULT < MULT):
        MULT = LIMIT_MULT

    m_len = len(matrices)

    current = [0] * m_len
    upper = [0] * m_len
    lower = [0] * m_len

    ok = [0] * m_len
    remaining = m_len

    # first we check which matrices can produce less than target*LIMIT_MULT hits in the first place
    for i in range(0, m_len):
        upper[i] = MOODS.tools.max_score(matrices[i], 4) - MOODS.tools.min_delta(matrices[i])/2
        
    scanner = MOODS.scan.Scanner(window_size)
    scanner.set_motifs(matrices, bg, upper)
    
    results = scanner.scan_max_hits(seq, LIMIT_MULT * target)
    counts = [len(n) for n in results]
    
    for i in range(0, m_len):
        hits = counts[i]
        # upper bound is a good limit
        if (hits >= target and hits < LIMIT_MULT * target):
            ok[i] = True
            current[i] = upper[i]
            remaining -= 1
            
        # even consensus sequences have too many hits
        # let's not return any hits for this matrix
        elif (hits >= LIMIT_MULT * target):
            current[i] = MOODS.tools.max_score(matrices[i], 4) + 1.0
            ok[i] = True
            remaining -= 1

    # search for the rest of matrices
    # we can assume that some threshold between the consensus score and
    # max score is what we want
    for i in range(0, m_len):
        if not ok[i]:
            current[i] = MOODS.tools.threshold_from_p(matrices[i], bg, p, 4)
            lower[i] = MOODS.tools.min_score(matrices[i], 4) - 1.0

    iteration = 0
    # binary search for all matrices in parallel
    while (remaining > 0):
        iteration += 1

        it_indices = []
        it_matrices = []
        it_thresholds = []

        for i in range(0, m_len):
            if not ok[i]:
                it_indices.append(i)
                it_matrices.append(matrices[i])
                it_thresholds.append(current[i])
        
        scanner.set_motifs(it_matrices, bg, it_thresholds)
        
        results = scanner.scan_max_hits(seq, LIMIT_MULT * target)
        counts = [len(n) for n in results]

        for j in range(0, len(it_indices)):
            i = it_indices[j]
            hits = counts[j]
            threshold = current[i]

            if (hits >= target and hits < MULT * target):
                ok[i] = True
                remaining -= 1

            # too many hits, INCREASE the threshold & lower bound
            elif (hits >= MULT*target):      
                current[i] = (upper[i] + threshold)/2.0
                lower[i] = threshold
        
            # too few this, DECREASE the threshold & upper bound
            else:   
                current[i] = (lower[i] + threshold)/2.0
                upper[i] = threshold
            
            # failsafe
            if ((current[i] == threshold or iteration == iterations) and not ok[i]):
                ok[i] = 1
                if (hits >= LIMIT_MULT * target):
                    current[i] = upper[i]
                    
                remaining -= 1

    scanner.set_motifs(matrices, bg, current)
        
    return scanner.scan(seq); 
                
    # return []

class ENSearch:
    """ Search for a L1 endonuclease motif given by a PWM from Jurka 1997.
    """

    def __init__(self, genome, genome_fasta_file, left_flank, right_flank):
        # Matrix from PMID 9050872 Jurka 1997
        self.matrix = [[60, 71, 279, 248, 238, 241],
                       [34, 37, 3, 4, 9, 14],
                       [43, 26, 32, 72, 72, 46],
                       [207, 210, 30, 20, 25, 43]]

        self.genome = genome
        self.fa = pysam.Fastafile(genome_fasta_file)
        self.left_flank = left_flank
        self.right_flank = right_flank
        # self.threshold = MOODS.tools.threshold_from_p(self.matrix, MOODS.tools.flat_bg(4), 0.2)

    @functools.lru_cache(maxsize=1024, typed=False)
    def pos_and_score(self, chrom, pos, te_strand):
        if pos is np.nan:
            return (np.nan, np.nan, np.nan)

        if te_strand == 1:
            # if te is in the + strand
            start = pos - self.left_flank
            end = pos + self.right_flank
        else:
            start = pos - self.right_flank
            end = pos + self.left_flank

        iv = self.genome.fit_interval(Interval(chrom, start, end))

        s = self.fa.fetch(iv.chrom, iv.start, iv.end).upper()

        if te_strand == 1:
            zeropos = pos - iv.start
        else:
            s = rc(s)
            zeropos = iv.end - pos
            # print("\n>> ", chrom, pos, te_strand,  file=sys.stderr, flush=True)

        # print(">> ", s,  file=sys.stderr, flush=True)

        def moods_results(s, matrix, left_flank):

            # This is a hack to work around a bug in MOODS. The last 3
            # parameters are: int iterations = 10, unsigned int MULT = 2,
            # size_t LIMIT_MULT = 10
            # Setting MULT and LIMIT_MULT to the size of the string fixes it

            results = scan_best_hits_dna(s, [matrix], 1, 10, len(s), len(s))
            if len(results) == 0:
                en_pos = - self.left_flank
                en_score = 0
                motif = ''
            else:
                # print(">> ", len(results[0]),  file=sys.stderr, flush=True)
                spos = results[0][0].pos
                motif = s[spos:spos + len(matrix[0])]
                en_pos = spos - zeropos
                en_score = int(results[0][0].score)

            return motif, en_pos, en_score

        return moods_results(s, self.matrix, self.left_flank)
