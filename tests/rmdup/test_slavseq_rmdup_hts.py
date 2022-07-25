#!/usr/bin/env python
__author__ = "Rohini Gadde"

import os
import filecmp
from .. import rmdup as _myself

class TestRmdup:
    
    # Created by running samtools view -s 123.01 on USD07_A1_S87_MDA.bam
    bam_in = _myself.__path__[0] + "/no_dedup.bam"

    # Created by running slavseq_rmdup_hts.pl, taking sample.bam as input
    bam_out_perl = _myself.__path__[0] + "/perl_dedup.bam"

    # Created by running slavseq_rmdup_hts.py (see below)
    bam_out_py = _myself.__path__[0] + "/py_dedup.bam"

    # Test prio_pair_rmdup function output
    def test_prio_pair_rmdup(self):
        from pyslavseq.rmdup.slavseq_rmdup_hts import prio_pair_rmdup

        os.environ['LC_ALL'] = "C"
        prio_pair_rmdup(self.bam_in, _myself.__path__[0] + "/py_all_fields.txt")
    
        check = filecmp.cmp(
            _myself.__path__[0] + "/perl_all_fields.txt", 
            _myself.__path__[0] + "/py_all_fields.txt", shallow=False)
        
        assert check
        
        os.remove(_myself.__path__[0] + "/py_all_fields.txt")

    # Compare headers of perl and python outputs
    def test_header(self):
        import subprocess

        if not os.path.exists(self.bam_out_py):
            # Change directory to git root directory to run script
            os.chdir(_myself.__path__[0])
            os.chdir("../../")

            subprocess.run(["pyslavseq/rmdup/slavseq_rmdup_hts.py", "-b", self.bam_in, "-o", self.bam_out_py], check=True)
            # os.system("pyslavseq/rmdup/slavseq_rmdup_hts.py -b " + bam_in + " -o " + bam_out_py)

        perl_header_path = _myself.__path__[0] + "/perl_header_sq.txt"
        py_header_path = _myself.__path__[0] + "/py_header_sq.txt"
        
        py_header_file = open(py_header_path, "w+")
        bam_py_file = open(self.bam_out_py, "r")

        p1 = subprocess.Popen(["samtools", "view", "-H"], stdin=bam_py_file, stdout=subprocess.PIPE)
        p2 = subprocess.run(["grep", "@SQ"], stdin=p1.stdout, stdout=py_header_file, check=True)
        
        bam_py_file.close()
        py_header_file.close()

        check_header = filecmp.cmp(perl_header_path, py_header_path, shallow=False)
        assert check_header

        os.remove(py_header_path)

    # Compare entire outputs of perl and python scripts
    def test_output(self):
        import subprocess

        if not os.path.exists(self.bam_out_py):
            # Change directory to git root directory to run script
            os.chdir(_myself.__path__[0])
            os.chdir("../../")

            subprocess.run(["pyslavseq/rmdup/slavseq_rmdup_hts.py", "-b", self.bam_in, "-o", self.bam_out_py], check=True)

        sam_out_perl = _myself.__path__[0] + "/perl_dedup.sam"
        sam_out_py = _myself.__path__[0] + "/py_dedup.sam"

        subprocess.run(["samtools", "view"], stdin=open(self.bam_out_py, "r"), stdout=open(sam_out_py, "w+"), check=True)

        check = filecmp.cmp(sam_out_perl, sam_out_py, shallow=False)

        assert check == True
        
        os.remove(self.bam_out_py)
        os.remove(sam_out_py)
