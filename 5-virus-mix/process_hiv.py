import os
import sys

import pysam

HIVG_DIR = sys.argv[1]
TRIM_SIZE = int(sys.argv[2])
TRIM_SIZE_e = int(sys.argv[3])

try:
    out_fasta = pysam.FastaFile(HIVG_DIR + '/' + "/out.fasta")
except IOError:
    sys.exit(1)
for ref in out_fasta.references:
    print ">%s" % ref
    seq = out_fasta.fetch(reference=ref)
    print seq[TRIM_SIZE : len(seq)-TRIM_SIZE_e].replace("-", "")
