"""Given a BAM, a contig, and an ending genomic position, aggressively call for
variants and generate a VCF."""
import sys

import numpy as np
import pysam

bam = pysam.AlignmentFile(sys.argv[1])
counts = np.array(bam.count_coverage(reference=sys.argv[2], start=0, end=int(sys.argv[3]), quality_threshold=0, read_callback='nofilter'))

COUNT_SENSITIVITY = 0
try:
    COUNT_SENSITIVITY = int(sys.argv[4])
except IndexError:
    pass

vcf_h = [
    "##fileformat=VCFv4.2",
]
vcf = []

sites = (counts > COUNT_SENSITIVITY).sum(axis=0)
for i, s in enumerate(sites):
    if s > 1:
        vcf.append([
            sys.argv[2],
            i+1,
            '.',
            'G',
            'A,T,C',
            0,
            '.',
            "INFO"
        ])


for r in vcf_h:
    print r
for r in vcf:
    print "\t".join([str(s) for s in r])
