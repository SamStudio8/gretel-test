"""Given a BAM, a contig, and an ending genomic position, aggressively call for
variants and generate a VCF."""
import sys

import numpy as np
import pysam
import argparse

parser = argparse.ArgumentParser('Aggressively call for variants and generate a VCF')
parser.add_argument('-b', help='bam file', required=True)
parser.add_argument('-r', help='reference sequence (contig) to match', required=True)
parser.add_argument('-s', help='start', type=int, default=0)
parser.add_argument('-e', help='end', type=int)
parser.add_argument('-d', help='minimum depth to call (default = 0)', type=int, default=0)
args = parser.parse_args()

bam = pysam.AlignmentFile(args.b)
counts = np.array(bam.count_coverage(reference=args.r, start=args.s, end=args.e, quality_threshold=0, read_callback='nofilter'))

COUNT_SENSITIVITY = args.d

vcf_h = [
    "##fileformat=VCFv4.2",
]
vcf = []

sites = (counts > COUNT_SENSITIVITY).sum(axis=0)
for i, s in enumerate(sites):
    if s > 1:
        vcf.append([
            args.r,
            i+1+args.s,
            '.',
            'G',
            'A,T,C',
            0,
            '.',
            "INFO"
        ])


for r in vcf_h:
    print(r)
for r in vcf:
    print("\t".join([str(s) for s in r]))
