"""Given a directory of UUIDs, each containing a file that ends in .sort.bam, and
a FASTA file of known haplotypes, create a heatmap of how the reads across each of
those BAMs, from each of the known genes, maps to the pseudo-reference."""
import os
import glob
import sys

import numpy as np
import pysam

import matplotlib.pyplot as plt

num_bams = 0
#for bam_fn in [f for f in os.listdir(sys.argv[1]) if os.path.join(sys.argv[1], f).endswith(".sort.bam")]:
for bam_fn in [f for f in glob.glob(sys.argv[1] + "/*.sort.bam")]:
    print bam_fn
    bam = pysam.AlignmentFile(bam_fn)
    if num_bams == 0:
        ref_name = bam.references[0]
        ref_len = bam.lengths[0]
        gene_names = list(pysam.FastaFile(sys.argv[2]).references)
        counts = np.zeros( (ref_len, len(gene_names)) )
    num_bams += 1

    for p_col in bam.pileup(reference=ref_name, start=0, end=ref_len):
        for pc_read in p_col.pileups:
            g = (pc_read.alignment.query_name.split("_%")[0].replace("HOOT_", ""))
            try:
                counts[p_col.pos][gene_names.index(g)] += 1
            except:
                print g, gene_names
                import sys; sys.exit(3)

#print(counts)
print num_bams
counts /= num_bams

plt.pcolor(counts.transpose())
plt.yticks([x+0.5 for x in range(len(gene_names))], gene_names)
plt.colorbar()
plt.axis([0, ref_len, 0, len(gene_names)])
plt.xticks(range(0, ref_len, 25))
plt.title("Read Coverage Rates\n" + sys.argv[1])
plt.show()
