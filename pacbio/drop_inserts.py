import pysam
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('bam', help='path to bam')
parser.add_argument('fasta', help='reference fasta')
parser.add_argument('contig', help='contig name')
parser.add_argument('snps', help='new line delimited 1-pos based snps')
args = parser.parse_args()

bam_fh = pysam.AlignmentFile(args.bam)
fasta_fh = pysam.FastaFile(args.fasta)
ref_seq = fasta_fh.fetch(reference=args.contig)

snps = [int(x)-1 for x in open(args.snps).readlines()] # snp positions are 1-pos, pysam is 0-pos

reads = bam_fh.fetch(args.contig)
for r in reads:
    aligned_tuples = r.get_aligned_pairs()
    seq = []

    if r.query_alignment_length <= 0.8 * len(ref_seq):
        # Drop CCS read if its smaller than 80% of the length of the full haplotype
        continue

    for at in aligned_tuples:
        qry = at[0] # position along CCS read
        ref = at[1] # position aginst reference haplotype

        if ref in snps:
            if not qry:
                # Deletion in haplotype, output no SNP here
                continue
            # else use the SNP as seen on the CCS query_sequence
            seq.append(r.query_sequence[qry])
        else:
            if not ref:
                # Insertion against reference, discard it
                continue
            # else this site is not an original Gretel SNP so fill it in with the reference
            # recall ref is pysam 0-pos
            seq.append(ref_seq[ref])

    print(">%s/compressed\n%s" % (r.query_name, "".join(seq)))
