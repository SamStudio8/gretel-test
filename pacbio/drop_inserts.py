import pysam
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('bam', help='path to bam')
parser.add_argument('fasta', help='reference fasta')
parser.add_argument('contig', help='contig name')
parser.add_argument('snps', help='new line delimited 1-pos based snps')
parser.add_argument('--onlysnps', action="store_true", help="only output bases on snp positions")
parser.add_argument('--useref', action="store_true", help="output reference bases at non-snp positions, not zmw")
parser.add_argument('--delchar', default='-', help="character to output for a deletion [default -]")
parser.add_argument('--padstart', action="store_true", help="pad the sequence with _ to maintain integrity of positions against the reference")
parser.add_argument('--bed', help="drop reads on a contig that start after or end before the provided co-ordinates [use 1-pos]")

args = parser.parse_args()

bed = {}
if args.bed:
    for line in open(args.bed):
        fields = line.strip().split()
        bed[fields[0]] = (int(fields[1]), int(fields[2]))

bam_fh = pysam.AlignmentFile(args.bam)
fasta_fh = pysam.FastaFile(args.fasta)
ref_seq = fasta_fh.fetch(reference=args.contig)

snps = [int(x)-1 for x in open(args.snps).readlines()] # snp positions are 1-pos, pysam is 0-pos

oob = 0
reads = bam_fh.fetch(args.contig)
for r in reads:
    aligned_tuples = r.get_aligned_pairs()
    seq = []

    if r.query_alignment_length <= 0.8 * len(ref_seq):
        # Drop CCS read if its smaller than 80% of the length of the full haplotype
        continue

    # Check if the read covers the area of interest
    if bed[args.contig]:
        if r.reference_start > bed[args.contig][0]-1:
            oob += 1
            continue
        elif r.reference_end < bed[args.contig][1]-1:
            oob += 1
            continue

    first = True
    for at in aligned_tuples:
        qry = at[0] # position along CCS read
        ref = at[1] # position aginst reference haplotype

        if ref in snps:
            if first and args.padstart:
                seq.extend(['_'] * len([snp for snp in snps if snp < ref]))
                first = False

            if not qry:
                # Deletion in ZMW read
                #seq.append(ref_seq[ref])
                seq.append(args.delchar)
                continue
            # else use the SNP as seen on the CCS query_sequence
            seq.append(r.query_sequence[qry])
        else:
            if not args.onlysnps:
                if not ref:
                    # Insertion against reference, discard it
                    continue
                if not qry:
                    # Deletion in ZMW
                    seq.append(args.delchar)
                    continue

                # else this site is not an original Gretel SNP so fill it in with the reference
                # recall ref is pysam 0-pos
                if args.useref:
                    seq.append(ref_seq[ref])
                else:
                    seq.append(r.query_sequence[qry])

    print(">%s/compressed\n%s" % (r.query_name, "".join(seq)))
sys.stderr.write("[INFO] %d ZMW dropped as OoB\n" % oob)
