import sys
import pysam
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--bam', help='path to bam of ZMW aligned to haplotypes', required=True)
parser.add_argument('--infasta', help='reference fasta', required=True)
parser.add_argument('--outfasta', help='output "cleaned" ZMW fasta')
parser.add_argument('--snps', help='new line delimited contig and 1-pos based snps, contig names are matched greedily', required=True)
parser.add_argument('--onlysnps', action="store_true", help="only output bases on snp positions")
parser.add_argument('--onlybed', action="store_true", help="dont just mark reads with --bed, actually drop them")
#parser.add_argument('--useref', action="store_true", help="output reference bases at non-snp positions, not zmw")
parser.add_argument('--delchar', default='-', help="character to output for a deletion [default -]")
parser.add_argument('--padstart', action="store_true", help="pad the sequence with _ to maintain integrity of positions against the reference")
parser.add_argument('--bed', help="mark reads on a contig that start after or end before the provided co-ordinates [use 1-pos]")
parser.add_argument('--qual', default=40, help="quality threshold below which to ignore ZMW basecalls [default 40]", type=int)
parser.add_argument('--counts', help="path to file containing ZMW name, ZMW pass count pairs")

args = parser.parse_args()

snps = {}
for line in open(args.snps):
    fields = line.strip().split()
    if fields[0] not in snps:
        snps[fields[0]] = set([])
    snps[fields[0]].add( int(fields[1])-1 )

bed = {}
if args.bed:
    for line in open(args.bed):
        fields = line.strip().split()
        bed[fields[0]] = (int(fields[1]), int(fields[2]))

pass_counts = {}
if args.counts:
    pass_counts = { l.split()[0]: int(l.split()[1]) for l in open(args.counts).readlines() }

bam_fh = pysam.AlignmentFile(args.bam)
fasta_fh = pysam.FastaFile(args.infasta)

outfasta_fh = None
if args.outfasta:
    outfasta_fh = open(args.outfasta, 'w')

oob_dropped = 0
for ghaplotype in bam_fh.references:
    reads = bam_fh.fetch(ghaplotype)
    ghap_ref = fasta_fh.fetch(reference=ghaplotype)

    for r in reads:
        matches = mismatches = deletions = qualmask = 0
        mismatch_pos = []
        mismatch_mut = []
        mismatch_rank = []
        oob = -1

        aligned_tuples = r.get_aligned_pairs()
        seq = []

        # Pick a suitable BED/SNP list
        short = None
        for b in snps:
            if ghaplotype.startswith(b):
                if short is not None:
                    sys.stderr.write("[OHNO] Haplotype name matches more than one entry in the SNP list... Taking the first...\n")
                short = b

        if short:
            if bed[short]:
                if r.reference_start > bed[short][0]-1:
                    oob += 1
                elif r.reference_end < bed[short][1]-1:
                    oob += 2
            if oob > -1:
                oob_dropped += 1
                continue

        first = True
        rank = 0
        for at in aligned_tuples:
            qry = at[0] # position along CCS read
            ref = at[1] # position aginst reference haplotype

            if ref in snps[short]:
                if first and args.padstart:
                    pad = len([snp for snp in snps[short] if snp < ref])
                    seq.extend(['_'] * pad)
                    rank=pad
                    first = False

                if not qry:
                    # Deletion in ZMW read
                    #seq.append(ref_seq[ref])
                    seq.append(args.delchar)
                    deletions += 1
                    continue

                # else use the SNP as seen on the CCS query_sequence
                if r.query_sequence[qry] != ghap_ref[ref]:
                    if int(r.query_qualities[qry]) >= args.qual:
                        mismatches += 1
                        mismatch_pos.append(ref+1) # 1-pos for display
                        mismatch_mut.append("%s>%s" % (ghap_ref[ref], r.query_sequence[qry]))
                        mismatch_rank.append(qry+1)
                    else:
                        qualmask += 1
                else:
                    matches+=1
                seq.append(r.query_sequence[qry])
                rank += 1

            else:
                if not args.onlysnps:
                    if not ref:
                        # Insertion against reference, discard it
                        continue
                    if not qry:
                        # Deletion in ZMW
                        seq.append(args.delchar)
                        deletions += 1
                        continue
                    if r.query_qualities[qry] >= args.qual:
                        seq.append(r.query_sequence[qry])
                    else:
                        seq.append('N')

        if args.outfasta:
            outfasta_fh.write(">%s/compressed\n%s\n" % (r.query_name, "".join(seq)))
        gregion = ghaplotype.split("_")[0]
        ghaplotype_n = ghaplotype.split("_")[1]
        ghaplotype_l = ghaplotype.split("__")[1]
        print("\t".join([str(x) for x in [
            gregion,
            ghaplotype_n,
            ghaplotype_l,
            r.query_name,
            pass_counts.get(r.query_name, 0),
            oob,
            matches,
            mismatches,
            ",".join([str(x) for x in mismatch_rank]),
            ",".join([str(x) for x in mismatch_pos]),
            ",".join(mismatch_mut),
            deletions,
            qualmask,
        ]]))
sys.stderr.write("[INFO] %d ZMW dropped as OoB\n" % oob_dropped)
