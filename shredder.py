import argparse
import sys
import random
import numpy as np

import pysam

parser = argparse.ArgumentParser(description="A very bad read generator.")
parser.add_argument("shred", type=int)
parser.add_argument("fasta")
parser.add_argument("--hits")

parser.add_argument("--subject", default=False, action='store_true')
parser.add_argument("--cover", default=3, type=int)
parser.add_argument("--error", default=0.0, type=float, help="percentage synthetic error")
parser.add_argument("--fa", default=False, action="store_true")
parser.add_argument("--sam", default=None)

parser.add_argument("--bed", default=None)

ARGS = parser.parse_args()
fasta = pysam.FastaFile(ARGS.fasta)

if ARGS.fa:
    readsym = ">"
else:
    readsym = "@HOOT_"


beds = {}
if ARGS.bed:
    beds = {f[0]: (int(f[1]), int(f[2])) for line in open(ARGS.bed).readlines() for f in [line.strip().split()]}

hits = {}
hit_file = None
if ARGS.hits:
    hit_file = open(ARGS.hits)

    for line in hit_file:
        fields = line.strip().split("\t")
        hits[fields[1]] = ( int(fields[8])-1, int(fields[9]) )

sam_h = [
    "@HD	VN:1.0	SO:unsorted",
    "@SQ	SN:HOOT	LN:%d" % max(fasta.lengths)
]
sam = []
def output_read(curr_seq, curr_name, start_1pos):
    print "%s%s" % (readsym, curr_name)
    curr_seq = "".join(curr_seq)
    print curr_seq
    if not ARGS.fa:
        print '+'
        print 'J' * len(curr_seq)
    sam.append([
        curr_name, 0, "HOOT", start_1pos, 42, str(len(curr_seq))+'M', '*', 0, len(curr_seq), curr_seq, "J"*(len(curr_seq)), "AS:i:0", "XN:i:0", "XM:i:0", "XO:i:0", "XG:i:0", "NM:i:0", "MD:Z:"+str(len(curr_seq)), "YT:Z:UU", "RG:Z:hoot"
    ])


for reference in fasta.references:
    if hit_file:
        try:
            seq = fasta.fetch(reference=fields[1], start=int(fields[8])-1, end=int(fields[9]))
        except KeyError:
            sys.stderr.write("Reference '%s' not in FASTA.\n" % fields[1])
            continue
        ref_extract = fields[1].split("|")[-2] #TODO ew
    else:
        seq = fasta.fetch(reference=reference, start=0, end=fasta.get_reference_length(reference))
        try:
            ref_extract = reference.split("|")[-2] #TODO ew
        except IndexError:
            ref_extract = reference

    if ARGS.shred == 0:
        ARGS.shred = len(shred)


    SEQ_START = 0
    SEQ_END = len(seq)
    SEQ_LEN = len(seq)
    if reference in beds:
        SEQ_START = beds[reference][0]
        SEQ_END = beds[reference][1]

        if SEQ_START > SEQ_END:
            SEQ_START, SEQ_END = SEQ_END, SEQ_START

        SEQ_LEN = len(seq[SEQ_START:SEQ_END])

    n_reads = int((ARGS.cover * (SEQ_LEN+ARGS.shred)) / ARGS.shred)

    for read_i, read_start in enumerate(np.random.randint( SEQ_START-(ARGS.shred/2), SEQ_END-(ARGS.shred/2), size=n_reads)):

        read_len = ARGS.shred
        if read_start < 0:
            read_len = ARGS.shred - abs(read_start)
            read_start = 0

        read_headbase = "%s_%%_%d" % (ref_extract, read_i)
        curr_name = "%s_%d" % (read_headbase, read_start)
        curr_seq = seq[read_start : read_start + read_len]

        if ARGS.error > 0:
            curr_seq = list(curr_seq)
            for ibase, base in enumerate(curr_seq):
                if random.uniform(0,100) < ARGS.error:
                    curr_seq[ibase] = random.choice(list(set(["A", "C", "G", "T"]) - set([base])))
        output_read(curr_seq, curr_name, read_start+1)



if ARGS.sam:
    sam_f = open(ARGS.sam, "w+")
    sam_h.extend(["@RG	ID:hoot	PL:Illumina	PU:Illumina	SM:hoot",
    "@PG	ID:shredder	PN:shredder	VN:1	CL:\"shredder4.py\""])
    sam_f.write("\n".join(sam_h))
    for r in sam:
        sam_f.write("\n"+ "\t".join([str(s) for s in r]))
    sam_f.close()
