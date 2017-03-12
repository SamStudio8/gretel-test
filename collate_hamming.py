"""Given:
    * A result directory
    * A CLUSTAL alignment of the known haplotypes
    * A VCF generated from the CLUSTAL alignment
    * A BED containing the known haplotypes names and start/ends
    * A data directory,
Produce, a table containing a recovery rate for each of the haplotypes in the BED file,
for each experiment in the data directory, given the results in the result directory."""
import os
import sys

import pysam

import scipy.spatial.distance as scipy_dist
import numpy as np

RESULT_DIR = sys.argv[1]
ALIGNED_FAN = sys.argv[2]

GOLD_VCFN = None
try:
    GOLD_VCFN = sys.argv[3]
except:
    pass

GENE_BEDN = None
try:
    GENE_BEDN = sys.argv[4]
except:
    pass

DATA_DIR = sys.argv[5]

# Read input haplotypes as aligned by clustal (or similar) and output to FASTA
# NOTE We expect sequences of equal length, deletions should be padded with -
in_haps_fh = pysam.FastaFile(ALIGNED_FAN)
in_haps = in_haps_fh.references

# Read gold reference VCF
try:
    vcf_h = open(GOLD_VCFN)
    ALL_SNPS = np.array([int(fields[1])-1 for fields in [line.split("\t") for line in vcf_h.readlines() if line[0] != "#"]])
    vcf_h.close()
except:
    pass

# Read BED
haps_bed = {}
try:
    BED = open(GENE_BEDN)
    for line in BED:
        fields = line.strip().split()
        haps_bed[fields[3]] = {
            "start": int(fields[1]),
            "end": int(fields[2]),
        }
except:
    pass

# Read crumbs
crumb_ranks = {}
try:
    for f_crumb in os.listdir(RESULT_DIR):
        if not f_crumb.endswith(".crumbs"):
            continue
        uuid = f_crumb.strip().replace(".crumbs", "")
        current_crumbs_uw = []
        current_crumbs_w = []
        crumb_ranks[uuid] = {}
        with open(RESULT_DIR + "/" + f_crumb) as fh_crumb:
            for line in fh_crumb:
                if line[0]=="#":
                    continue
                fields = line.strip().split("\t")
                current_crumbs_uw.append(abs(float(fields[2])))
                current_crumbs_w.append(abs(float(fields[1])))
        crumb_ranks[uuid]["uw"] = np.argsort(np.argsort(current_crumbs_uw))
        crumb_ranks[uuid]["w"] = np.argsort(np.argsort(current_crumbs_w))
except Exception as e:
    print e
    import sys; sys.exit(4)

print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
    "uuid",
    "in_hname",
    "out_hname",
    "recovery",
    "snptot",
    "rank_w",
    "rank_uw",
    "reads_propdrop",
    "n_outhaps",
))

rank = open("ranks", "w")
rank.write("uuid\tin_hname\torder\trank\n")
for out_fn in [f for f in os.listdir(RESULT_DIR) if f.endswith("out.fasta")]:
    uuid = out_fn.replace(".out.fasta", "")
    sys.stderr.write("%s\n" % uuid)

    uuid_ranks = []

    # Get the output haplotypes
    out_haps_fn = RESULT_DIR + "/" + out_fn
    if not os.path.exists(out_haps_fn) or os.path.getsize(out_haps_fn) == 0:
        # No result
        out_haps_fh = None
        out_haps = []
    else:
        out_haps_fh = pysam.FastaFile(out_haps_fn)
        out_haps = out_haps_fh.references

    uuid_vcf_fh = open(RESULT_DIR + "/" + uuid + ".raw.vcf")
    SUB_ALL_SNPS = np.array([int(fields[1])-1 for fields in [line.split("\t") for line in uuid_vcf_fh.readlines() if line[0] != "#"]])
    uuid_vcf_fh.close()

    read_stats = {}
    for in_hname in in_haps:
        read_stats[in_hname] = {"dropped": 0, "not_dropped":0, "total": 0}

    if os.path.exists(DATA_DIR + "/" + uuid + "/unaligned.fq"):
        fq = open(DATA_DIR + "/" + uuid + "/unaligned.fq")
        for line in fq:
            if line.startswith("@HOOT"):
                g = (line.strip().split("_%")[0].split("_", 1)[1])
                if g in read_stats:
                    read_stats[g]["dropped"] += 1
                    read_stats[g]["total"] += 1
                else:
                    print g
                    sys.exit(5);
        fq.close()

    if os.path.exists(DATA_DIR + "/" + uuid + "/out.sam"):
        sam = open(DATA_DIR + "/" + uuid + "/out.sam")
        for line in sam:
            if line.startswith("HOOT"):
                g = (line.strip().split("_%")[0].split("_", 1)[1])
                if g in read_stats:
                    read_stats[g]["not_dropped"] += 1
                    read_stats[g]["total"] += 1
                else:
                    sys.exit(5);
        sam.close()

    best_score_d = {}
    for in_hname in in_haps:
        best_score_d[in_hname] = {
            "hamming": 1.0,
            "out_hname": None,
            "hamming_prop": 0.0,
            "snp_n": 0,
            "rank_w_0": -1,
            "rank_uw_0": -1,
        }
        in_seq = in_haps_fh.fetch(reference=in_hname)

        if in_hname in haps_bed:
            start = haps_bed[in_hname]["start"]
            end = haps_bed[in_hname]["end"]
        else:
            start = 0
            end = len(in_seq)

        for out_hname in out_haps:
            out_seq = out_haps_fh.fetch(reference=out_hname)

            try:
                indices = np.take(ALL_SNPS, np.where( (ALL_SNPS >= start) & (ALL_SNPS < end) ))[0]
                np_inseq = np.array(list(in_seq))[indices]
                np_outseq = np.array(list(out_seq))[indices]
                distance = scipy_dist.hamming(np_inseq, np_outseq)
            except ValueError:
                print uuid, in_hname, out_hname
                sys.exit(2)

            if distance < best_score_d[in_hname]["hamming"]:
                best_score_d[in_hname]["hamming"] = distance
                best_score_d[in_hname]["hamming_prop"] = (1.0 - distance) * 100
                best_score_d[in_hname]["out_hname"] = out_hname
                best_score_d[in_hname]["snp_n"] = len(SUB_ALL_SNPS)
                best_score_d[in_hname]["rank_uw_0"] = crumb_ranks[uuid]["uw"][int(out_hname.split("__")[0])]
                best_score_d[in_hname]["rank_w_0"] = crumb_ranks[uuid]["w"][int(out_hname.split("__")[0])]

    for in_hname in sorted(in_haps):
        print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                uuid,
                in_hname.split(".")[0], # Drop accession version
                best_score_d[in_hname]["out_hname"],
                best_score_d[in_hname]["hamming_prop"],
                best_score_d[in_hname]["snp_n"],
                best_score_d[in_hname]["rank_w_0"],
                best_score_d[in_hname]["rank_uw_0"],
                float(read_stats[in_hname]["dropped"]) / read_stats[in_hname]["total"] * 100,
                len(out_haps),
        ))
        uuid_ranks.append([in_hname, best_score_d[in_hname]["rank_w_0"]])

    for hap_order, hap_meta in enumerate(sorted(uuid_ranks, key=lambda x: x[1])):
        if hap_meta[1] >= 0:
            rank.write("%s\t%s\t%d\t%d\n" % (uuid, hap_meta[0], hap_order+1, hap_meta[1]))
rank.close()

