import json
import os
import sys

import pysam
import scipy.spatial.distance as scipy_dist
import numpy as np

RESULT_DIR = sys.argv[1]
UUID_MAP_F = sys.argv[2]

ROOT_DIR = sys.argv[3]

uuids = {}
for line in open(UUID_MAP_F):
    fields = line.strip().split("\t")
    uuids[fields[0]] = {
        "read_len": int(fields[1]),
        "coverage": int(fields[2]),
        "run": int(fields[3]),
        "tree_i": int(fields[5]),
        "diversity": float(fields[8]),
        "tree_uuid": fields[6],
        "haps_uuid": fields[7],
    }

print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
    "hapn", "snpn", "run", "hapi", "recovery", "master_dist", "best_it", "cov", "readsize", "uuid", "dropped", "snptot", "propdrop", "tree",
))

for out_fn in [f for f in os.listdir(RESULT_DIR) if f.endswith("out.fasta")]:
    out_uuid = out_fn.replace(".out.fasta", "")

    master_fn = ROOT_DIR + "/masters/" + uuids[out_uuid]["tree_uuid"] + ".master.fa"
    master_h = pysam.FastaFile(master_fn)
    master = master_h.fetch(reference="1")

    in_haps_fn = ROOT_DIR + "/haps/" + uuids[out_uuid]["haps_uuid"] + ".fa"
    in_haps = pysam.FastaFile(in_haps_fn)

    NO_OUT = False

    out_haps_fn = RESULT_DIR + "/" + out_fn
    out_haps = None
    if not os.path.exists(out_haps_fn) or os.path.getsize(out_haps_fn) == 0:
        NO_OUT = True
    else:
        out_haps = pysam.FastaFile(out_haps_fn)


    try:
        vcf_fn = ROOT_DIR + "/haps/" + uuids[out_uuid]["haps_uuid"] + ".fa.vcf"
        #vcf_h = open(RESULT_DIR + "/" + out_fn.replace("out.fasta", "raw.vcf"))
        vcf_h = open(vcf_fn)
        snps = [int(fields[1])-1 for fields in [line.split("\t") for line in vcf_h.readlines() if line[0] != "#"]]
        vcf_h.close()
    except IOError as e:
        vcf_h.close()
        sys.stderr.write(str(e))
        import sys; sys.exit(3);

    best_scores = {}
    diffs = []

    for in_ref in sorted(in_haps.references):
        best_scores[in_ref] = [1.0, "-1__-1", 0]
        if NO_OUT:
            continue

        for out_ref in sorted(out_haps.references):
            in_seq = in_haps.fetch(reference=in_ref)
            out_seq = out_haps.fetch(reference=out_ref)

            try:
                np_inseq = np.array(list(in_seq))[snps]
                np_outseq = np.array(list(out_seq))[snps]
                distance = scipy_dist.hamming(np_inseq, np_outseq)

            except ValueError:
                print len(in_seq), len(out_seq), metahaplome, run
                import sys; sys.exit(2)

            if distance < best_scores[in_ref][0]:
                ref_dist = scipy_dist.hamming(list(in_seq), list(out_seq))
                best_scores[in_ref] = [distance, out_ref, ref_dist]
                diffs = [i for i in range(len(in_seq)) if in_seq[i]!=out_seq[i]]

    for in_ref in sorted(best_scores):
        in_ref_i = ord(in_ref) - 65
        print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
            5,
            uuids[out_uuid]["diversity"],
            uuids[out_uuid]["run"],
            in_ref,
            (1.0 - best_scores[in_ref][0])*100,
            (best_scores[in_ref][2])*100,
            int(best_scores[in_ref][1].split("__")[0])+1,
            uuids[out_uuid]["coverage"],
            uuids[out_uuid]["read_len"],
            out_uuid,
            0,
            len(snps),
            0.0,
            uuids[out_uuid]["tree_i"],
        ))

    in_haps.close()
    master_h.close()

    if not NO_OUT:
        out_haps.close()

