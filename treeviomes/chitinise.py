import chitin
import uuid
import os

exp_uuid = str(uuid.uuid4())
try:
    os.mkdir(exp_uuid)
except OSError:
    import sys; sys.exit(1)

c = chitin.Chitin()
c.skip_integrity = True
c.ignore_dot=True
c.ignore_parents=True

HAPS = [5]

ofh = open(exp_uuid + "/"+ exp_uuid + ".manifest.chitin", "w")
for tree_i in range(5):
    for divergence in [0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 5, 10, 25]:
        divergence = divergence/100.0

        tree_uuid = str(uuid.uuid4())
        command_set = c.parse_script("scripts/make_tree_seqs.sh", None, max(HAPS), divergence, tree_uuid, exp_uuid)
        cr = c.super_handle(command_set)

        for n_hap in HAPS:
            haps_uuid = str(uuid.uuid4())
            command_set = c.parse_script("scripts/make_haps.sh", n_hap, tree_uuid, haps_uuid, exp_uuid)
            c.super_handle(command_set)

            for coverage in [3, 5, 7, 10, 25, 50, 75, 100, 250, 500]:
                for read_size in [50, 75, 100, 150, 250, 500]:
                    # Sample reads 10 times and call gretel
                    for i in range(10):
                        run_uuid = str(uuid.uuid4())
                        ofh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                            run_uuid,
                            read_size,
                            coverage,
                            i,
                            n_hap,
                            tree_i,
                            tree_uuid,
                            haps_uuid,
                            divergence,
                        ))


                        #command_set = c.parse_script("scripts/make_reads.sh", read_size, run_uuid, haps_uuid, coverage, exp_uuid)
                        #c.super_handle(command_set)
