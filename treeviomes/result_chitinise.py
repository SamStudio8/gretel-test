import chitin
import os
import sys

uuid_to_hap = {}
MANIFEST = open(sys.argv[2])
for line in MANIFEST:
    fields = line.strip().split("\t")
    uuid_to_hap[fields[0]] = fields[7]

c = chitin.Chitin(MAX_PROC=12)
c.skip_integrity=True
c.show_stderr=True
c.ignore_parents=True
c.ignore_dot=True

DATA = sys.argv[3]

for result_f in [f for f in os.listdir(sys.argv[1]) if f.endswith("out.fasta")]:
    result_uuid = result_f.split(".")[0]

    haps_fa = DATA + '/haps/' + uuid_to_hap[result_uuid] + '.fa'
    commands = c.parse_script("TREE.sh", os.path.abspath(sys.argv[1]), result_uuid, haps_fa)
    c.execute(commands)
