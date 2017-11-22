import chitin
import os

import sys
fasta = os.path.abspath(sys.argv[1])

c = chitin.Chitin(MAX_PROC=4)
c.skip_integrity = True
c.ignore_dot=True
c.ignore_parents=True
c.show_stderr=True

exp = chitin.util.register_experiment(os.path.abspath("data"), create_dir=True)

for read_size in [50, 100, 150, 250]:
    for cov in [3, 5, 7, 10, 25, 50]:
        for i in range(100):

            run = chitin.util.register_run(exp.uuid, meta={
                "read_set": i,
                "read_size": read_size,
                "read_coverage": cov,
            }, create_dir=True)
            commands = c.parse_script("PREPARE.sh", read_size, run.get_path(), cov, fasta)
            c.execute(commands, run=run.uuid)

