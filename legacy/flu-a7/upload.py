import chitin
import sys

chitin.util.copy_experiment_archive(sys.argv[1], "bert", new_root="/ibers/ernie/scratch/msn/flu-a7/", manifest=True)
