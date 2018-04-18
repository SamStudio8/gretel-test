"""Given a manifest file A, and a file B containing additional information to melt onto A,
produce a file containing each line of B, with the associated line in A that shares
the same UUID appended"""
import os
import sys

FILE_AN = open(sys.argv[1])
FILE_BN = open(sys.argv[2])
UUID_COL_A = 0
UUID_COL_B = 0

HEADER_A = FILE_AN.readline().strip()
HEADER_B = FILE_BN.readline().strip()

uuid_to_filea_lines = {}
for line in FILE_AN:
    line = line.strip()
    fields = line.split("\t")
    uuid = fields[UUID_COL_A]

    if uuid not in uuid_to_filea_lines:
        uuid_to_filea_lines[uuid] = line
    else:
        import sys; sys.exit(1)

print ("%s\t%s" % (HEADER_B, HEADER_A))

for line in FILE_BN:
    line = line.strip()
    fields = line.split("\t")
    uuid = fields[UUID_COL_B]

    print("%s\t%s" % (
        line,
        uuid_to_filea_lines[uuid]
    ))
