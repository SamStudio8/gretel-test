import sys
from readfq import readfq # thanks heng

fa_fh = readfq(open(sys.argv[1]))
ls_fh = open(sys.argv[2])
threshold = int(sys.argv[3])

count_map = {}
for line in ls_fh:
    fields = line.strip().split('\t')
    count = fields[0]
    primary = fields[1].split(',')[0]

    if primary not in count_map:
        count_map[primary] = int(count)
    else:
        print("oh no")
        sys.exit(1)

for name, seq, qual in fa_fh:
    if count_map.get(name, 1) >= threshold:
        print('>%s_%d\n%s' % (name, count_map.get(name, 1), seq))
