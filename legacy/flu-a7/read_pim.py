LIMIT = 99.0

for line in open("pim"):
    if line[0] == "#" or len(line.strip()) == 0:
        continue

    fields = line.split()
    i = int(fields[0][:-1])-1
    name = fields[1]
    lower_diag_r = [float(x) for x in fields[2 : 2+i]]
    for j, identity in enumerate(lower_diag_r):
        if identity >= LIMIT:
            print line.split()[1], i, j, min(lower_diag_r)
            break

