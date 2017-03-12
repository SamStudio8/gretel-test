import sys
import numpy as np


series = {}

import matplotlib.cm as cm

f = open(sys.argv[1])
for line in f:
    fields = line.strip().split("\t")

    if fields[1] not in series:
        series[fields[1]] = {
            "x": [],
            "y": [],
        }

    likelihood = fields[0].split("__")[1]
    print "%s\t%s\t%s" % (sys.argv[2], line.strip(), likelihood)
    series[fields[1]]["x"].append(float(fields[11]))
    series[fields[1]]["y"].append(float(likelihood))

import matplotlib.pyplot as plt
colors = iter(cm.rainbow(np.linspace(0, 1, len(series))))
for s in series:
    plt.scatter(series[s]["x"],series[s]["y"], color=next(colors))
plt.show()


