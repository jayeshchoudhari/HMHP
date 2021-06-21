import numpy as np
import sys

f =open("sampleSetIds.txt", "w")

setIds = np.random.choice(range(int(sys.argv[1])), int(0.01*int(sys.argv[1])))

st = '\n'.join([str(i) for i in setIds])
f.write(st)

f.close()

