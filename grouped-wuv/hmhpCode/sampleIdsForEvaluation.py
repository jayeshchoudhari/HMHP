import numpy as np

f =open("sampleSetIds.txt", "w")

setIds = np.random.choice(range(1000000), 10000)

st = '\n'.join([str(i) for i in setIds])
f.write(st)

f.close()

