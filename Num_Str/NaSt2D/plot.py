import numpy as np

all = np.fromfile("source/wave.out", dtype=float, count=1)
print(all)
