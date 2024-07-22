import numpy as np

arr = np.random.rand(338, 423, 71998)
print(arr.nbytes/1e9)
# takes 78.7GB
print(arr.mean())