import numpy as np
import time

t0 = time.time()
d = np.load('particle_detection/data/particles_eleanorlong.npz')
data = d['particles']
print(data.shape)
t1 = time.time()

t2 = time.time()
data2 = np.load('particle_detection/data/particles_eleanorlong.npy')
print(data2.shape)
t3 = time.time()

print(f'npz:{t1-t0:.1f}s, npy:{t3-t2:.1f}s')