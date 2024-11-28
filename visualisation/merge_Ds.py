import common
import numpy as np

datas = []

SOURCE = 'f_first_first'
CROSSOVER = 0.6

files = common.files_from_argv('box_counting/data/', 'counted_')
assert len(files) == 2

data0 = common.load(f'visualisation/data/Ds_from_{SOURCE}_{files[1]}.npz')
data1 = common.load(f'visualisation/data/Ds_from_{SOURCE}_{files[0]}.npz')
k0 = data0['ks']
k1 = data1['ks']

assert k0.min() < CROSSOVER
assert k0.max() > CROSSOVER
assert k1.min() < CROSSOVER
assert k1.max() > CROSSOVER


assert data0['max_time_hours'] > data1['max_time_hours']
assert (particle_diameter := data0['particle_diameter']) == data1['particle_diameter']

data0_end_index   = np.argmax(k0 >  CROSSOVER)
data1_start_index = np.argmax(k1 >= CROSSOVER)
print(k0)
# print
print(data0_end_index, data1_start_index)
assert data1_start_index > 0
assert data0_end_index   > 0
assert data1_start_index < k1.size
assert data0_end_index   < k0.size

def concat(arr0, arr1):
    print(data0_end_index, data1_start_index)
    d0 = arr0[:data0_end_index]
    d1 = arr1[data1_start_index:]
    print(d0.shape, d1.shape)
    # ratio_at_join
    # print(d0)
    # print(d0, d1)
    return np.concatenate((d0, d1), axis=0)

ks     = concat(k0,              k1)
Ds     = concat(data0['Ds'],     data1['Ds'])
D_uncs = concat(data0['D_uncs'], data1['D_uncs'])

assert ks.shape == Ds.shape
assert Ds.shape == D_uncs.shape

dataout = dict(
    Ds                = Ds,
    D_uncs            = D_uncs,
    ks                = ks,
    particle_diameter = particle_diameter,
    max_time_hours    = data1['max_time_hours']
)
common.save_data(f'visualisation/data/Ds_from_{SOURCE}_{files[1]}_merged.npz', **dataout)