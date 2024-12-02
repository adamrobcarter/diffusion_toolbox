import common
import numpy as np

datas = []

files = common.files_from_argv('box_counting/data/', 'counted_')
assert len(files) == 2

data0 = common.load(f'box_counting/data/counted_{files[0]}.npz')
data1 = common.load(f'box_counting/data/counted_{files[1]}.npz')

assert data0['time_step'] < data1['time_step']

assert np.all(data0['box_sizes'] == data1['box_sizes'])
assert data0['pack_frac_given'] == data1['pack_frac_given']

t0 = data0['time_step'] * np.arange(data0['N2_mean'].shape[1])
t1 = data1['time_step'] * np.arange(data1['N2_mean'].shape[1])

if data0['time_step'] == 0.5 and data1['time_step'] == 8:
    CROSSOVER_TIME = 304
if data0['time_step'] == 0.5 and data1['time_step'] == 16:
    CROSSOVER_TIME = 992

# we need to make sure we choose a time for the crossover that's in both datasets
assert CROSSOVER_TIME in t0
assert CROSSOVER_TIME in t1
data0_end_index   = np.argmax(t0 == CROSSOVER_TIME)
data1_start_index = np.argmax(t1 == CROSSOVER_TIME)
assert data1_start_index > 0
assert data0_end_index   > 0
assert t0[data0_end_index] == t1[data1_start_index]

def concat(arr0, arr1):
    d0 = arr0[:, :data0_end_index]
    d1 = arr1[:, data1_start_index:]
    # ratio_at_join
    # print(d0)
    return np.concatenate((d0, d1), axis=1)

t       = concat(t0[np.newaxis, :], t1[np.newaxis, :]).squeeze() # the newaxis is cause concat wants a 2D array, the squeeze reverses it
N2_mean = concat(data0['N2_mean'],  data1['N2_mean'] )
N2_std  = concat(data0['N2_std'],   data1['N2_std']  )
# N_var   = concat(data0['N_var'][np.newaxis, :],   data1['N_var'][np.newaxis, :]  ).squeeze() # the newaxis is cause concat wants a 2D array, the squeeze reverses it
assert t[1] == data0['time_step']
assert np.unique(t[1:]-t[:-1])[0] == data0['time_step']
assert np.unique(t[1:]-t[:-1])[1] == data1['time_step']

def mean_and_ratio(a0, a1):
    ratio = (a0 / a1).mean()
    avg   = (a0 + a1) / 2
    return ratio, avg

N_mean_ratio, N_mean_avg = mean_and_ratio(data1['N_mean'], data0['N_mean'])
assert 0.99 < N_mean_ratio < 1.01, f'N_mean_ratio = {N_mean_ratio}'

# N_var_ratio, N_var_avg = mean_and_ratio(data1['N_var'], data0['N_var'])
# print('N var ratio', N_var_ratio)
# assert 0.99 < N_var_ratio < 1.05, f'N_var_ratio = {N_var_ratio}'

N_var_mod_ratio, N_var_mod_avg = mean_and_ratio(data1['N_var_mod'], data0['N_var_mod'])
assert 0.99 < N_var_mod_ratio < 1.26, f'N_var_mod_ratio = {N_var_mod_ratio}'

# N_var_mod_ratio, N_var_mod_avg = mean_and_ratio(data1['N_var_mod'], data0['N_var_mod'])
# assert 0.99 < N_var_mod_ratio < 1.01, f'N_var_mod_ratio = {N_var_mod_ratio}'

##############################################
L_CROSSOVER = 3
L_long  = data1 ['box_sizes']
L_short = data0['box_sizes']

assert data1['max_time_hours'] > data0['max_time_hours']
assert (particle_diameter := data1['particle_diameter']) == data0['particle_diameter']

data1_start_index = np.argmax(L_long  >  L_CROSSOVER)
data0_end_index  = np.argmax(L_short >= L_CROSSOVER)
# print(k_long)
# print
# print(data1_end_index, data0_start_index)
assert data0_end_index  > 0
assert data1_start_index > 0
assert data0_end_index  < L_short.size
assert data1_start_index < L_long.size

def concat_L(arr_long, arr_short):
    if len(arr_short.shape) == 2:
        d_short = arr_short[:, :data0_end_index]
        d_long  = arr_long [:, data1_start_index:]
        return np.concatenate((d_short, d_long), axis=1)
    else:
        d_short = arr_short[:data0_end_index]
        d_long  = arr_long [data1_start_index:]
        return np.concatenate((d_short, d_long), axis=0)

N_var = concat_L(data1['N_var'], data0['N_var'])
print(data0['N_var'])
print(data1['N_var'])
print(N_var)

dataout = dict(
    t                 = t,
    # fitting_points    = points,
    N2_mean           = N2_mean,
    N2_std            = N2_std,
    box_sizes         = data0['box_sizes'],
    sep_sizes         = data0['sep_sizes'],
    num_boxes         = data0['num_boxes'],
    N_mean            = N_mean_avg,
    N_var             = N_var,
    N_var_mod         = N_var_mod_avg,
    N_var_mod_std     = None,
    time_step         = data0['time_step'], # we use data0 because time_step is used for the first point fit
    pack_frac         = data0.get('pack_frac'),
    particle_diameter = data0.get('particle_diameter'),
    computation_time  = None,
    box_coords        = None,
    pack_frac_given   = data0.get('pack_frac_given'),
    window_size_x     = data0['window_size_x'],
    window_size_y     = data0['window_size_y'],
    max_time_hours    = data1['max_time_hours']
)
common.save_data(f'box_counting/data/counted_{files[1]}_merged.npz', **dataout)