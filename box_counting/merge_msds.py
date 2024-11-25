import common
import numpy as np

datas = []

files = common.files_from_argv('box_counting/data/', 'counted_')
assert len(files) == 2

data0 = common.load(f'box_counting/data/counted_{files[0]}.npz')
data1 = common.load(f'box_counting/data/counted_{files[1]}.npz')

assert data0['time_step'] < data1['time_step']

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
assert t[1] == data0['time_step']
assert np.unique(t[1:]-t[:-1])[0] == data0['time_step']
assert np.unique(t[1:]-t[:-1])[1] == data1['time_step']

# we now need to calculate fitting points
# as the points are now weirdly spaced
# we work out the fraction of points that should be taken from each dataset
target_num_points = 1000
d0_span = np.log10(CROSSOVER_TIME / t[1])
d1_span = np.log10(t[-1] / CROSSOVER_TIME)
print('spans', d0_span, d1_span, 'ratio', d0_span/d1_span)
num_points_0 = int(target_num_points * d0_span / (d0_span + d1_span))
num_points_1 = int(target_num_points * d1_span / (d0_span + d1_span))
print('nums', num_points_0, num_points_1, 'ratio', num_points_0/num_points_1)
points_0 = common.exponential_integers(1, data0_end_index,        num=num_points_0) - 1
points_1 = common.exponential_integers(data1_start_index, t.size, num=num_points_1) - 1
points = np.concatenate((points_0, points_1))
t[points] # check the slice works
print('ratios', points_0.size/points_1.size, num_points_0/num_points_1)

def mean_and_ratio(a0, a1):
    ratio = (a0 / a1).mean()
    avg   = (a0 + a1) / 2
    return ratio, avg

N_mean_ratio, N_mean_avg = mean_and_ratio(data1['N_mean'], data0['N_mean'])
assert 0.99 < N_mean_ratio < 1.01, f'N_mean_ratio = {N_mean_ratio}'

N_var_ratio, N_var_avg = mean_and_ratio(data1['N_var'], data0['N_var'])
assert 0.99 < N_var_ratio < 1.04, f'N_var_ratio = {N_var_ratio}'

N_var_mod_ratio, N_var_mod_avg = mean_and_ratio(data1['N_var_mod'], data0['N_var_mod'])
assert 0.99 < N_var_mod_ratio < 1.2, f'N_var_mod_ratio = {N_var_mod_ratio}'

# N_var_mod_ratio, N_var_mod_avg = mean_and_ratio(data1['N_var_mod'], data0['N_var_mod'])
# assert 0.99 < N_var_mod_ratio < 1.01, f'N_var_mod_ratio = {N_var_mod_ratio}'

dataout = dict(
    t                 = t,
    # fitting_points    = points,
    N2_mean           = N2_mean,
    N2_std            = N2_std,
    box_sizes         = data0['box_sizes'],
    sep_sizes         = data0['sep_sizes'],
    num_boxes         = data0['num_boxes'],
    N_mean            = N_mean_avg,
    N_var             = N_var_avg,
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