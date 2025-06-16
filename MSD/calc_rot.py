import common
import MSD.MSD
import numpy as np
import time

# no parrelelisation currently

def go(file):
    data = common.load(f'particle_linking/data/trajs_{file}.npz')
    particles = data['particles']
    
    
    time_column = data.get('dimension', 2)
    t = np.unique(particles[:, time_column])
    t_interval = t[1:] - t[:-1]
    assert np.all(np.isclose(t_interval, t_interval[0])) # need isclose in case time is irrational


    t0 = time.time()

    # max_time_origins = 4000

    # t = [0, 0.5, 15.5, 16, 16.5, 32, 64, 128, 256, 512, 1024]
    # t = [0, 0.5]
    t = np.unique(particles[:, data.get('dimension', 2)])
    t_indices = common.exponential_indices(t, num=120)[:80] # remove the big t ones, otherwise it's rather slow
    t = t[t_indices]
    t = [0, *t]
    assert t[0] == 0
    print('times', np.unique(particles[:, data.get('dimension', 2)]))

    # this is hack where we move phi, theta, psi to where x, y, z were
    particles[:, [0, 1, 2]] = particles[:, [5, 6, 7]]
    particles = particles[:, [0, 1, 2, 3, 4]] # remove extra columns

    msd, msd_unc = MSD.MSD.calc_mixt_new_xyz(particles, t, data.get('dimension', 2))
    # msd, msd_unc = MSD.MSD.calc_mixt_new(particles, [0, 0.5, 16, 32, 64])
    # msd, msd_unc = MSD.MSD.calc_mixt(particles, [0, 1, 16, 512], max_time_origins)

    assert np.all(msd[:, 0]) == 0, f'msd[0] = {msd[:, 0]}'


    t1 = time.time()

    print('msd', msd[0], msd[1])
    # print(msd[1]/(4*data['time_step']))

    assert np.all(msd[:, 1]) > 0

    assert common.nanfrac(msd[1:]) == 0, f'nanfrac(msd[1:]) = {common.nanfrac(msd[1:])}'

    common.save_data(f'MSD/data/msd_rot_{file}', msd=msd, msd_unc=msd_unc, t=t,
        particle_diameter=data.get('particle_diameter'), pack_frac_given=data.get('pack_frac_given'), pack_frac=data.get('pack_frac'),
        pixel_size=data.get('pixel_size'), window_size_x=data.get('window_size_x'), window_size_y=data.get('window_size_y'),
        computation_time=t1-t0, dimension=data.get('dimension', 2),
    )

if __name__ == '__main__':
    for file in common.files_from_argv('particle_linking/data/', 'trajs_'):
        go(file)