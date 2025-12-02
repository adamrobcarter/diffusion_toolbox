import common
import MSD.MSD
import numpy as np
import time
import MSD.show

# no parrelelisation currently

def go(file, quiet=False):
    data = common.load(f'particle_linking/data/trajs_{file}.npz', quiet=quiet)
    particles = data['particles']

    
    time_column = data.get('dimension', 2)
    assert particles.size
    t = np.unique(particles[:, time_column])
    assert len(t)
    t_interval = np.diff(t)
    consistent_timesteps = np.all(np.isclose(t_interval, t_interval[0])) # need isclose in case time is irrational

    t0 = time.time()

    # if 'mixt' in file or 'shear' in file or 'sim_porous' in file:
    if not consistent_timesteps:
        print('calculating mixt')
        # max_time_origins = 4000

        # t = [0, 0.5, 15.5, 16, 16.5, 32, 64, 128, 256, 512, 1024]
        # t = [0, 0.5]
        t_indices = common.exponential_indices(t, 120) # remove the big t ones, otherwise it's rather slow
        t = t[t_indices]
        t = [0, *t]
        assert t[0] == 0

        msd, msd_unc = MSD.MSD.calc_mixt_new(particles, t, dimension=data.get('dimension', 2))
        # msd, msd_unc = MSD.MSD.calc_mixt_new(particles, [0, 0.5, 16, 32, 64])
        # msd, msd_unc = MSD.MSD.calc_mixt(particles, [0, 1, 16, 512], max_time_origins)

        assert msd[0] == 0, f'msd[0] = {msd[0]}'

    elif particles.size > 1e6 or True:

        timestep = data['time_step']
        
        if t_interval[0] != 1:
            print('converting time')
            particles[:, time_column] /= t_interval[0]
            # if the timestep was irrational these might not be integer, so we fix that
            particles[:, time_column] = np.round(particles[:, time_column])

            timestep *= t_interval[0] # wtf is this?

        # print('consider "or True" you dummy')
        if not quiet: print('Calculating incrementally')
        # we need the data sorted by ID (col 4) for this
        msd, msd_unc = MSD.MSD.calc_incremental(particles, data.get('dimension', 2))

        t = np.arange(0, msd.size) * timestep

    else:
        msd, msd_unc = MSD.MSD.calc(particles)
        assert False
    t1 = time.time()

    if not quiet: print('msd', msd[0], msd[1])
    # print(msd[1]/(4*data['time_step']))

    assert msd[1] > 0

    assert common.nanfrac(msd[1:]) == 0, f'nanfrac(msd[1:]) = {common.nanfrac(msd[1:])}'

    
    assert t[0] == 0
    assert t[1] == timestep

    common.save_data(f'MSD/data/msd_{file}', msd=msd, msd_unc=msd_unc, t=t,
        particle_diameter=data.get('particle_diameter'), pack_frac_given=data.get('pack_frac_given'), pack_frac=data.get('pack_frac'),
        pixel_size=data.get('pixel_size'), window_size_x=data.get('window_size_x'), window_size_y=data.get('window_size_y'),
        computation_time=t1-t0, T=data.get('T'),
        quiet=quiet
    )

    # MSD.show.go(file, quiet=quiet)

if __name__ == '__main__':
    for file in common.files_from_argv('particle_linking/data/', 'trajs_'):
        go(file)