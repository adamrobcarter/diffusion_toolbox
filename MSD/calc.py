import common
import MSD.MSD
import numpy as np
import time

# no parrelelisation currently

for file in common.files_from_argv('particle_linking/data/', 'trajs_'):
    data = common.load(f'particle_linking/data/trajs_{file}.npz')
    particles = data['particles']

    t0 = time.time()
    if np.prod(particles.shape) > 1e7 or True:
        print('consider "or True" you dummy')
        print('Calculating incrementally')
        # we need the data sorted by ID (col 4) for this
        msd, msd_unc = MSD.MSD.calc_incremental(particles)
    else:
        msd, msd_unc = MSD.MSD.calc(particles)
    t1 = time.time()

    print('msd', msd[0], msd[1])
    print(msd[1]/(4*data['time_step']))

    common.save_data(f'MSD/data/msd_{file}', msd=msd, msd_unc=msd_unc, time_step=data['time_step'],
                     computation_time=t1-t0)