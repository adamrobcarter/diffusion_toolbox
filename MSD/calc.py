import common
import MSD.MSD
import numpy as np
import time

# no parrelelisation currently

for file in common.files_from_argv('particle_linking/data/', 'trajs_'):
    data = common.load(f'particle_linking/data/trajs_{file}.npz')
    particles = data['particles']

    # print('width', particles[:, 1].max(), particles[:, 0].max())

    t0 = time.time()
    if particles.size > 1e6 or True:
        # print('consider "or True" you dummy')
        print('Calculating incrementally')
        # we need the data sorted by ID (col 4) for this
        msd, msd_unc = MSD.MSD.calc_incremental(particles)
    else:
        msd, msd_unc = MSD.MSD.calc(particles)
    t1 = time.time()

    print('msd', msd[0], msd[1])
    print(msd[1]/(4*data['time_step']))

    assert common.nanfrac(msd[1:]) == 0, f'nanfrac(msd[1:]) = {common.nanfrac(msd[1:])}'

    common.save_data(f'MSD/data/msd_{file}', msd=msd, msd_unc=msd_unc, time_step=data['time_step'],
        particle_diameter=data.get('particle_diameter'), pack_frac_given=data.get('pack_frac_given'),
        pixel_size=data.get('pixel_size'), window_size_x=data.get('window_size_x'), window_size_y=data.get('window_size_y'),
        computation_time=t1-t0
    )