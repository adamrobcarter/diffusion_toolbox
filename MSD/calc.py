import common
import MSD.MSD
import numpy as np
import time

# no parrelelisation currently

def go(file):
    data = common.load(f'particle_linking/data/trajs_{file}.npz')
    particles = data['particles']

    t0 = time.time()

    if 'mixt' in file:
        # max_time_origins = 4000

        # t = [0, 0.5, 15.5, 16, 16.5, 32, 64, 128, 256, 512, 1024]
        t = [0, 0.5]

        msd, msd_unc = MSD.MSD.calc_mixt_new(particles, t)
        # msd, msd_unc = MSD.MSD.calc_mixt_new(particles, [0, 0.5, 16, 32, 64])
        # msd, msd_unc = MSD.MSD.calc_mixt(particles, [0, 1, 16, 512], max_time_origins)

        assert msd[0] == 0, f'msd[0] = {msd[0]}'

    elif particles.size > 1e6 or True:
        # print('consider "or True" you dummy')
        print('Calculating incrementally')
        # we need the data sorted by ID (col 4) for this
        msd, msd_unc = MSD.MSD.calc_incremental(particles)

        t = np.arange(0, msd.size) * data['time_step']

    else:
        msd, msd_unc = MSD.MSD.calc(particles)
    t1 = time.time()

    print('msd', msd[0], msd[1])
    print(msd[1]/(4*data['time_step']))

    assert msd[1] > 0

    assert common.nanfrac(msd[1:]) == 0, f'nanfrac(msd[1:]) = {common.nanfrac(msd[1:])}'

    common.save_data(f'MSD/data/msd_{file}', msd=msd, msd_unc=msd_unc, t=t,
        particle_diameter=data.get('particle_diameter'), pack_frac_given=data.get('pack_frac_given'), pack_frac=data.get('pack_frac'),
        pixel_size=data.get('pixel_size'), window_size_x=data.get('window_size_x'), window_size_y=data.get('window_size_y'),
        computation_time=t1-t0
    )

if __name__ == '__main__':
    for file in common.files_from_argv('particle_linking/data/', 'trajs_'):
        go(file)