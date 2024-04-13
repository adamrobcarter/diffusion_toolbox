import common
import MSD.MSD
import numpy as np

for file in common.files_from_argv('particle_linking/data/', 'trajs_'):
    data = common.load(f'particle_linking/data/trajs_{file}.npz')
    particles = data['particles']
    msd, msd_unc = MSD.MSD.calc(particles)

    np.savez(f'MSD/data/msd_{file}', msd=msd, msd_unc=msd_unc, time_step=data['time_step'])