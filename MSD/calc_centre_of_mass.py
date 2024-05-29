import common
import MSD.MSD
import numpy as np
import tqdm

for file in common.files_from_argv('particle_linking/data/', 'trajs_'):
    data = common.load(f'particle_linking/data/trajs_{file}.npz')
    particles = data['particles']

    msds = []
    msd_uncs = []

    data_ = MSD.MSD.reshape(particles)
    groupsizes = [1, 4, 16, 64, 256]

    msds = np.full((len(groupsizes), data_.shape[1]), np.nan)
    msd_uncs = np.full_like(msds, np.nan)

    for group_index, groupsize in enumerate(tqdm.tqdm(groupsizes)):
        msd, msd_unc = MSD.MSD.calc_centre_of_mass(data_, groupsize)
        msds[group_index, :] = msd
        msd_uncs[group_index, :] = msd_unc

    common.save_data(f'MSD/data/msd_centre_of_mass_{file}', msds=msds, msd_uncs=msd_uncs, groupsizes=groupsizes, time_step=data['time_step'])