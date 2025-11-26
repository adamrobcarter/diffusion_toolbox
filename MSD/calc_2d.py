import common
import MSD.MSD
import numpy as np
import time
import MSD.show

# no parrelelisation currently
# calculates only over x and y

# there is a lot of repeated code in here and calc.py :(

def go(file):
    data = common.load(f'particle_linking/data/trajs_{file}.npz')
    particles = data['particles']
    
    
    time_column = data.get('dimension', 2)
    frames = np.unique(particles[:, time_column])
    frame_interval = np.diff(frames)
    consistent_timesteps = np.all(np.isclose(frame_interval, frame_interval[0])) # need isclose in case time is irrational


    t0 = time.time()

    timestep = data['time_step']

    # if 'mixt' in file or 'shear' in file:
    if not consistent_timesteps:
        # max_time_origins = 4000

        # t = [0, 0.5, 15.5, 16, 16.5, 32, 64, 128, 256, 512, 1024]
        # t = [0, 0.5]
        frames = np.unique(particles[:, data.get('dimension', 2)])
        t_indices = common.exponential_indices(frames, num=120)[:80] # remove the big t ones, otherwise it's rather slow
        frames = frames[t_indices]
        frames = [0, *frames]
        assert frames[0] == 0
        print('times', np.unique(particles[:, data.get('dimension', 2)]))

        msd, msd_unc = MSD.MSD.calc_mixt_new_xyz(particles, frames, data.get('dimension', 2))
        # msd, msd_unc = MSD.MSD.calc_mixt_new(particles, [0, 0.5, 16, 32, 64])
        # msd, msd_unc = MSD.MSD.calc_mixt(particles, [0, 1, 16, 512], max_time_origins)

        assert np.all(msd[:, 0]) == 0, f'msd[0] = {msd[:, 0]}'

    else:
        
        if frame_interval[0] != 1:
            particles[:, time_column] /= frame_interval[0]
            # if the timestep was irrational these might not be integer, so we fix that
            particles[:, time_column] = np.round(particles[:, time_column])

        msd, msd_unc = MSD.MSD.calc_incremental_xyz(particles, data.get('dimension', 2))

    msd_2d = msd[0, :] + msd[1, :]
    msd_2d_unc = msd_unc[0, :] + msd_unc[1, :]

    t1 = time.time()

    assert np.all(msd[:, 1]) > 0

    assert common.nanfrac(msd[1:]) == 0, f'nanfrac(msd[1:]) = {common.nanfrac(msd[1:])}'

    t = frames * timestep
    assert t[0] == 0
    assert t[1] == timestep

    common.save_data(f'MSD/data/msd_{file}_2d', msd=msd_2d, msd_unc=msd_2d_unc, t=t,
        particle_diameter=data.get('particle_diameter'), pack_frac_given=data.get('pack_frac_given'), pack_frac=data.get('pack_frac'),
        pixel_size=data.get('pixel_size'), window_size_x=data.get('window_size_x'), window_size_y=data.get('window_size_y'),
        computation_time=t1-t0, dimension=data.get('dimension', 2),
    )

if __name__ == '__main__':
    for file in common.files_from_argv('particle_linking/data/', 'trajs_'):
        go(file)
        
        MSD.show.go(f'{file}_2d')