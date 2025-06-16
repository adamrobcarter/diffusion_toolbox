import numpy as np
import common
import time, os
import tqdm
import subprocess
import itertools

def go(infile, outfile, nth_timestep=1, max_time=None,
       extra_source_file=None, skiprows=None, multiply_time_by=None, **kwargs):

    print(f'loading npz, last modified {common.get_last_modified_time(infile)} ago')
    t0 = time.time()
    data = common.load(infile)
    particles = data['particles']
    t1 = time.time()
    print(particles.shape, particles.dtype,  f'loaded {common.format_bytes(particles.nbytes)} in {t1-t0:.1f}s at {common.format_bytes(particles.nbytes/(t1-t0))}/sec')

    TIME_COLUMN = 3

    if np.isnan(particles).any():
        print('NAN FOUND')
        nans = np.isnan(particles).any(axis=1)
        particles = particles[~nans, :]
    
    print('xy min max', particles[:, 0].min(), particles[:, 0].max(), particles[:, 1].min(), particles[:, 1].max())
    
    times = np.unique(particles[:, TIME_COLUMN])
    num_timesteps = times.size
    last_timestep = times[-1]

    particles[:, TIME_COLUMN] -= particles[:, TIME_COLUMN].min()
    assert particles[:, TIME_COLUMN].min() == 0

    frame_time_deltas = particles[1:, TIME_COLUMN] - particles[:-1, TIME_COLUMN]
    assert np.all(frame_time_deltas >= 0)

    if particles.shape[1] == 9: # we have quaternion for rotation
        q_w = particles[:, 5]
        q_x = particles[:, 6]
        q_y = particles[:, 7]
        q_z = particles[:, 8]
        # convert to euler angles
        phi   = np.arctan2(2*(q_w*q_x + q_y*q_z), 1 - 2*(q_x**2 + q_y**2))
        theta = np.arcsin(2*(q_w*q_y - q_x*q_z))
        psi   = np.arctan2(2*(q_w*q_z + q_x*q_y), 1 - 2*(q_y**2 + q_z**2))
        particles[:, 5] = phi
        particles[:, 6] = theta
        particles[:, 7] = psi
        particles = particles[:, [0, 1, 2, 3, 4, 5, 6, 7]] # remove extra column

    newdata = common.copy_not_particles(data) # copy so we can modify
    newdata['particles'] = particles

    common.save_data(f'particle_linking/data/trajs_{outfile}.npz',
        **newdata,
        max_time_hours=round(last_timestep/60/60, 2),
        source_file=infile, extra_source_file=extra_source_file,
        dimension=3,
        **kwargs
    )

    # end_timestep = data[:, 2].max() // 8
    # data_small = data[data[:, 2] < end_timestep, :]
    # common.save_data(f'particle_linking/data/trajs_{outfile}_div8.npz',
    #     particles=data_small,
    #     time_step=dt, particle_diameter=particle_diameter,
    #     max_time_hours=round(end_timestep/60/60, 2),
    #     source_file=infile, extra_source_file=extra_source_file,
    #     **kwargs
    # )
        
    # end_timestep = data[:, 2].max() // 64
    # data_small = data_small[data_small[:, 2] < end_timestep, :]
    # common.save_data(f'particle_linking/data/trajs_{outfile}_div64.npz',
    #     particles=data_small,
    #     time_step=dt, particle_diameter=particle_diameter,
    #     max_time_hours=round(end_timestep/60/60, 2),
    #     source_file=infile,extra_source_file=extra_source_file,
    #     **kwargs
    # )
    # print()
    

def go_mesu_shear(filepath, suffix='', skip_rsync=False, **kwargs):
    filename = filepath.split('/')[-1]
    if not skip_rsync:
        rsync_command = ['rsync', f'cartera@login.mesu.sorbonne-universite.fr:{filepath}', f'raw_data/mesu/{filename}', '--progress']
        print('launching', ' '.join(rsync_command))
        subprocess.run(rsync_command, check=True)

    internalname = filename.split('.npz')[0]

    go(
        f'raw_data/mesu/{filename}',
        f'sim_{internalname}{suffix}',
        extra_source_file=filepath,
    )

if __name__ == '__main__':
    # go_mesu_shear('/store/cartera/shear/shear0.0_RHS.npz')
    # go_mesu_shear('/store/cartera/shear/shear0.0_theta45_RHS.npz')

    # go_mesu_shear('/store/cartera/shear/shear0.0_T0_theta45_RHS.npz')
    # go_mesu_shear('/store/cartera/shear/shear0.0_T0_theta15_RHS.npz')
    # go_mesu_shear('/store/cartera/shear/shear0.0_T0_theta30_RHS.npz')
    
    
    
    # for nblobs in [2562]:
    # # for nblobs in [42, 162, 642, 2562]:
    #     for theta in [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85]:
    #         theta = int(theta)
    #         go_mesu_shear(f'/store/cartera/shear/shear0.0_T0_theta{theta}_RHS_nblobs{nblobs}.npz')
    
    # for nblobs in [42, 162, 642, 2562]:
    #     for theta in [0, 2, 4, 6, 8, 10, 12, 14]:
    #         theta = int(theta)
    #         go_mesu_shear(f'/store/cartera/shear/shear0.0_T0_theta{theta}_RHS_nblobs{nblobs}.npz')

    # for shear in [0.0, 0.25, 0.5, 0.75, 1.0]:
    #     go_mesu_shear(f'/store/cartera/shear/shear{shear}_T0_RHS_nblobs642.npz')

    # go_mesu_shear('/store/cartera/shear/shear0.0_theta0_RHS_nblobs642.npz')
    # go_mesu_shear('/store/cartera/shear/shear0.0_theta0_RHS_nblobs42.npz')
    
    # for nblobs in [42, 162, 642, 2562]:
    #     go_mesu_shear(f'/store/cartera/shear/shear0.0_theta0_RFD_nblobs{nblobs}.npz')

        
    # for n_blobs in [42, 162]:
    #     for T in [296, 148]:
    #         go_mesu_shear(f'/store/cartera/shear/shear0.0_T{T}_theta0_RFD_nblobs{n_blobs}.npz')


        
    # t_max   = [10000]
    # t_save  = [0.01]
    # dt      = [0.01]
    # # dt      = [0.01, 0.005]
    # shear   = [0.080357]
    # # shear   = [0.0, 0.2, 0.4, 0.6, 0.8]
    # gravity = [True]
    # wall    = [True]
    # T       = [296]
    # # method = 'RFDp'
    # # method  = ['RHS', 'RFDp', 'RFDm']
    # method  = ['RHS']
    # suffix  = ['']
    # # n_blobs = [42, 162]
    # n_blobs = [162]
    # theta   = [10]
    # # theta   = [0, 2, 4, 6, 8, 10]


    # for combo in itertools.product(t_max, t_save, dt, shear, gravity, wall, theta, T, method, suffix, n_blobs):
    #     t_max, t_save, dt, shear, gravity, wall, theta, T, method, suffix, n_blobs = combo
    #     go_mesu_shear(f'/store/cartera/shear/shear{shear}_T{T}_theta{theta}_{method}_nblobs{n_blobs}_dt{dt}_tmax{t_max}.npz')


    # testing rotational diffuision
    # go_mesu_shear(f'/store/cartera/shear/shear0_T296_nograv_theta0_nowall_EMRFD_nblobs162_dt0.01_tmax1000.npz', time_step=0.01) # time_step not needed in future
    # go_mesu_shear(f'/store/cartera/shear/shear0.080357_T296_theta10_EMRFD_nblobs162_dt0.01_tmax1000.npz', time_step=0.01) # time_step not needed in future
    # go_mesu_shear(f'/store/cartera/shear/shear0_T296_theta0_EMRFD_nblobs162_dt0.01_tmax1000.npz', time_step=0.01) # time_step not needed in future

    # monolayer
    go_mesu_shear('/store/cartera/shear/shear0_T296_theta10_L50_phi0.5_Trap_nblobs42_dt0.01_tmax10.npz')
    go_mesu_shear('/store/cartera/shear/shear0_T296_theta10_L50_phi0.01_Trap_nblobs42_dt0.01_tmax10.npz')