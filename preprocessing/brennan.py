"""
Test docstring for preprocessing.brennan
"""

import numpy as np
import common
import time, os
import tqdm
import subprocess
import multiprocessing
import functools

# eleanorlong is
# x min max 0.6982900415999999 293.91458112
# y min max 0.7108168319999999 367.6492512
# num_timesteps 71997.0

# def read_file_skip_timesteps(file, nth_timesteps, max_frames, expected_particles_per_frame):
#     arraysize = 0
#     if max_frames:
#         progress = tqdm.tqdm(total=max_frames)
#     else:
#         progress = tqdm.tqdm()

#     # num_lines = 0

#     with open(file) as f:
#         for line in f:
#             # num_lines += 1

#             if nth_timesteps == 1:
#                 yield line

#                 if max_frames == None:
#                     progress.update()

#                 else:
#                     # we can use the density to guess the current frame number
#                     # which will make it quicker as we won't have to parse the data
#                     # actually it doesn't seem to be quicker rip
#                     # t = num_lines / expected_particles_per_frame
                    
#                     # 4.04it/s by the end

#                     # old option - is this slower?
#                     parts = line.split()
#                     t = int(parts[2])
#                     progress.n = t # don't do progress.refresh here! it's very slow
                    
#                     if t > max_frames:
#                         break

#                 arraysize += 4 * 4
            
#             else:
#                 parts = line.split()
#                 t = int(parts[2])

#                 if t % nth_timesteps == 0:
#                     yield line
#                     arraysize += 4 * 4
                    
#                 if max_frames == None:
#                     progress.update()
#                 else:
#                     progress.n = t # don't do progress.refresh here! it's very slow
                    
#                     if t > max_frames:
#                         break

#             progress.set_description(f'loaded {common.format_bytes(arraysize)}')

# def read_file_tqdm(file):
#     arraysize = 0
#     with open(file) as f:
#         for line in (progress := tqdm.tqdm(f, total=1e10)):
#             yield line
            
#             arraysize += 4 * 4
#             progress.set_description(f'loaded {common.format_bytes(arraysize)}')


def go(infile, outfile, L, pack_frac_given, particle_diameter, dt=None, nth_timestep=1, max_time=None,
       extra_source_file=None, skiprows=None, quiet=False, multiply_time_by=None):
    print(f'loading raw file, last modified {common.get_last_modified_time(infile)} ago')

    density = 4/np.pi * pack_frac_given / particle_diameter**2
    expected_particles_per_frame = density * L**2
    if max_time:
        max_frame = max_time * dt
        max_rows = int((max_frame + 1) * expected_particles_per_frame)
        if not quiet: print('max_rows', max_rows)

        max_items = max_rows * 3
    else:
        max_rows = None
        max_items = -1

    t0 = time.time()
    # data = np.loadtxt(infile)
    # data = np.loadtxt(infile, dtype=np.float32)
    # data = np.loadtxt(read_file_skip_timesteps(infile, nth_timesteps, max_frames, expected_particles_per_frame), dtype=np.float32)
    if infile.endswith('.txt'):
        data = np.loadtxt(infile, max_rows=max_rows, dtype=np.float32)
    elif infile.endswith('.bin'):
        data = np.fromfile(infile, dtype=np.float32, count=max_items) # there is a parameter to this function for not loading the whole array
        data = data.reshape((-1, 4)) # was 3 but now we save ID
    t1 = time.time()
    if not quiet: print(data.shape, data.dtype,  f'loaded {common.format_bytes(data.nbytes)} in {t1-t0:.1f}s at {common.format_bytes(data.nbytes/(t1-t0))}/sec')

    if dt == None:
        assert multiply_time_by == None
        assert 'mixt' not in infile

        # dt = None, times = [0, 0.5, 1, 1.5, ...]
        times = np.unique(data[:, 2])
        dt = times[1] - times[0]
        multiply_time_by = 1/dt
        # dt = 0.5, times = [0, 1, 2, 3, ...]
        
    if multiply_time_by:
        if not quiet: print('doing multiply time by')
        data[:, 2] *= multiply_time_by
    
    # discard the non-nth step rows
    if nth_timestep > 1:
        if not quiet: print('doing nth-time')
        data = data[data[:, 2] % nth_timestep == 0, :]
        dt *= nth_timestep
        # we gotta readjust the timestep number now we discarded some timesteps
        data[:, 2] /= nth_timestep

    if skiprows:
        #### there is probably a way of doing this with loadtxt/fromfile that would be quicker
        if not quiet: print(f'doing skiprows ({skiprows})')
        if not quiet: print(data.shape)
        data = data[skiprows:, :]
        if not quiet: print(data.shape)
        # now need to reset the time column
        data[:, 2] -= data[:, 2].min()
        if not quiet: print()

    # we should delete the last timestep cause it may be incomplete
    if not quiet: print('trimming end')
    last_timestep = data[:, 2].max()
    data = data[data[:, 2] != last_timestep, :]
    
    if not quiet: print('xy min max', data[:, 0].min(), data[:, 0].max(), data[:, 1].min(), data[:, 1].max())
    
    # print(f'loaded in {t1-t0:.0f}s. shape', data.shape, common.arraysize(data))
    # num_timesteps = data[:, 2].max()+1
    times = np.unique(data[:, 2])
    num_timesteps = times.size
    if not quiet: print(times, num_timesteps)
    if not quiet: print(f'{num_timesteps:.0f} timesteps, {data[:, 2].max()*dt/60/60:.1f} hours')
    assert num_timesteps > 30

    density = common.calc_density(data, L, L, dimension=2)
    print('density', density)
    pack_frac_calced = common.calc_pack_frac(data, particle_diameter, L, L, dimension=2)
    if not quiet: print(f'pack_frac_calced={pack_frac_calced:.4f}, pack_frac_given={pack_frac_given:.4f}')
    if data.shape[0] > 1000 * 50:
        assert np.isclose(pack_frac_calced, pack_frac_given, rtol=0.1), f'pack frac calced {pack_frac_calced}, given {pack_frac_given}'

    data[:, 2] -= data[:, 2].min()
    assert data[:, 2].min() == 0

    if '_pot' in infile:
        if not quiet: print('shifting')
        data[:, 0] += L/2
        data[:, 1] += L/2
        if not quiet: print('xy min max', data[:, 0].min(), data[:, 0].max(), data[:, 1].min(), data[:, 1].max())
    else:
        if not quiet: print('periodic wrapping')
        # the data should already be inside a periodic box
        # but sometimes just a small number of points seem to be slightly outside the border
        # maybe because in the sim the saving is done before the periodic wrapping or something?
        data[:, 0] = data[:, 0] % L
        data[:, 1] = data[:, 1] % L

        assert L == L # for the mo

        if not quiet: print('x,y min,max', data[:, 0].min(), data[:, 0].max(), data[:, 1].min(), data[:, 1].max())
        keep = ( data[:, 0] >= 0 ) & ( data[:, 0] <= L ) &  ( data[:, 1] >= 0 ) & ( data[:, 1] <= L ) # I think they should be "<" not "<=" but for some reason this is failing on sim_nohydro_034_L1280
        assert keep.sum() == keep.size

    common.save_data(f'particle_detection/data/particles_{outfile}.npz',
        particles=data,
        time_step=dt, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
        window_size_x=L, window_size_y=L, max_time_hours=round(last_timestep*dt/60/60, 2),
        source_file=infile, density=density, extra_source_file=extra_source_file,
        # quiet=quiet,
    )
    
    if not '_pot' in infile:
        data_unwrap = common.periodic_unwrap(data, 2, [0, 1], [L, L], quiet=quiet)
        common.save_data(f'particle_detection/data/particles_{outfile}_unwrap.npz',
            particles=data_unwrap,
            time_step=dt, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
            window_size_x=L, window_size_y=L, max_time_hours=round(last_timestep*dt/60/60, 2),
            source_file=infile, density=density, extra_source_file=extra_source_file,
            quiet=quiet,
        )
        common.save_data(f'particle_linking/data/trajs_{outfile}_unwrap.npz',
            particles=data_unwrap,
            time_step=dt, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
            window_size_x=L, window_size_y=L, max_time_hours=round(last_timestep*dt/60/60, 2),
            source_file=infile, density=density, extra_source_file=extra_source_file,
            quiet=quiet,
        )

    # this used to be here but we can not possibly get the trajs straight from the simulation
    # without either unwrapping them, or cutting them at the periodic box crossings
    # if data.shape[1] == 4:
    #     common.save_data(f'particle_linking/data/trajs_{outfile}.npz',
    #         particles=data,
    #         time_step=dt, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
    #         window_size_x=L, window_size_y=L, max_time_hours=round(last_timestep*dt/60/60, 2),
    #         source_file=infile, density=density, extra_source_file=extra_source_file,
    #     )
    
    if data.size > 5e7 or True:
        end_timestep = data[:, 2].max() // 8
        data_small = data[data[:, 2] < end_timestep, :]
        data_small_unwrap = data_unwrap[data_unwrap[:, 2] < end_timestep, :]

        common.save_data(f'particle_detection/data/particles_{outfile}_div8.npz',
            particles=data_small,
            time_step=dt, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
            window_size_x=L, window_size_y=L, max_time_hours=round(end_timestep*dt/60/60, 2),
            source_file=infile, density=density, extra_source_file=extra_source_file,
            quiet=quiet,
        )

        if not '_pot' in infile:
            common.save_data(f'particle_detection/data/particles_{outfile}_unwrap_div8.npz',
                particles=data_small_unwrap,
                time_step=dt, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
                window_size_x=L, window_size_y=L, max_time_hours=round(end_timestep*dt/60/60, 2),
                source_file=infile, density=density, extra_source_file=extra_source_file,
                quiet=quiet,
            )
            common.save_data(f'particle_linking/data/trajs_{outfile}_unwrap_div8.npz',
                particles=data_small_unwrap,
                time_step=dt, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
                window_size_x=L, window_size_y=L, max_time_hours=round(end_timestep*dt/60/60, 2),
                source_file=infile, density=density, extra_source_file=extra_source_file,
                quiet=quiet,
            )
        # this used to be here but we can not possibly get the trajs straight from the simulation
        # without either unwrapping them, or cutting them at the periodic box crossings
        # if data.shape[1] == 4:
        #     common.save_data(f'particle_linking/data/trajs_{outfile}_div8.npz',
        #         particles=data_small,
        #         time_step=dt, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
        #         window_size_x=L, window_size_y=L, max_time_hours=round(end_timestep*dt/60/60, 2),
        #         source_file=infile, density=density, extra_source_file=extra_source_file,
        #     )

            
        end_timestep = data[:, 2].max() // 64
        data_small = data_small[data_small[:, 2] < end_timestep, :]
        data_small_unwrap = data_unwrap[data_unwrap[:, 2] < end_timestep, :]

        common.save_data(f'particle_detection/data/particles_{outfile}_div64.npz',
            particles=data_small,
            time_step=dt, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
            window_size_x=L, window_size_y=L, max_time_hours=round(end_timestep*dt/60/60, 2),
            source_file=infile, density=density, extra_source_file=extra_source_file,
            quiet=quiet,
        )

        if not '_pot' in infile:
            common.save_data(f'particle_detection/data/particles_{outfile}_unwrap_div64.npz',
                particles=data_small_unwrap,
                time_step=dt, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
                window_size_x=L, window_size_y=L, max_time_hours=round(end_timestep*dt/60/60, 2),
                source_file=infile, density=density, extra_source_file=extra_source_file,
                quiet=quiet,
            )
            common.save_data(f'particle_linking/data/trajs_{outfile}_unwrap_div64.npz',
                particles=data_small_unwrap,
                time_step=dt, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
                window_size_x=L, window_size_y=L, max_time_hours=round(end_timestep*dt/60/60, 2),
                source_file=infile, density=density, extra_source_file=extra_source_file,
                quiet=quiet,
            )
        # this used to be here but we can not possibly get the trajs straight from the simulation
        # without either unwrapping them, or cutting them at the periodic box crossings
        # if data.shape[1] == 4:
        #     common.save_data(f'particle_linking/data/trajs_{outfile}_div8.npz',
        #         particles=data_small,
        #         time_step=dt, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
        #         window_size_x=L, window_size_y=L, max_time_hours=round(end_timestep*dt/60/60, 2),
        #         source_file=infile, density=density, extra_source_file=extra_source_file,
        #     )
    print()
    
# 0.34
# go('raw_data/brennan/spec_softetakt_long_run_dtau_0.025_nsave_2.suspension_phi_0.34_L_320_modified.txt', 'brennan_hydro_034', 320, 320)
# go('raw_data/brennan/noHydro2D_Leim_run_dt_0.0125_nsave_40.suspension_phi_0.34_L_640_modified.txt', 'brennan_nohydro_034', 640, 320)
# 0.66
# go('/data2/acarter/Spectral_Sophie_Boxes/data/spec_softetakt_long_run_dtau_0.025_nsave_4.suspension_phi_0.66_L_288_modified.txt', 'brennan_hydro_066', 288, 288)
# go('/data2/acarter/Spectral_Sophie_Boxes/data/noHydro2D_Leim_run_dt_9.765625e5_nsave_256_long.suspension_phi_0.66_L_640_eq_modified.txt', 'brennan_nohydro_066', 640, 288)

# go('/data2/acarter/sim/RigidMultiblobsWall/Lubrication/Lubrication_Examples/Monolayer/data/nohydro2D_L1280_dt4.suspension_phi_0.34_L_1280_modified.txt', 'sim_nohydro_034_L1280', 1280, 1280, 4)
# go('/data2/acarter/sim/RigidMultiblobsWall/Lubrication/Lubrication_Examples/Monolayer/data/nohydro2D_L640_dt1.suspension_phi_0.34_L_640_modified.txt',   'sim_nohydro_034_L640',  640,  640,  1)

# from brennnan just before paper
# go('raw_data/brennan/spec_softetakt_dt_0.25_nsave_2.suspension_phi_0.016_L_320_modified.txt', 'brennan_hydro_002_L320', 320, 320, 0.5, 0.016, 2.972)
# go('raw_data/brennan/spec_softetakt_dt_0.25_nsave_2.suspension_phi_0.016_L_640_modified.txt', 'brennan_hydro_002_L640', 640, 640, 0.5, 0.016, 2.972)
# go('raw_data/brennan/spec_softetakt_dt_0.25_nsave_2.suspension_phi_0.114_L_320_modified.txt', 'brennan_hydro_011_L320', 320, 320, 0.5, 0.114, 2.972)
# go('raw_data/brennan/spec_softetakt_dt_0.25_nsave_64.suspension_phi_0.016_L_320_modified.txt', 'brennan_hydro_002_L320_longer', 320, 320, 16, 0.016, 2.972)

# from anubis
# go('raw_data/anubis/nohydro2D_L640_dt0.5_s2.972.suspension_phi0.016_L640_s2.972.bin', 'sim_nohydro_002_L640',        640, 640, 0.5, 0.016, 2.972)
# go('raw_data/anubis/nohydro2D_L640_dt16_s2.972.suspension_phi0.016_L640_s2.972.bin',  'sim_nohydro_002_L640_longer', 640, 640, 16,  0.016, 2.972)
# go('raw_data/anubis/nohydro2D_L160_dt0.5_s2.972.suspension_phi0.114_L160_s2.972.bin', 'sim_nohydro_011_L160',        160, 160, 0.5, 0.114, 2.972)
# go('raw_data/anubis/nohydro2D_L160_dt16_s2.972.suspension_phi0.114_L160_s2.972.bin',  'sim_nohydro_011_L160_longer', 160, 160, 16,  0.114, 2.972)
# from anubis after first paper submission
# go('raw_data/anubis/nohydro2D_L1280_dt64_s2.972.suspension_phi0.114_L1280_s2.972.bin',  'sim_nohydro_011_L1280_longer', 1280, 1280, 64,  0.114, 2.972)

datas = [
    # L  dt  phi
    # (320, 1, 0.01),
    # (640, 1, 0.01),
    # (1280, 4, 0.01),
    # (160, 1, 0.01),
    # (160, 4, 0.1),
    # (320, 4, 0.1),
    # (640, 8, 0.1, '', 1, None), # this was the one that was done again
    # (640, 0.5, 0.1, '_short'),
    # (320, 8, 0.1, '', 1, 24*60*60), # this is 24h
    # (320, 8, 0.01, '', 1, 24*60*60),
    # (640, 8, 0.01, '', 1, 24*60*60),
    # (320, 8, 0.1, '', 1, 24*60*60),
    # (640, 8, 0.1, '', 1, 24*60*60),
    # (544, 8, 0.1, '', 1, 24*60*60),
    # (544, 0.5, 0.1, '_dt2', 4, 24*60*60),
    # (800, 8, 0.02, '', 1, 24*60*60),
    # (320, 8, 0.114, '_long', 1, None, 2.972),
    # (640, 0.5, 0.10, '', 1, None, 2.79),
    # (640, 8, 0.10, '_longer', 2, None, 2.79),

    # (320, 0.5, 0.114, '',        1, None, 2.972), # sim_nohydro_011_L320
    # (320, 0.5, 0.016, '',        1, None, 2.972), # sim_nohydro_002_L320
    # (320, 16, 0.114, '_longer',  1, None, 2.972), # sim_nohydro_011_L320_longer
    # (320, 32, 0.114, '_longest', 1, None, 2.972), # sim_nohydro_011_L320_longest
    # (320, 16, 0.016, '_longer',  1, None, 2.972), # sim_nohydro_002_L320_longer
    # (640, 0.5, 0.016, '',        1, None, 2.972), # sim_nohydro_002_L640
    # (640, 0.5, 0.114, '',        1, None, 2.972), # sim_nohydro_011_L640
    # (640, 16, 0.114, '_longer',  1, None, 2.972), # sim_nohydro_011_L640_longer

    # after first submission of paper
    # (1280, 0.5, 0.114, '',        1, None, 2.972), # sim_nohydro_011_L1280
    # (1280, 0.5, 0.016, '',        1, None, 2.972), # sim_nohydro_002_L1280
    # (1280, 64, 0.016, '_longer',        1, None, 2.972), # sim_nohydro_011_L1280 # not yet run but this line is here for when it is
    
    # reimport
    # (544, 0.5, 0.1, '',        1, None, 2.79), # sim_nohydro_011_L320
]

# if __name__ == '__main__':
#     for L, dt, phi, suffix, nth_timestep, max_time, particle_diameter in datas:
#         phistr = f'{phi*100:.0f}'.zfill(3)
#         # go(f'/data2/acarter/sim/RigidMultiblobsWall/Lubrication/Lubrication_Examples/Monolayer/data/nohydro2D_L{L}_dt{dt}.suspension_phi_{phi}_L_{L}_modified.txt', f'sim_nohydro_{phistr}_L{L}{suffix}',  L, L, dt, phi, nth_timestep, max_time)
        
#         if particle_diameter == 2.972:
#             # new ones
#             file = f'/data2/acarter/sim/RigidMultiblobsWall/Lubrication/Lubrication_Examples/Monolayer/data/nohydro2D_L{L}_dt{dt}_s2.972.suspension_phi{phi}_L{L}_s2.972.bin'
#         else:
#             # old ones
#             file = f'/data2/acarter/sim/RigidMultiblobsWall/Lubrication/Lubrication_Examples/Monolayer/data/nohydro2D_L{L}_dt{dt}.suspension_phi_{phi}_L_{L}_modified.txt'
        
#         go(
#             infile = file,
#             outfile = f'sim_nohydro_{phistr}_L{L}{suffix}',
#             L=L, 
#             dt=dt, 
#             pack_frac_given=phi, nth_timestep=nth_timestep, max_time=max_time, particle_diameter=particle_diameter
#         )


# new times (not frames) method
datas2 = [
    # (640, 16, 0.5, 0.114, '_mixt', 1, None, 2.972),
    # (1280, 64, 0.5, 0.114, '_mixt', 1, None, 2.972),
]

# if __name__ == '__main__':
#     for L, t1, t2, phi, suffix, nth_timestep, max_time, particle_diameter in datas2:
#         phistr = f'{phi*100:.0f}'.zfill(3)
#         # go(f'/data2/acarter/sim/RigidMultiblobsWall/Lubrication/Lubrication_Examples/Monolayer/data/nohydro2D_L{L}_dt{dt}.suspension_phi_{phi}_L_{L}_modified.txt', f'sim_nohydro_{phistr}_L{L}{suffix}',  L, L, dt, phi, nth_timestep, max_time)
#         go(
#             # new ones:
#         f'/data2/acarter/sim/RigidMultiblobsWall/Lubrication/Lubrication_Examples/Monolayer/data/nohydro2D_L{L}_t{t1}_{t2}.suspension_phi{phi}_L{L}_s2.972.bin',
#         # old ones:
#         #    f'/data2/acarter/sim/RigidMultiblobsWall/Lubrication/Lubrication_Examples/Monolayer/data/nohydro2D_L{L}_dt{dt}.suspension_phi_{phi}_L_{L}_modified.txt',
#         f'sim_nohydro_{phistr}_L{L}{suffix}',
#         L=L,
#         dt=1,
#         pack_frac_given=phi, nth_timestep=nth_timestep, max_time=max_time, particle_diameter=particle_diameter)


##### mesu

def go_mesu(filepath, suffix='', skip_rsync=False, quiet=False, **kwargs):

    filename = filepath.split('/')[-1]
    if not skip_rsync:
        rsync_command = ['rsync', f'cartera@login.mesu.sorbonne-universite.fr:{filepath}', f'raw_data/mesu/{filename}', '--progress']
        if not quiet: print('launching', ' '.join(rsync_command))
        subprocess.run(rsync_command, check=True)

    filename_no_ext = filename.split('.bin')[0]

    L   = int  (filename_no_ext.split('_L' )[1].split('_')[0])
    phi = float(filename_no_ext.split('phi')[1].split('_')[0])
    hydro = filename_no_ext.split('_')[0]

    if 'sigma' in filename_no_ext:
        particle_diameter = float(filename_no_ext.split('_sigma')[1].split('_')[0])
    else:
        particle_diameter = 2.972

    print(f'phi {phi}')
    
    go(
        f'raw_data/mesu/{filename}',
        f'sim_{hydro}_{phi}_L{L}{suffix}',
        L=L,
        pack_frac_given=phi,
        particle_diameter=particle_diameter,
        extra_source_file=filepath,
        quiet=quiet,
        **kwargs
    )

processes = []
process_names = []

def go_mesu_subprocess(*args, **kwargs):
    # kwargs['quiet'] = True
    # tqdm.__init__ = functools.partialmethod(tqdm.__init__, disable=True)
    # os.environ['TQDM_DISABLE'] = '1'
    task = multiprocessing.Process(target=go_mesu, args=args, kwargs=kwargs)
    processes.append(task)
    process_names.append(args[0][28:])
    task.start()

    # chatgpt reccomends using concurrent.futures.ProcessPoolExecutor or multiprocessing.Pool for error handling

if __name__ == '__main__':
    pass
    # hydro/nohydro plot
    # go_mesu('/store/cartera/2d_monolayer/hydro_phi0.114_L640.bin', '_mixt')
    # go_mesu('/store/cartera/2d_monolayer/nohydro_phi0.114_L640.bin', '_mixt', skiprows=int(90e6)),
    # go_mesu('/store/cartera/2d_monolayer/hydro_t0.5_phi0.114_L640.bin', '', multiply_time_by=1/0.5, dt=0.5)
    # go_mesu('/store/cartera/2d_monolayer/hydro_t16_phi0.114_L640.bin', suffix='_longer', multiply_time_by=1/16, dt=16)
    # go_mesu('/store/cartera/2d_monolayer/hydro_t0.5_phi0.114_L320.bin', '', multiply_time_by=1/0.5, dt=0.5)
    # go_mesu('/store/cartera/2d_monolayer/hydro_t16_phi0.114_L320.bin', suffix='_longer', multiply_time_by=1/16, dt=16)

    # go_mesu('/store/cartera/2d_monolayer/nohydro_phi0.016_L1280.bin', '_mixt')
    # go_mesu('/store/cartera/2d_monolayer/nohydro_phi0.114_L1280.bin', '_mixt')
    # go_mesu('/store/cartera/2d_monolayer/nohydro_phi0.114_L1280.bin', '_mixtmesu')
    # go_mesu('/store/cartera/2d_monolayer/nohydro_phi0.114_L320.bin', '_mixt')

    # go_mesu('/store/cartera/2d_monolayer/nohydro_test_mixt_phi0.114_L320.bin', '_test_mixt')
    # go_mesu('/store/cartera/2d_monolayer/nohydro_test_singlet_phi0.114_L320.bin', '_test_singlet')

    # go_mesu('/store/cartera/2d_monolayer/hydro_t0.5_pot_phi0.114_L320.bin', suffix='_pot', multiply_time_by=1/0.5, dt=0.5)
    # go_mesu('/store/cartera/2d_monolayer/hydro_t16_pot_phi0.114_L320.bin', suffix='_pot_longer', multiply_time_by=1/16, dt=16)
    # go_mesu('/store/cartera/2d_monolayer/hydro_t0.5_pot_phi0.114_L640.bin', suffix='_pot', multiply_time_by=1/0.5, dt=0.5)
    # go_mesu('/store/cartera/2d_monolayer/hydro_t16_pot_phi0.114_L640.bin', suffix='_pot_longer', multiply_time_by=1/16, dt=16)
    # go_mesu('/store/cartera/2d_monolayer/hydro_t0.5_pot_phi0.114_L1280.bin', suffix='_pot', multiply_time_by=1/0.5, dt=0.5, nth_timestep=4)

    # go_mesu('/store/cartera/2d_monolayer/nohydro_t0.5_pot_phi0.114_L640.bin', suffix='_pot', multiply_time_by=1/0.5, dt=0.5)

    # missing one from before
    # go_mesu('/store/cartera/2d_monolayer/nohydro_t0.5_phi0.114_L1280.bin', suffix='_pot', multiply_time_by=1/0.5, dt=0.5, skip_rsync=True)

    # nohydro xy pot conf for fkt windowing
    # go_mesu('/store/cartera/2d_monolayer/nohydro_t0.5_pot_phi0.114_L640.bin', suffix='_pot', multiply_time_by=1/0.5, dt=0.5)
    # go_mesu('/store/cartera/2d_monolayer/nohydro_t16_pot_phi0.114_L640.bin', suffix='_pot_longer', multiply_time_by=1/16, dt=16) # unfinished
    # go_mesu('/store/cartera/2d_monolayer/nohydro_t0.5_pot_phi0.016_L640.bin', suffix='_pot', multiply_time_by=1/0.5, dt=0.5)
    # go_mesu('/store/cartera/2d_monolayer/nohydro_t16_pot_phi0.016_L640.bin', suffix='_pot_longer', multiply_time_by=1/16, dt=16)
    # go_mesu('/store/cartera/2d_monolayer/nohydro_t16_longest_pot_phi0.016_L640.bin', suffix='_pot_longest', multiply_time_by=1/16, dt=16)

    # others for fig S6
    # go_mesu('/store/cartera/2d_monolayer/nohydro_t0.5_pot_phi0.114_L320.bin', suffix='_pot', multiply_time_by=1/0.5, dt=0.5)
    # go_mesu('/store/cartera/2d_monolayer/nohydro_t0.5_pot_phi0.114_L1280.bin', suffix='_pot', multiply_time_by=1/0.5, dt=0.5)
    # go_mesu('/store/cartera/2d_monolayer/nohydro_t16_pot_phi0.114_L320.bin', suffix='_pot_longer', multiply_time_by=1/16, dt=16)
    # go_mesu('/store/cartera/2d_monolayer/nohydro_t64_pot_phi0.114_L1280.bin', suffix='_pot_longer')
    
    # for PNV
    # go_mesu('/store/cartera/2d_monolayer/nohydro_t0.5_short_phi0.01_L320.bin', suffix='_short', multiply_time_by=1/0.5, dt=0.5)
    # go_mesu('/store/cartera/2d_monolayer/nohydro_t0.5_short_phi0.1_L320.bin', suffix='_short', multiply_time_by=1/0.5, dt=0.5)
    # go_mesu('/store/cartera/2d_monolayer/nohydro_t0.5_short_phi0.3_L320.bin', suffix='_short', multiply_time_by=1/0.5, dt=0.5)
    # go_mesu('/store/cartera/2d_monolayer/nohydro_t0.5_short_phi0.5_L320.bin', suffix='_short', multiply_time_by=1/0.5, dt=0.5)
    # go_mesu('/store/cartera/2d_monolayer/nohydro_short_phi0.2_L320.bin', suffix='_short', multiply_time_by=1/0.5, dt=0.5)
    # go_mesu('/store/cartera/2d_monolayer/nohydro_short_phi0.2_L320.bin', suffix='_short', multiply_time_by=1/0.5, dt=0.5)
    # for phi in [0.01, 0.1, 0.2]:
    #     go_mesu_subprocess(f'/store/cartera/2d_monolayer/nohydro_short_phi{phi}_L320.bin', suffix='_short')
    
    # zconf
    # go_mesu('/store/cartera/2d_monolayer/hydro_t0.5_pot_zconf_phi0.114_L640.bin', suffix='_pot_zconf', skip_rsync=True)

    # Eleanor drift vs theta interacting
    # go_mesu('/store/cartera/2d_monolayer/hydro_t0.5_theta10_phi0.02_L640.bin', suffix='_theta10')
    # go_mesu('/store/cartera/2d_monolayer/hydro_t0.5_theta10_phi0.04_L640.bin', suffix='_theta10')
    # go_mesu('/store/cartera/2d_monolayer/hydro_t0.5_theta10_phi0.06_L640.bin', suffix='_theta10')
    # go_mesu('/store/cartera/2d_monolayer/hydro_t0.5_theta10_phi0.08_L640.bin', suffix='_theta10')
    # go_mesu('/store/cartera/2d_monolayer/hydro_t0.5_theta10_phi0.1_L640.bin', suffix='_theta10')
    # for phi in [0.02, 0.04, 0.06, 0.08, 0.1]:
    #     go_mesu_subprocess(f'/store/cartera/2d_monolayer/nohydro_t0.5_theta10_phi{phi}_L640.bin', suffix='_theta10', skip_rsync=True)

    # CoM (first one also used by PNV)
    # go_mesu('/store/cartera/2d_monolayer/nohydro_short_phi0.1_L320.bin', suffix='_short', multiply_time_by=1/0.5, dt=0.5, skip_rsync=True)
    # go_mesu('/store/cartera/2d_monolayer/nohydro_short_phi0.1_L160.bin', suffix='_short', multiply_time_by=1/0.5, dt=0.5, skip_rsync=True)
    # go_mesu('/store/cartera/2d_monolayer/nohydro_short_phi0.1_L640.bin', suffix='_short', multiply_time_by=1/0.5, dt=0.5, skip_rsync=True)
    # go_mesu('/store/cartera/2d_monolayer/nohydro_short_phi0.1_L1280.bin', suffix='_short', multiply_time_by=1/0.5, dt=0.5)
    # for L in np.logspace(np.log10(100), np.log10(1000), num=10):
    #     L = int(L)
    #     go_mesu(f'/store/cartera/2d_monolayer/nohydro_com_phi0.1_L{L}.bin', suffix='_CoM', multiply_time_by=1/0.5, dt=0.5)
    # for L in np.logspace(np.log10(10), np.log10(1000), num=20):
    #     L = int(L)
    #     go_mesu_subprocess(f'/store/cartera/2d_monolayer/nohydro_t1e4_phi0.1_L{L}.bin', suffix='_t1e4')
    
    # with multiprocessing.Pool(cores) as pool:
    #     task = functools.partial(calc_centre_of_mass_onepoint_for_single_timeorigin, groupsize)

        # all_data = []

        # for time_origin in tqdm.trange(num_timesteps-1, desc='preparing data'):
        #     data_these_timesteps = get_data_at_timestep(data, time_origin)
        #     all_data.append(data_these_timesteps)
        #     print('size', common.arraysize(data_these_timesteps, mult=len(all_data)))

        # results = list(tqdm.tqdm(pool.imap(task, get_data_at_timesteps()), total=num_time_origins, desc=f'N={str(groupsize).ljust(3)}'))

        
    # for phi in [0.1]:
    for phi in [0.02, 0.04, 0.06, 0.08, 0.1]:
        go_mesu(f'/store/cartera/2d_monolayer/hydro_t0.5_10m_sigma5.944_theta10_phi{phi}_L1280.bin', suffix='_theta0_10m_sigma2x')
    #     go_mesu_subprocess  (f'/store/cartera/2d_monolayer/nohydro_t0.5_10m_theta10_phi{phi}_L640.bin', suffix='_theta10_10m')
    #     # go_mesu_subprocess(f'/store/cartera/2d_monolayer/hydro_t0.5_1h_theta0_phi{phi}_L640.bin', suffix='_theta0_1h', skip_rsync=True)
    

    # for L in np.logspace(np.log10(10), np.log10(1000), num=5):
    #     go_mesu(f'/store/cartera/2d_monolayer/hydro_t1e4_phi0.1_L{int(L)}.bin', suffix='_theta0_1e4')

    print(f'{len(processes)} tasks')
    for i, task in enumerate(tqdm.tqdm(processes, desc='processes')):
        task.join() # block
        print(f'task {i} : {process_names[i]} finished')

