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


def go(infile, outfile, L, max_time=None,
       extra_source_file=None, quiet=False, fix_timestep=1, metadata={}, binary_metadata={}):
    print(f'loading raw file, last modified {common.get_last_modified_time(infile)} ago')

    metadata['dimension'] = 3

    dt = metadata['t_save'] * fix_timestep
    particle_diameter = metadata['a'] * 2
    pack_frac_given = metadata['phi']

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
    
    assert infile.endswith('.bin')

    data_ryker = np.fromfile(infile, dtype=binary_metadata['dtype'], count=max_items) # there is a parameter to this function for not loading the whole array
    data_ryker = data_ryker.reshape((-1, 1+metadata['N_colloids']*3))
    print(f'data_ryker is {data_ryker.shape}, {data_ryker.dtype}, {common.arraysize(data_ryker)}')
    t1 = time.time()
    if not quiet: print(data_ryker.shape, data_ryker.dtype,  f'loaded {common.format_bytes(data_ryker.nbytes)} in {t1-t0:.1f}s at {common.format_bytes(data_ryker.nbytes/(t1-t0))}/sec')

    print('reshaping')

    num_expected_timesteps = data_ryker.shape[0]
    # assert binary_metadata['n_rows'] == data_ryker.shape[0], f"binary_metadata['n_rows'] = {binary_metadata['n_rows']} != data_ryker.shape[0] = {data_ryker.shape[0]}"

    num_colloids = metadata['N_colloids']
    assert num_colloids == binary_metadata['N']

    print(f'creating data, will be {(num_colloids*num_expected_timesteps, 5)} {num_colloids * num_expected_timesteps * 5 * 4/1e9:.1f}GB')
    assert data_ryker.nbytes + num_colloids * num_expected_timesteps * 5 * 4 < 100e9
    time.sleep(1) # so the printing buffer can be cleared
    data = np.full((num_colloids*num_expected_timesteps, 5), np.nan, dtype=np.float32)

    for t in tqdm.trange(num_expected_timesteps):
        xyz_this_timestep = data_ryker[t, 1:]
        starting_row_i = t*num_colloids
        data[starting_row_i:starting_row_i+num_colloids, [0, 1, 2]] = xyz_this_timestep.reshape((num_colloids, 3))
        data[starting_row_i:starting_row_i+num_colloids, 3] = t
        data[starting_row_i:starting_row_i+num_colloids, 4] = np.arange(num_colloids)
    #     for particle in range(metadata['N_colloids']):
    #         row_i = t*metadata['N_colloids'] + particle
    #         data[row_i, [0, 1, 2]]
    #         data_ryker[t, [1+3*particle, 1+3*particle+1, 1+3*particle+2]]
    #         data[row_i, [0, 1, 2]] = data_ryker[t, [1+3*particle, 1+3*particle+1, 1+3*particle+2]]
    #         data[row_i, 3] = t
    #         data[row_i, 4] = particle

    # # discard the non-nth step rows
    # if nth_timestep > 1:
    #     if not quiet: print('doing nth-time')
    #     data = data[data[:, time_column] % nth_timestep == 0, :]
    #     dt *= nth_timestep
    #     # we gotta readjust the timestep number now we discarded some timesteps
    #     data[:, time_column] /= nth_timestep

    # if skiprows:
    #     #### there is probably a way of doing this with loadtxt/fromfile that would be quicker
    #     if not quiet: print(f'doing skiprows ({skiprows})')
    #     if not quiet: print(data.shape)
    #     data = data[skiprows:, :]
    #     if not quiet: print(data.shape)
    #     # now need to reset the time column
    #     data[:, time_column] -= data[:, time_column].min()
    #     if not quiet: print()

    time_column = 3
    id_column = 4

    if fix_timestep > 1:
        data[:, time_column] /= fix_timestep

    data[:, time_column] -= data[:, time_column].min()
    assert data[:, time_column].min() == 0

    last_timestep = data[:, time_column].max()
    assert np.isfinite(last_timestep)
    
    if not quiet: print('xy min max', data[:, 0].min(), data[:, 0].max(), data[:, 1].min(), data[:, 1].max())
    
    # print(f'loaded in {t1-t0:.0f}s. shape', data.shape, common.arraysize(data))
    # num_timesteps = data[:, time_column].max()+1
    times, num_particles_at_frame = np.unique(data[:, time_column], return_counts=True)
    num_timesteps = times.size
    assert num_timesteps == num_expected_timesteps
    # if not quiet: print(times, num_timesteps)
    if not quiet: print(f'{num_timesteps:.0f} timesteps, {last_timestep*dt/60/60:.1f} hours')
    # assert num_timesteps > 30

    max_particles_at_frame = num_particles_at_frame.max()
    assert max_particles_at_frame < 1.5 * num_particles_at_frame.mean(), f'max particles at frame {max_particles_at_frame} avg particles at frame {num_particles_at_frame.mean():.1f}'


    density = common.calc_density(data, L, L, dimension=3)
    print('density', density)
    pack_frac_calced = common.calc_pack_frac(data, particle_diameter, L, L, dimension=3)
    if not quiet: print(f'pack_frac_calced={pack_frac_calced:.4f}, pack_frac_given={pack_frac_given:.4f}')
    if data.shape[0] > 1000 * 50:
        assert np.isclose(pack_frac_calced, pack_frac_given, rtol=0.1), f'pack frac calced {pack_frac_calced}, given {pack_frac_given}'


    if 'NBody' in metadata['solver_name']:

        outside = (data[:, 0] < 0) | (data[:, 0] > L) | (data[:, 0] < 1) | (data[:, 1] > L)
        print(f'{outside.sum()/outside.size:.1%} of particles under potential influence')

        # want to find how many particles ever feel the influence of the potential
        num_ever_outside = np.unique(data[outside, id_column]).size
        print(f'{num_ever_outside/num_particles_at_frame.mean():.1%} ever feel influence of potential')

    elif metadata['solver_name'] == 'DPStokes':
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
        **metadata
        # quiet=quiet,
    )
    
    do_unwrapping = not ('NBody' in metadata['solver_name'])

    if do_unwrapping:
        data_unwrap = common.periodic_unwrap(data, 3, [0, 1], [L, L], quiet=quiet)
        common.save_data(f'particle_detection/data/particles_{outfile}_unwrap.npz',
            particles=data_unwrap,
            time_step=dt, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
            window_size_x=L, window_size_y=L, max_time_hours=round(last_timestep*dt/60/60, 2),
            source_file=infile, density=density, extra_source_file=extra_source_file,
            quiet=quiet,
            **metadata
        )
        common.save_data(f'particle_linking/data/trajs_{outfile}_unwrap.npz',
            particles=data_unwrap,
            time_step=dt, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
            window_size_x=L, window_size_y=L, max_time_hours=round(last_timestep*dt/60/60, 2),
            source_file=infile, density=density, extra_source_file=extra_source_file,
            quiet=quiet,
            **metadata
        )
    else:
        common.save_data(f'particle_linking/data/trajs_{outfile}.npz',
            particles=data,
            time_step=dt, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
            window_size_x=L, window_size_y=L, max_time_hours=round(last_timestep*dt/60/60, 2),
            source_file=infile, density=density, extra_source_file=extra_source_file,
            quiet=quiet,
            **metadata
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
        end_timestep = data[:, time_column].max() // 8
        data_small = data[data[:, time_column] < end_timestep, :]

        common.save_data(f'particle_detection/data/particles_{outfile}_div8.npz',
            particles=data_small,
            time_step=dt, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
            window_size_x=L, window_size_y=L, max_time_hours=round(end_timestep*dt/60/60, 2),
            source_file=infile, density=density, extra_source_file=extra_source_file,
            quiet=quiet,
            **metadata
        )

        if do_unwrapping:
            data_small_unwrap = data_unwrap[data_unwrap[:, 2] < end_timestep, :]
            common.save_data(f'particle_detection/data/particles_{outfile}_unwrap_div8.npz',
                particles=data_small_unwrap,
                time_step=dt, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
                window_size_x=L, window_size_y=L, max_time_hours=round(end_timestep*dt/60/60, 2),
                source_file=infile, density=density, extra_source_file=extra_source_file,
                quiet=quiet,
                **metadata
            )
            common.save_data(f'particle_linking/data/trajs_{outfile}_unwrap_div8.npz',
                particles=data_small_unwrap,
                time_step=dt, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
                window_size_x=L, window_size_y=L, max_time_hours=round(end_timestep*dt/60/60, 2),
                source_file=infile, density=density, extra_source_file=extra_source_file,
                quiet=quiet,
                **metadata
            )
        else:
            common.save_data(f'particle_linking/data/trajs_{outfile}_div8.npz',
                particles=data_small,
                time_step=dt, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
                window_size_x=L, window_size_y=L, max_time_hours=round(end_timestep*dt/60/60, 2),
                source_file=infile, density=density, extra_source_file=extra_source_file,
                quiet=quiet,
                **metadata
            )

            
        end_timestep = data[:, time_column].max() // 64
        data_small = data_small[data_small[:, 2] < end_timestep, :]

        common.save_data(f'particle_detection/data/particles_{outfile}_div64.npz',
            particles=data_small,
            time_step=dt, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
            window_size_x=L, window_size_y=L, max_time_hours=round(end_timestep*dt/60/60, 2),
            source_file=infile, density=density, extra_source_file=extra_source_file,
            quiet=quiet,
            **metadata
        )

        if do_unwrapping:
            data_small_unwrap = data_unwrap[data_unwrap[:, 2] < end_timestep, :]
            common.save_data(f'particle_detection/data/particles_{outfile}_unwrap_div64.npz',
                particles=data_small_unwrap,
                time_step=dt, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
                window_size_x=L, window_size_y=L, max_time_hours=round(end_timestep*dt/60/60, 2),
                source_file=infile, density=density, extra_source_file=extra_source_file,
                quiet=quiet,
                **metadata
            )
            common.save_data(f'particle_linking/data/trajs_{outfile}_unwrap_div64.npz',
                particles=data_small_unwrap,
                time_step=dt, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
                window_size_x=L, window_size_y=L, max_time_hours=round(end_timestep*dt/60/60, 2),
                source_file=infile, density=density, extra_source_file=extra_source_file,
                quiet=quiet,
                **metadata
            )
        else:
            common.save_data(f'particle_linking/data/trajs_{outfile}_div64.npz',
                particles=data_small,
                time_step=dt, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
                window_size_x=L, window_size_y=L, max_time_hours=round(end_timestep*dt/60/60, 2),
                source_file=infile, density=density, extra_source_file=extra_source_file,
                quiet=quiet,
                **metadata
            )
            
    print()
    

import json

def go_mesu(directory, suffix='', skip_rsync=False, quiet=False, small=False, short=False, **kwargs):

    assert directory.endswith('/')
    if small:
        filepath = f'{directory}colloids_small{small}.bin'
        kwargs['fix_timestep'] = small
    elif short:
        filepath = f'{directory}colloids_short{short}.bin'
        suffix = f'_short{suffix}'
    else:
        filepath = f'{directory}colloids.bin'
        
    filepath_metadata = f'{directory}params.json'
    filepath_binary_metadata = f'{directory}binary_metadata.json'

    filename = filepath.split('/')[-2]
    json_filepath_local = f'raw_data/mesu/{filename}{suffix}.params.json'
    binary_json_filepath_local = f'raw_data/mesu/{filename}{suffix}.binary.json'
    colloids_filepath_local = f'raw_data/mesu/{filename}{suffix}.bin'
    kwargs['extra_source_file']=filepath

    if not skip_rsync:
        rsync_command = ['rsync', f'cartera@login.mesu.sorbonne-universite.fr:{filepath}', colloids_filepath_local, '--progress']
        if not quiet: print('launching', ' '.join(rsync_command))
        subprocess.run(rsync_command, check=True)

        rsync_command = ['rsync', f'cartera@login.mesu.sorbonne-universite.fr:{filepath_metadata}', json_filepath_local, '--progress']
        if not quiet: print('launching', ' '.join(rsync_command))
        subprocess.run(rsync_command, check=True)

        rsync_command = ['rsync', f'cartera@login.mesu.sorbonne-universite.fr:{filepath_binary_metadata}', binary_json_filepath_local, '--progress']
        if not quiet: print('launching', ' '.join(rsync_command))
        subprocess.run(rsync_command, check=True)
        
    go_after_rsync(suffix, kwargs, json_filepath_local, binary_json_filepath_local, colloids_filepath_local, quiet)

def go_after_rsync(suffix, kwargs, json_filepath_local, binary_json_filepath_local, colloids_filepath_local, quiet=False):
    with open(json_filepath_local) as json_file:
        metadata = json.load(json_file)
    print(metadata)
        
    with open(binary_json_filepath_local) as binary_json_file:
        binary_metadata = json.load(binary_json_file)
    print(binary_metadata)

    L   = int  (metadata['Lx'])
    assert metadata['Lx'] == metadata['Ly']
    phi = metadata['phi']
    # hydro = filename_no_ext.split('_')[0]
    hydro = 'hydro_' + metadata['solver_name'].lower()

    print(f'phi {phi}')

    
    go(
        colloids_filepath_local,
        f'ld_{hydro}_{phi}_L{L}{suffix}',
        L=L,
        quiet=quiet,
        metadata=metadata,
        binary_metadata=binary_metadata,
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
    
    # # for phi in [0.1]:
    # for phi in [0.02, 0.04, 0.06, 0.08, 0.1]:
    #     go_mesu(f'/store/cartera/2d_monolayer/hydro_t0.5_10m_sigma5.944_theta10_phi{phi}_L1280.bin', suffix='_theta0_10m_sigma2x')
    # #     go_mesu_subprocess  (f'/store/cartera/2d_monolayer/nohydro_t0.5_10m_theta10_phi{phi}_L640.bin', suffix='_theta10_10m')
    # #     # go_mesu_subprocess(f'/store/cartera/2d_monolayer/hydro_t0.5_1h_theta0_phi{phi}_L640.bin', suffix='_theta0_1h', skip_rsync=True)
    

    # # for L in np.logspace(np.log10(10), np.log10(1000), num=5):
    # #     go_mesu(f'/store/cartera/2d_monolayer/hydro_t1e4_phi0.1_L{int(L)}.bin', suffix='_theta0_1e4')

    # print(f'{len(processes)} tasks')
    # for i, task in enumerate(tqdm.tqdm(processes, desc='processes')):
    #     task.join() # block
    #     print(f'task {i} : {process_names[i]} finished')

    # for solver in ['DPStokes', 'NBody']:
    #     for dt in [50, 100, 150]:
    #         go_mesu(f'/store/cartera/libmobility_diffusion/solver_{solver}_N_7638_L_640_dt_{dt}_run_0/', suffix=f'_dt{dt}')

    # go_mesu('/store/cartera/libmobility_diffusion/solver_NBody_N_30552_L_1280_dt_125_t_604800_64_run_0/')
    # go_mesu('/store/cartera/libmobility_diffusion/solver_DPStokes_N_30552_L_1280_dt_125_t_604800_64_run_0/')
    # go_mesu('/store/cartera/libmobility_diffusion/solver_NBody_N_122205_L_2560_dt_125_t_604800_64_run_0/')
    # go_mesu('/store/cartera/libmobility_diffusion/solver_DPStokes_N_122205_L_2560_dt_125_t_604800_64_run_0/')
    

    # 5% pack frac
    # go_mesu('/store/cartera/libmobility_diffusion/solver_NBody_N_53599_L_2560_dt_125_t_604800_64_run_0/')

    # vacf
    # go_mesu('/store/cartera/libmobility_diffusion/solver_NBody_N_7638_L_640_dt_25_t_100_0.025_run_0/', suffix='_t100_0.025')
    # don't need this much data, 10s probably would be fine

    # libmobility paper
    # go_mesu('/store/cartera/libmobility_diffusion/solver_NBody_N_122205_L_2560_dt_125_t_3600_1_run_0/', suffix='_t1h_1')
    # go_mesu('/store/cartera/libmobility_diffusion/solver_NBody_open_N_122205_L_2560_dt_20_t_3600_1_run_0/', suffix='_t1h_1')
    # go_mesu('/store/cartera/libmobility_diffusion/solver_NBody_open_N_122205_L_2560_dt_20_t_28800_64_run_1/', suffix='_t8h_64')
    # libmobility bigboi
    go_after_rsync('', {}, 'raw_data/brennan/large_diffusion_data/params_partial.json', 'raw_data/brennan/large_diffusion_data/binary_metadata_partial.json', '/data2/acarter/toolbox/raw_data/brennan/large_diffusion_data/colloids_partial.bin')