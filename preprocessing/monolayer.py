import numpy as np
import common
import time, os
import tqdm
import subprocess
import multiprocessing
import functools
import common

def go(infile, outfile, L, pack_frac_given, particle_diameter, nth_timestep=1,
       extra_source_file=None, skiprows=None, quiet=False, multiply_time_by=None):
    print(f'loading raw file, last modified {common.get_last_modified_time(infile)} ago')

    expected_density = 4/np.pi * pack_frac_given / particle_diameter**2
    expected_particles_per_frame = expected_density * L**2

    time_column = 3
    id_column = 4
    time_column = 3

    assert infile.endswith('.npz')

    t0 = time.time()
    
    all_data = np.load(infile)
    data = all_data['particles']
    metadata = common.copy_not_particles(all_data)
    metadata['dimension'] = 3
    
    t1 = time.time()
    if not quiet: print(data.shape, data.dtype,  f'loaded {common.format_bytes(data.nbytes)} in {t1-t0:.1f}s at {common.format_bytes(data.nbytes/(t1-t0))}/sec')

    saved_time_step = all_data['saved_time_step']
    # if dt == None:
    #     assert multiply_time_by == None
    #     assert 'mixt' not in infile

    #     # dt = None, times = [0, 0.5, 1, 1.5, ...]
    #     times = np.unique(data[:, time_column])
    #     dt = times[1] - times[0]
    #     multiply_time_by = 1/dt
        # dt = 0.5, times = [0, 1, 2, 3, ...]
        
    # if multiply_time_by:
    #     if not quiet: print('doing multiply time by')
    #     data[:, time_column] *= multiply_time_by
    
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

    # we should delete the last timestep cause it may be incomplete
    # if not quiet: print('trimming end')
    last_timestep = data[:, time_column].max()
    # data = data[data[:, time_column] != last_timestep, :]
    
    if not quiet: print('xy min max', data[:, 0].min(), data[:, 0].max(), data[:, 1].min(), data[:, 1].max())
    
    # print(f'loaded in {t1-t0:.0f}s. shape', data.shape, common.arraysize(data))
    # num_timesteps = data[:, time_colum].max()+1
    times = np.unique(data[:, time_column])
    num_timesteps = times.size
    if not quiet: print(times, num_timesteps)
    if not quiet: print(f'{num_timesteps:.0f} timesteps, {data[:, time_column].max()*saved_time_step/60/60:.1f} hours')
    assert num_timesteps > 30

    density = common.calc_density(data, L, L, dimension=3)
    print('density', density)
    pack_frac_calced = common.calc_pack_frac(data, particle_diameter, L, L, dimension=3)
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
        time_step=saved_time_step, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
        window_size_x=L, window_size_y=L, max_time_hours=round(last_timestep*saved_time_step/60/60, 2),
        source_file=infile, density=density, extra_source_file=extra_source_file,
        **metadata
        # quiet=quiet,
    )
    # check that it saved properly
    np.load(f'particle_detection/data/particles_{outfile}.npz')
    
    if not '_pot' in infile:
        data_unwrap = common.periodic_unwrap(data, 3, [0, 1], [L, L], quiet=quiet)
        common.save_data(f'particle_detection/data/particles_{outfile}_unwrap.npz',
            particles=data_unwrap,
            time_step=saved_time_step, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
            window_size_x=L, window_size_y=L, max_time_hours=round(last_timestep*saved_time_step/60/60, 2),
            source_file=infile, density=density, extra_source_file=extra_source_file,
            quiet=quiet,
            **metadata
        )
        common.save_data(f'particle_linking/data/trajs_{outfile}_unwrap.npz',
            particles=data_unwrap,
            time_step=saved_time_step, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
            window_size_x=L, window_size_y=L, max_time_hours=round(last_timestep*saved_time_step/60/60, 2),
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
    
    end_timestep = data[:, time_column].max() // 8
    data_small = data[data[:, time_column] < end_timestep, :]
    assert data_small.size

    common.save_data(f'particle_detection/data/particles_{outfile}_div8.npz',
        particles=data_small,
        time_step=saved_time_step, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
        window_size_x=L, window_size_y=L, max_time_hours=round(end_timestep*saved_time_step/60/60, 2),
        source_file=infile, density=density, extra_source_file=extra_source_file,
        quiet=quiet,
        **metadata
    )

    if not '_pot' in infile:
        data_small_unwrap = data_unwrap[data_unwrap[:, 2] < end_timestep, :]
        assert data_small_unwrap.size
        common.save_data(f'particle_detection/data/particles_{outfile}_unwrap_div8.npz',
            particles=data_small_unwrap,
            time_step=saved_time_step, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
            window_size_x=L, window_size_y=L, max_time_hours=round(end_timestep*saved_time_step/60/60, 2),
            source_file=infile, density=density, extra_source_file=extra_source_file,
            quiet=quiet,
            **metadata
        )
        common.save_data(f'particle_linking/data/trajs_{outfile}_unwrap_div8.npz',
            particles=data_small_unwrap,
            time_step=saved_time_step, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
            window_size_x=L, window_size_y=L, max_time_hours=round(end_timestep*saved_time_step/60/60, 2),
            source_file=infile, density=density, extra_source_file=extra_source_file,
            quiet=quiet,
            **metadata
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

        
    end_timestep = data[:, time_column].max() // 64
    data_small = data_small[data_small[:, time_column] < end_timestep, :]

    common.save_data(f'particle_detection/data/particles_{outfile}_div64.npz',
        particles=data_small,
        time_step=saved_time_step, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
        window_size_x=L, window_size_y=L, max_time_hours=round(end_timestep*saved_time_step/60/60, 2),
        source_file=infile, density=density, extra_source_file=extra_source_file,
        quiet=quiet,
        **metadata
    )

    if not '_pot' in infile:
        data_small_unwrap = data_unwrap[data_unwrap[:, 2] < end_timestep, :]
        common.save_data(f'particle_detection/data/particles_{outfile}_unwrap_div64.npz',
            particles=data_small_unwrap,
            time_step=saved_time_step, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
            window_size_x=L, window_size_y=L, max_time_hours=round(end_timestep*saved_time_step/60/60, 2),
            source_file=infile, density=density, extra_source_file=extra_source_file,
            quiet=quiet,
            **metadata
        )
        common.save_data(f'particle_linking/data/trajs_{outfile}_unwrap_div64.npz',
            particles=data_small_unwrap,
            time_step=saved_time_step, particle_diameter=particle_diameter, pack_frac_given=pack_frac_given, pack_frac=pack_frac_calced,
            window_size_x=L, window_size_y=L, max_time_hours=round(end_timestep*saved_time_step/60/60, 2),
            source_file=infile, density=density, extra_source_file=extra_source_file,
            quiet=quiet,
            **metadata
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
    


# def go_mesu_npz(filepath, suffix='', skip_rsync=False, quiet=False, **kwargs):

#     filename = filepath.split('/')[-1]
#     if not skip_rsync:
#         rsync_command = ['rsync', f'cartera@login.mesu.sorbonne-universite.fr:{filepath}', f'raw_data/mesu/{filename}', '--progress']
#         if not quiet: print('launching', ' '.join(rsync_command))
#         subprocess.run(rsync_command, check=True)

#     filename_no_ext = filename.split('.bin')[0]

#     L   = int  (filename_no_ext.split('_L' )[1].split('_')[0])
#     phi = float(filename_no_ext.split('phi')[1].split('_')[0])
#     hydro = filename_no_ext.split('_')[0]

#     if 'sigma' in filename_no_ext:
#         particle_diameter = float(filename_no_ext.split('_sigma')[1].split('_')[0])
#     else:
#         particle_diameter = 2.972

#     print(f'phi {phi}')
    
#     go(
#         f'raw_data/mesu/{filename}',
#         f'sim_{hydro}_{phi}_L{L}{suffix}',
#         L=L,
#         pack_frac_given=phi,
#         particle_diameter=particle_diameter,
#         extra_source_file=filepath,
#         quiet=quiet,
#         **kwargs
#     )

processes = []
process_names = []

def go_mesu_subprocess(*args, **kwargs):
    # kwargs['quiet'] = True
    # tqdm.__init__ = functools.partialmethod(tqdm.__init__, disable=True)
    # os.environ['TQDM_DISABLE'] = '1'
    task = multiprocessing.Process(target=go_mesu_npz, args=args, kwargs=kwargs)
    processes.append(task)
    process_names.append(args[0][28:])
    task.start()

    # chatgpt reccomends using concurrent.futures.ProcessPoolExecutor or multiprocessing.Pool for error handling


##### mesu npz ######
def go_mesu_npz(filepath, suffix='', skip_rsync=False, quiet=False, **kwargs):

    filename = filepath.split('/')[-1]
    if not skip_rsync:
        rsync_command = ['rsync', f'cartera@login.mesu.sorbonne-universite.fr:{filepath}', f'raw_data/mesu/{filename}', '--progress']
        if not quiet: print('launching', ' '.join(rsync_command))
        subprocess.run(rsync_command, check=True)

    filename_no_ext = filename.split('.npz')[0]

    L   = int  (filename_no_ext.split('_L' )[1].split('_')[0])
    phi = float(filename_no_ext.split('phi')[1].split('_')[0])
    hydro = filename_no_ext.split('_')[0]

    if 'sigma' in filename_no_ext:
        particle_diameter = float(filename_no_ext.split('_sigma')[1].split('_')[0])
    else:
        particle_diameter = 2.972

    print(f'phi {phi}')

    outfile = f'sim_{hydro}_{phi}_L{L}'
    if '_pot' in filename_no_ext:
        outfile += '_pot'
    outfile += suffix
    # outfile = f'sim_{filename_no_ext}'
    
    go(
        f'raw_data/mesu/{filename}',
        outfile,
        L=L,
        pack_frac_given=phi,
        particle_diameter=particle_diameter,
        extra_source_file=filepath,
        quiet=quiet,
        **kwargs
    )
    
# for L in np.logspace(np.log10(100), np.log10(2000), num=10)[[-1]]:
#     L = int(L)
#     # go_mesu_npz(f'/store/cartera/2d_monolayer/hydro_t1e3_0.5_pot_phi0.114_L{L}.npz', suffix='_t1e3_0.5_pot')
#     go_mesu_npz(f'/store/cartera/2d_monolayer/hydro_t1e3_0.5_phi0.114_L{L}.npz', suffix='_t1e3_0.5')
    
# for phi in [0.02, 0.04, 0.06, 0.08, 0.1]:
#     go_mesu_npz(f'/store/cartera/2d_monolayer/hydro_t1e2_0.5_pot_phi{phi}_L500.npz', suffix='_t1e2_0.5_pot')
#     go_mesu_npz(f'/store/cartera/2d_monolayer/hydro_t1e2_0.5_phi{phi}_L500.npz', suffix='_t1e2_0.5')

def make_chunks_internal(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]
def make_chunks(lst, n):
    return list(make_chunks_internal(lst, n))
    
if __name__ == '__main__':
    # apparently not having the __name__ check can lead to "a process in the process pool was terminated abruptly while the future was running or pending"
    import functools, concurrent.futures
    with concurrent.futures.ProcessPoolExecutor() as executor:
        functions = []

        def go_mesu_subprocess2(*args, **kwargs):
            bound = functools.partial(go_mesu_npz, *args, **kwargs)
            functions.append(bound)
                       
        for L, t_max, t_max_name in [
            # ( 320, 60*60*16, '16h'),
            (1280, 60*60* 1,  '1h'),
            # ( 320, 60*4, '4m'),
            # ( 320, 60*1, '1m'),
        ]:
            if True:
                for seed in range(1, 11):
                    if seed != 8:
                        continue
                    # if t_max != 60*15:
                    #     continue
                    # if L != 320:
                    #     continue
                    # go_mesu_npz(f'/store/cartera/2d_monolayer/nohydro_t{t_max_name}_0.5_seed{seed}_phi0.114_L{L}.npz', suffix=f'_seed{seed}_t{t_max_name}_0.5')#, skip_rsync=True)
                    # go_mesu_npz(
                    #     f'/store/cartera/2d_monolayer/nohydro_phi0.114_L{L}_seed{seed}_t{t_max_name}_0.5.npz',
                    #     suffix=f'_seed{seed}_t{t_max_name}_0.5',
                    #     # skip_rsync=True
                    # )


        for chunk_i, chunk in enumerate(chunks := make_chunks(functions, 10)):
            futures = [executor.submit(func) for func in chunk]
            for future in concurrent.futures.as_completed(futures):
                # try:
                    result = future.result()  # Raises exception from subprocess if any
                    # print(result)
                # except Exception as e:
                #     raise RuntimeError(f"Subprocess raised: {e}")
            print(f'finished chunk {chunk_i+1}/{len(chunks)}')


    # go_mesu_npz(
    #     f'/store/cartera/2d_monolayer/hydro_phi0.114_L320_pot_seed1_t1h_0.5.npz',
    #     suffix=f'_seed1_t1h_0.5',
    #     skip_rsync=True
    # )
    # go_mesu_npz(
    #     f'/store/cartera/2d_monolayer/hydro_phi0.114_L320_seed1_t1h_0.5.npz',
    #     suffix=f'_seed1_t1h_0.5',
    #     skip_rsync=True
    # )

    # vacf
    go_mesu_npz(
        '/store/cartera/2d_monolayer/nohydro_phi0.01_L640_t2s_0.0005.npz',
        suffix=f'_t2s_0.0005',
        # skip_rsync=True
    )