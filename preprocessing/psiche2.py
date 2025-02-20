import sys
# import preprocessing.edf_read_dirty
import numpy as np
import common
# import matplotlib.pyplot as plt
import h5py
import os
# import gc
import tqdm
import time

timestep_map = {
    '0014_vermiculite_1_2um-dil10': 0.5,
    '0165_Muscovite50_100um_silice4um': 0.2,
    '0073_Muscovite100_200um_radios2': 0.5,
    '0045_muscovite_100_200um': 0.5,
    '0185_RCPSi125': 1,
    '0022_vermiculite_1_2um': 1,
    '0044_muscovite_100_200um': 0.5,
    '0187_RCPSi125_musc2-10': 1, # guessed!!!!!!!
    '0048_muscovite_100_200um_silice4um_interface': 0.5,
    '0057_Glass106um_Radios2': 0.5, # guessed !!!!!!!
    '0189_H2O_muscovite2-10um_fall_200mn': 0.2,
    '0047_muscovite_100-200um_time15min-silice4um': 0.5,
    '0162_Muscovite50_100um': 1,
    '0020_silice7p75um_sediment_x6p8_z15': 0.5,
    '0049_muscovite__100_200um_silice4um_long': 0.5,
    '0023_silice4um': 1,
    '0021_silice1p7um': 1,
    '0080_muscovite100-200um_silice1p7um_NaOH': 2,
    '0101_RCP_PMMA_KI0_1M': 0.1,
    '0163_Muscovite50_100um_silice4um': 1,
    '0092_muscovite50_100um_silice4um_z19apres': 2, # guessed!!!
    '0024_silice4um_sediment': 1,
    '0011_H2O_glass_cap': 1, # guessed
    '0010_emptyglasscap': 1, # guessed
    '0009_test_flat_long': 1, # guessed
}

endframe_map = {
    'psiche074': 96,
}

skips = ['0154', '0015', '0016', '0046',
          '0003', '0004', '0005', '0006', '0007', '0008',
        #   '0009',
        #   '0010',
        #   '0011',
          '0012', '0013',
          '0018', # no timestep recorded
          '0019', '0091', # bin
          '0166', '0092', # tomo
          '0164', # abandoned
          ]

def preprocess(directory_path, directory_name, destination_filename, destination_desc):
    print(f'processing {directory_name}')
    t0 = time.time()
    
    raw_file = f'{directory_name}.nxs'



    for skip_start in skips:
        if directory_name.startswith(skip_start):
            print(f'skipping {directory_name}')
            return

    if 'texp100ms' in directory_name:
        time_step = 0.1
    elif 'exp0p2s' in directory_name or '200ms' in directory_name:
        time_step = 0.2
    elif 'exp500ms' in directory_name or 'exptime500m' in directory_name.lower():
        time_step = 0.5
    elif 'exp200ms' in directory_name:
        time_step = 0.2
    elif 'exp1s' in directory_name or 'exptime1' in directory_name or 'Time1s' in directory_name:
        time_step = 1
    elif 'exp2s' in directory_name or 'expTime2s' in directory_name or 'exptime2s' in directory_name:
        time_step = 2
    elif 'exp4s' in directory_name:
        time_step = 4
    elif 'exp5s' in directory_name or 'expTime5s' in directory_name or 'exptime5s' in directory_name:
        time_step = 5
    elif 'texp8s' in directory_name:
        time_step = 8
    elif 'texp16s' in directory_name:
        time_step = 16
    elif directory_name in timestep_map:
        time_step = timestep_map[directory_name]
    else:
        time_step = float(input('I could not find the timestep. Please enter it (in seconds): '))
    # print(f'time step = {time_step} + 0.012 s')
    time_step += 0.012

    # load raw
    with h5py.File(f'{directory_path}/{raw_file}', 'r') as f:
        obj = f['flyscan_00001/scan_data/orca_image']
        # print(obj.shape, obj.dtype)

        if obj.shape[0] > 1000 and False:
            proj = np.full((obj.shape[0]//2, obj.shape[1], obj.shape[2]), np.nan, dtype=np.float16)
            print(proj.shape, proj.dtype, proj.nbytes/1e9)

            for i in tqdm.trange(0, obj.shape[0]//2):
                proj[i, :, :] = obj[i*2, :, :]

            time_step *= 2

        else:
            final_frame = obj.shape[0]
            if destination_filename == 'psiche089':
                final_frame = 531
            proj = np.full((final_frame, obj.shape[1], obj.shape[2]), np.nan, dtype=np.float16)
            
            for i in range(final_frame): # possibly because obj is not a real numpy object
                proj[i, :, :] = obj[i, :, :]

            # proj = obj.astype(np.float16)
            # print(proj.shape, proj.dtype, f'{proj.nbytes/1e9:.1f}GB')

        del obj

    if end := endframe_map.get(destination_filename):
        proj = proj[:end+1, :, :]
        print('final frame', end)

    # load refs and dark
    with h5py.File(f'{directory_path}/pre_ref.nxs', 'r') as f:
        refA = f['ref_2d'][:] # flat field - no sample, beam on

    with h5py.File(f'{directory_path}/post_ref.nxs', 'r') as f:
        refB = f['ref_2d'][:] # flat field - no sample, beam on

    with h5py.File(f'{directory_path}/post_dark.nxs', 'r') as f:
        dark = f['dark_2d'][:] # dark - beam off


    # print('subtracting dark 1')
    refA -= dark
    # print('subtracting dark 2')
    refB -= dark
    # print('subtracting dark 3')
    proj -= dark
    del dark


    # print('changing dtypes')
    refA = refA.astype(np.uint32)
    refB = refB.astype(np.uint32)

    # print('averaging ref')
    ref = (refA + refB)/2
    if destination_filename == 'psiche089':
        ref = refA
    del refA, refB
    
    # print('checking no negative')
    assert np.all(ref > 0)

    # print('dividing by ref')
    # stack = proj / ref
    for i in range(proj.shape[0]):
        proj[i] = proj[i] / ref
    del ref

    # print('getting stats')
    # print(proj.dtype, proj.shape, proj.min(), proj.mean(), proj.max())

    # print('finishing')
    # proj = proj.astype(np.float16)

    # print('up is down')
    proj = proj[:, ::-1, :]

    if destination_filename in ['psiche086', 'psiche087', 'psiche088']:
        proj = proj[:95, :, :]
        print('removed later frames')

    particle_diameter = None
    if 'silice4' in directory_name:
        particle_diameter = 4
    if 'silice1p7' in directory_name:
        particle_diameter = 1.7
    if 'silice7p75' in directory_name:
        particle_diameter = 7.75
    if 'silice0p96' in directory_name:
        particle_diameter = 0.96
    if 'silice0p7' in directory_name:
        particle_diameter = 0.7
    if 'PS2' in directory_name:
        particle_diameter = 2
    if 'PS1' in directory_name:
        particle_diameter = 1

    to_save = dict(
        NAME=destination_desc,
        particle_diameter=particle_diameter
    )

    common.save_data(f'preprocessing/data/stack_{destination_filename}.npz',
        stack=proj, pixel_size=0.325, time_step=time_step,
        **to_save)

    n = proj.shape[0] // 50
    common.save_data(f'preprocessing/data/stack_{destination_filename}_small.npz',
        stack=proj[::n], pixel_size=0.325, time_step=time_step*n, nth_frame=n,
        **to_save)

    # print(f'done in {time.time()-t0:.0f}s')


files = list(os.scandir('/data2/acarter/psiche/PSICHE_0624'))

do = ['psiche029', 'psiche030', 'psiche031', 'psiche032', 'psiche033', 'psiche034',
      'psiche035', 'psiche036', 'psiche037', 'psiche038', 'psiche039']
# do = []

if len(do):
    print('WARNING NOT DOING ALL')

for f in tqdm.tqdm(files):
    try:
        if f.is_dir():
            # print(f.name, f.path)

            if 'tomo' in f.name.lower():
                # print('skipping tomo')
                continue

            if f.name.startswith('_') or f.name == 'slurmdir':
                continue

            internal_name = f'psiche{f.name[1:4]}'
            internal_desc = f.name[5:].replace('_', ' ')
            # print(internal_desc)

            if (len(do)) and not (internal_name in do):
                continue

            do.remove(internal_name)

            preprocess(f.path, f.name, internal_name, internal_desc)
            # f.path: /data2/acarter/psiche/PSICHE_0624/0064_PS3um_Au_exp500ms
            # f.name: 0064_PS3um_Au_exp500ms
            # internal_name: psiche064
            # internal_desc: PS3um Au exp500ms
            
            # break
    except Exception as err:
        print(f'failed on {f.name}')
        print(err)

if len(do):
    left = ' '.join(do)
    raise Exception(f'Did not find all items: {left}')