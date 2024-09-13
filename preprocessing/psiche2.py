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
}

skips = ['0154', '0015', '0016', '0046',
          '0003', '0004', '0005', '0006', '0007', '0008', '0009',
          '0010', '0011', '0012', '0013',
          '0166', # tomo
          '0164', # abandoned
          ]

def preprocess(directory_path, directory_name, destination_filename, destination_desc):
    print('processing', directory_name)
    t0 = time.time()
    
    raw_file = f'{directory_name}.nxs'



    for skip_start in skips:
        if directory_name.startswith(skip_start):
            print(f'skipping {directory_name}')
            return

    if 'texp100ms' in directory_name:
        time_step = 0.1
    elif 'exp500ms' in directory_name or 'exptime500m' in directory_name:
        time_step = 0.5
    elif 'exp200ms' in directory_name:
        time_step = 0.2
    elif 'exp1s' in directory_name or 'exptime1' in directory_name:
        time_step = 1
    elif 'exp2s' in directory_name or 'expTime2s' in directory_name or 'exptime2s' in directory_name:
        time_step = 2
    elif 'exp5s' in directory_name or 'expTime5s' in directory_name or 'exptime5s' in directory_name:
        time_step = 5
    elif 'texp16s' in directory_name:
        time_step = 16
    elif directory_name in timestep_map:
        time_step = timestep_map[directory_name]
    else:
        time_step = float(input('I could not find the timestep. Please enter it (in seconds): '))
    print(f'time step = {time_step} + 0.012 s')
    time_step += 0.012

    return

    # load raw
    with h5py.File(f'{directory_path}/{raw_file}', 'r') as f:
        obj = f['flyscan_00001/scan_data/orca_image']
        print(obj.shape, obj.dtype)

        if obj.shape[0] > 1000 and False:
            proj = np.full((obj.shape[0]//2, obj.shape[1], obj.shape[2]), np.nan, dtype=np.float16)
            print(proj.shape, proj.dtype, proj.nbytes/1e9)

            for i in tqdm.trange(0, obj.shape[0]//2):
                proj[i, :, :] = obj[i*2, :, :]

            time_step *= 2

        else:
            proj = np.full(obj.shape, np.nan, dtype=np.float16)
            for i in tqdm.trange(obj.shape[0]): # possibly because obj is not a real numpy object
                proj[i, :, :] = obj[i, :, :]

            # proj = obj.astype(np.float16)
            print(proj.shape, proj.dtype, f'{proj.nbytes/1e9:.1f}GB')

    # load refs and dark
    with h5py.File(f'{directory_path}/pre_ref.nxs', 'r') as f:
        refA = f['ref_2d'][:]

    with h5py.File(f'{directory_path}/post_ref.nxs', 'r') as f:
        refB = f['ref_2d'][:]

    with h5py.File(f'{directory_path}/post_dark.nxs', 'r') as f:
        dark = f['dark_2d'][:]


    print('subtracting dark 1')
    refA -= dark
    print('subtracting dark 2')
    refB -= dark
    print('subtracting dark 3')
    proj -= dark
    del dark


    print('changing dtypes')
    refA = refA.astype(np.uint32)
    refB = refB.astype(np.uint32)

    print('averaging ref')
    ref = (refA + refB)/2
    del refA, refB
    
    print('checking no negative')
    assert np.all(ref > 0)

    print('dividing by ref')
    # stack = proj / ref
    for i in range(proj.shape[0]):
        proj[i] = proj[i] / ref
    del ref

    # print('getting stats')
    # print(proj.dtype, proj.shape, proj.min(), proj.mean(), proj.max())

    # print('finishing')
    # proj = proj.astype(np.float16)

    print('up is down')
    proj = proj[:, ::-1, :]

    common.save_data(f'preprocessing/data/stack_{destination_filename}.npz',
        stack=proj, pixel_size=0.325, time_step=time_step,
        NAME=destination_desc)

    n = proj.shape[0] // 50
    common.save_data(f'preprocessing/data/stack_{destination_filename}_small.npz',
        stack=proj[::n], pixel_size=0.325, time_step=time_step*n, nth_frame=n,
        NAME=destination_desc)

    print(f'done in {time.time()-t0:.0f}s')

for f in os.scandir('/data2/acarter/psiche/PSICHE_0624'):
    if f.is_dir():
        print(f.name, f.path)

        if 'tomo' in f.name.lower():
            continue

        internal_name = f'psiche{f.name[1:4]}'
        internal_desc = f.name[5:].replace('_', ' ')
        print(internal_desc)

        preprocess(f.path, f.name, internal_name, internal_desc)
        # break