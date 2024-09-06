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

for file in sys.argv[1:]:
    print('processing', file)
    t0 = time.time()

    raw_file = None
    for filename in os.listdir(f'raw_data/psiche/{file}/'):
        if filename.endswith('.nxs') and filename not in ['post_dark.nxs', 'pre_ref.nxs', 'post_ref.nxs']:
            print('found', filename)
            raw_file = filename

            # get name
            with open(f'raw_data/psiche/{file}/NAME', 'w') as namefile:
                file_label = filename[5:-4]
                namefile.write(file_label)



    if 'exp500ms' in file_label:
        time_step = 0.5
    elif 'exp200ms' in file_label:
        time_step = 0.2
    elif 'exp1s' in file_label:
        time_step = 1
    elif 'exp2s' in file_label:
        time_step = 2
    elif 'exp5s' in file_label:
        time_step = 5
    else:
        time_step = float(input('I could not find the timestep. Please enter it (in seconds): '))
    print(f'time step = {time_step} + 0.012 s')
    time_step += 0.012

    # load raw
    with h5py.File(f'raw_data/psiche/{file}/{raw_file}', 'r') as f:
        obj = f['flyscan_00001/scan_data/orca_image']
        print(obj.shape)

        if obj.shape[0] > 1000:
            proj = np.full((obj.shape[0]//2, obj.shape[1], obj.shape[2]), np.nan, dtype=np.float16)
            print(proj.shape, proj.dtype, proj.nbytes/1e9)

            for i in tqdm.trange(0, obj.shape[0]//2):
                proj[i, :, :] = obj[i*2, :, :]

            time_step *= 2

        else:
            proj = np.full(obj.shape, np.nan, dtype=np.float16)
            print(proj.shape, proj.dtype, proj.nbytes/1e9)

            for i in tqdm.trange(obj.shape[0]):
                proj[i, :, :] = obj[i, :, :]

    # load refs and dark
    with h5py.File(f'raw_data/psiche/{file}/pre_ref.nxs', 'r') as f:
        refA = f['ref_2d'][:]

    with h5py.File(f'raw_data/psiche/{file}/post_ref.nxs', 'r') as f:
        refB = f['ref_2d'][:]

    with h5py.File(f'raw_data/psiche/{file}/post_dark.nxs', 'r') as f:
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

    common.save_data(f'preprocessing/data/stack_{file}.npz', stack=proj, pixel_size=0.325, time_step=time_step, we_processed=True)

    n = proj.shape[0] // 50
    common.save_data(f'preprocessing/data/stack_{file}_small.npz', stack=proj[::n], pixel_size=0.325, time_step=time_step*n, we_processed=True, nth_frame=n)

    print(f'done in {time.time()-t0:.0f}s')