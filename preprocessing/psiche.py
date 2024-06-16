import sys
import preprocessing.edf_read_dirty
import numpy as np
import common
import h5py
import os
import sys
import DDM.static_fourier
"""
for file in sys.argv[1:]:
    print('processing', file)

    for filename in os.listdir(f'raw_data/psiche/{file}/'):
        if filename.endswith('.nxs') and filename not in ['post_dark.nxs', 'pre_ref.nxs', 'post_ref.nxs']:

            # with h5py.File(f'raw_data/psiche/{file}/{filename}', 'r') as f:
            #     scan = f['flyscan_00001']

            with open(f'raw_data/psiche/{file}/NAME', 'w') as namefile:
                file_label = filename[5:-4]
                namefile.write(file_label)
                    
                # filecontents = f['flyscan_00001']
                # scan_data = filecontents['scan_data']
                # raw_stack = scan_data['orca_image'][:]
                # print(raw_stack.shape)
            
    if file_label.endswith('500ms'):
        time_step = 0.5
    elif file_label.endswith('200ms'):
        time_step = 0.2
    elif file_label.endswith('1s'):
        time_step = 1
    elif file_label.endswith('2s'):
        time_step = 2
    elif file_label.endswith('5s'):
        time_step = 5
    else:
        time_step = float(input('I could not find the timestep. Please enter it (in seconds): '))
    print(f'time step = {time_step} + 0.012 s')
    time_step += 0.012

    data = preprocessing.edf_read_dirty.vol_read_virtual(f'raw_data/psiche/{file}/cor.vol')
    print('converting to float16')
    data = data.astype(np.float16)

    print('transposing')
    stack = np.transpose(data, axes=(2, 1, 0)) # x,y,t to t,y,x

    stack = stack[:, ::-1, :]

    pixel_size = 0.325


    # print('>1 fraction',   (stack > 1) .sum()/stack.shape[0])
    # print('zero fraction', (stack == 0).sum()/stack.size)
    # print('nan per frame', np.isnan(stack).sum()/stack.shape[0])
    # print('inf per frame', np.isinf(stack).sum()/stack.shape[0])
    print('min, mean, max', stack.min(), stack.mean(), stack.max())

    common.save_data(f'preprocessing/data/stack_{file}.npz', stack=stack, pixel_size=pixel_size, time_step=time_step)

    DDM.static_fourier.do_static_fourier(file, stack, pixel_size)

    # bins = np.logspace(np.log10(stack.min()), np.log10(stack.max()), 200)
    # bins = np.linspace(stack.min(), stack.max(), 200)
    # fig_hist, ax_hist = plt.subplots(1, 1)
    # ax_hist.hist(stack.flatten(), bins=bins)
    # del stack
    # # ax_hist.loglog()
    # ax_hist.semilogy()
    # common.save_fig(fig_hist, f'preprocessing/figures_png/hist_{file}.png')

    # stack = stack[:, 100:-100, 100:-100]
    # print('cropping edge')

    # stack = np.copy(stack)
    # stack[stack > 1] = 1
    # print('removed >1')
    

print('done')