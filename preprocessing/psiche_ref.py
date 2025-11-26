import sys
import preprocessing.edf_read_dirty
import numpy as np
import common
import matplotlib.pyplot as plt
import h5py

if __name__ == '__main__':
    for file in sys.argv[1:]:
        print('processing', file)

        # with h5py.File(f'raw_data/psiche/{file}/raw.nxs', 'r') as f:
        #     filecontents = f['flyscan_00001']
        #     scan_data = filecontents['scan_data']
        #     del filecontents
        #     raw_stack = scan_data['orca_image'][:]
        #     del scan_data
        #     print(raw_stack.nbytes / 1e9)

            # print(a['scan_log'])
            # print(filecontents['title'][()])
            # print(filecontents['configuration'][()])
            # print(filecontents['experiment_identifier'][()])
            # print(a['start_time'])
            # print(a['end_time'])
            # print(a['duration'])
            # print(a['run_cycle'])
            # print(a['User'])
            # print(a['PSICHE'])

        with h5py.File(f'raw_data/psiche/{file}/pre_ref.nxs', 'r') as f:
            refA = f['ref_0'][:]
            refA_mean = np.median(refA, axis=0)
            print(refA.nbytes / 1e9)
            del refA

        with h5py.File(f'raw_data/psiche/{file}/post_ref.nxs', 'r') as f:
            refB = f['ref_0'][:]
            refB_mean = np.median(refB, axis=0)
            print(refB.nbytes / 1e9)
            del refB

        with h5py.File(f'raw_data/psiche/{file}/post_dark.nxs', 'r') as f:
            dark0 = f['dark_0'][:]
            dark_mean = np.median(dark0, axis=0)
            print(dark0.nbytes / 1e9)
            del dark0

        ref_avg = (refA_mean + refB_mean)/2 
        stack = ref_avg - dark_mean
        # stack = (raw_stack - dark_mean) / (ref_avg - dark_mean)

        # print('setting negative values to zero')
        # stack[stack < 0] = 0

        print(stack.dtype, stack.shape, stack.min(), stack.mean(), stack.max())



        bins = np.logspace(np.log10(stack.min()), np.log10(stack.max()), 200)
        bins = np.linspace(stack.min(), stack.max(), 200)
        fig_hist, ax_hist = plt.subplots(1, 1)
        ax_hist.hist(stack.flatten(), bins=bins)
        # ax_hist.loglog()
        ax_hist.semilogy()
        common.save_fig(fig_hist, f'preprocessing/figures_png/hist2_{file}.png')

        stack = dark_mean
        bins = np.logspace(np.log10(stack.min()), np.log10(stack.max()), 200)
        bins = np.linspace(stack.min(), stack.max(), 200)
        fig_hist, ax_hist = plt.subplots(1, 1)
        ax_hist.hist(stack.flatten(), bins=bins)
        # ax_hist.loglog()
        ax_hist.semilogy()
        common.save_fig(fig_hist, f'preprocessing/figures_png/hist_dark_{file}.png')

        stack = ref_avg
        bins = np.logspace(np.log10(stack.min()), np.log10(stack.max()), 200)
        bins = np.linspace(stack.min(), stack.max(), 200)
        fig_hist, ax_hist = plt.subplots(1, 1)
        ax_hist.hist(stack.flatten(), bins=bins)
        # ax_hist.loglog()
        ax_hist.semilogy()
        common.save_fig(fig_hist, f'preprocessing/figures_png/hist_ref_{file}.png')

        # stack = raw_stack
        # bins = np.logspace(np.log10(stack.min()), np.log10(stack.max()), 200)
        # bins = np.linspace(stack.min(), stack.max(), 200)
        # fig_hist, ax_hist = plt.subplots(1, 1)
        # ax_hist.hist(stack.flatten(), bins=bins)
        # # ax_hist.loglog()
        # ax_hist.semilogy()
        # common.save_fig(fig_hist, f'preprocessing/figures_png/hist_stack_raw_{file}.png')

        # with h5py.File(f'raw_data/psiche/{file}/post_ref.nxs', 'r') as f:
        #     filecontents = f['flyscan_00001']
        #     data = filecontents['scan_data']
        #     raw_stack = data['orca_image'][:]
        #     print(raw_stack.shape)

        # with h5py.File(f'raw_data/psiche/{file}/post_dark.nxs', 'r') as f:
        #     filecontents = f['flyscan_00001']
        #     data = filecontents['scan_data']
        #     raw_stack = data['orca_image'][:]
        #     print(raw_stack.shape)


        # data = preprocessing.edf_read_dirty.vol_read_virtual(f'raw_data/psiche/{file}/cor.vol')
        # stack = np.transpose(data, axes=(2, 0, 1))

        # print('>1 fraction',   (stack > 1) .sum()/stack.shape[0])
        # print('zero fraction', (stack == 0).sum()/stack.size)

        # print('min, mean, max', stack.min(), stack.mean(), stack.max())

        # bins = np.logspace(np.log10(stack.min()), np.log10(stack.max()), 200)
        # bins = np.linspace(stack.min(), stack.max(), 200)
        # fig_hist, ax_hist = plt.subplots(1, 1)
        # ax_hist.hist(stack.flatten(), bins=bins)
        # # ax_hist.loglog()
        # common.save_fig(fig_hist, f'preprocessing/figures_png/hist_{file}.png')

        # # stack = stack[:, 100:-100, 100:-100]
        # # print('cropping edge')

        # # stack = np.copy(stack)
        # # stack[stack > 1] = 1
        # # print('removed >1')
        
        common.save_data(f'preprocessing/data/stack_{file}_2.npz', stack=stack, 
                pixel_size=0.325, time_step=1.06)

    print('done')