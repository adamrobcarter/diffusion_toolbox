import common
import numpy as np
import warnings
import tqdm
import numpy.ma

MAX_TIME_ORIGINS = 100000

if __name__ == '__main__':
    for file in common.files_from_argv('box_counting/data/', 'counted_counts_'):
        data = common.load(f'box_counting/data/counted_counts_{file}.npz')
        counts = data['counts']
        print(counts.shape)
        num_timesteps = counts.shape[3]
        num_box_sizes = counts.shape[0]

        d_frames = common.exponential_integers(1, num_timesteps, 70) - 1
        assert np.all(d_frames < num_timesteps)

        use_every_nth_frame = max(int(num_timesteps / MAX_TIME_ORIGINS), 1)
        if use_every_nth_frame > 1:
            warnings.warn(f'Using every {use_every_nth_frame}th frame as a time origin. Eventually you may want to use every frame')
        else:
            print('using every frame as a time origin')
        time_origins = np.arange(0, num_timesteps, use_every_nth_frame)

        NtN0_avgs = np.full((num_box_sizes, len(time_origins), len(d_frames)), np.nan)
        N = np.full(num_box_sizes, np.nan)
        N_sq = np.full(num_box_sizes, np.nan)

        for box_size_index in tqdm.trange(num_box_sizes):
            these_counts = counts[box_size_index, :, :, :]
            
            # couunts comes padded by nans for different box sizes, so we remove these
            # this is a horrible way of doing that, but i cba
            if np.isnan(these_counts[0, :, 0]).any():
                num_boxes_y = np.argmax(np.isnan(these_counts[0, :, 0]))
                these_counts = these_counts[:, :num_boxes_y, :]
            if np.isnan(these_counts[:, 0, 0]).any():
                num_boxes_x = np.argmax(np.isnan(these_counts[:, 0, 0]))
                these_counts = these_counts[:num_boxes_x, :, :]
                
            assert these_counts.size
            assert np.isnan(these_counts).sum() == 0

            for time_origin_index, time_origin in enumerate(time_origins):
                N0 = these_counts[:, :, time_origin]
                used_t = time_origin + d_frames
                used_t = used_t[used_t < num_timesteps]
                Nt = these_counts[:, :, used_t]

                NtN0 = Nt[:, :, :] * N0[:, :, np.newaxis]
                assert np.isnan(NtN0).sum() == 0
                assert NtN0.size

                NtN0_avg = NtN0.mean(axis=(0, 1)) # avg over boxes in x and y
                assert np.isnan(NtN0_avg).sum() == 0

                NtN0_avgs[box_size_index, time_origin_index, :NtN0_avg.size] = NtN0_avg

            N[box_size_index] = these_counts.mean()
            N_sq[box_size_index] = (these_counts**2).mean()


        final_NtN0 = np.nanmean(NtN0_avgs, axis=1) # average over time origins
        # final_NtN0_unc = np.nanmean(NtN0_avgs, axis=1) # average over time origins
        # print(final_NtN0)

        common.save_data(f'box_counting/data/NtN0_{file}.npz',
            NtN0=final_NtN0, t=d_frames*data['time_step'], box_sizes=data['box_sizes'],
            N=N, N_sq=N_sq,
            particle_diameter=data.get('particle_diameter'),
            pixel_size=data.get('pixel_size'), window_size_x=data['window_size_x'], window_size_y=data['window_size_y'],
            NAME=data.get('NAME'), channel=data.get('channel'),
        )



