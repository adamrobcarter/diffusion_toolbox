import common
import os
import numpy as np
import tqdm
import matplotlib.pyplot as plt
import preprocessing.stack_movie

def load_tiffs_from_folder(folder):

    files = os.listdir(folder)

    stack = np.full((len(files), 2048, 2048), np.nan)

    for file_i, file in enumerate(tqdm.tqdm(files, desc='loading tiffs')):
        tiff = common.load_tif(f'{folder}/{file}')

        stack[file_i, :, :] = tiff

    return stack


flat_folder   = r'C:\Users\Adam\Documents\faxtor\SAMPLE_0003\MEASUREMENT_0007\PCO_EDGE'
dark_folder   = r'C:\Users\Adam\Documents\faxtor\SAMPLE_0003\MEASUREMENT_0006\PCO_EDGE'
sample_folder = r'C:\Users\Adam\Documents\faxtor\SAMPLE_0003\MEASUREMENT_0008\PCO_EDGE'

flats = load_tiffs_from_folder(flat_folder)
darks = load_tiffs_from_folder(dark_folder)
samples = load_tiffs_from_folder(sample_folder)

print('making figure')
fig, ax = plt.subplots(1, 1)
ax.plot(flats  .mean(axis=(1, 2)), label='flats')
ax.plot(darks  .mean(axis=(1, 2)), label='darks')
ax.plot(samples.mean(axis=(1, 2)), label='samples')
ax.legend()
common.save_fig(fig, f'preprocessing/figures_png/intensity_time_faxtor_test.png')

# preprocessing.stack_movie.save_array_movie(flats  [:20, :, :], 1, 1, 'faxtor_test', 'preprocessing/figures_png/movie_faxtor_test_flats.gif', 2048, 2048)
# preprocessing.stack_movie.save_array_movie(samples[:20, :, :], 1, 1, 'faxtor_test', 'preprocessing/figures_png/movie_faxtor_test_samples.gif', 2048, 2048)

flats = flats[1:]
darks = darks[1:]
samples = samples[1:]

print('doing processing')
flat = flats.mean(axis=0)
dark = darks.mean(axis=0)

print(f'flat mean {flats.mean()} {flats.dtype}')
print(f'dark mean {darks.mean()} {darks.dtype}')
print(f'sample mean {samples.mean()} {samples.dtype}')

print(f'flat mean  [:100, :100] {flats  [:, :100, :100].mean()} {flats.dtype}')
print(f'dark mean  [:100, :100] {darks  [:, :100, :100].mean()} {darks.dtype}')
print(f'sample mean[:100, :100] {samples[:, :100, :100].mean()} {samples.dtype}')

stack = np.log( (samples - dark) / (flat - dark) )

# stack = np.log(samples / flat)

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4))

vmin = min(flats[0, :, :].min(), samples[0, :, :].min())
vmax = max(flats[0, :, :].max(), samples[0, :, :].max())
ax1.imshow(flats[0, :, :], vmin=vmin, vmax=vmax)
ax2.imshow(samples[0, :, :], vmin=vmin, vmax=vmax)
ax3.imshow(samples[0, :, :]/flats[0, :, :])
common.save_fig(fig, f'preprocessing/figures_png/flat_sample_faxtor_test.png')

# stack = np.log(flats)

# preprocessing.stack_movie.save_array_movie(stack[:20, :, :], 1, 1, 'faxtor_test', 'preprocessing/figures_png/movie_faxtor_test_stack.gif', 2048, 2048)

stack_new = np.full((4, 2048, 2048), np.nan)
stack_new [0, :, :] = flats[0, :, :]
stack_new [1, :, :] = flats[-1, :, :]
stack_new [2, :, :] = samples[0, :, :]
stack_new [3, :, :] = samples[-1, :, :]
preprocessing.stack_movie.save_array_movie(stack_new[:, :, :], 1, 1, 'faxtor_test', 'preprocessing/figures_png/movie_faxtor_test_flat_sample.gif', 2048, 2048)


common.save_data('preprocessing/data/stack_faxtor_test.npz', stack=stack, pixel_size=1, time_step=1)