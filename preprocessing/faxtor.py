import common
import os
import numpy as np
import tqdm
import matplotlib.pyplot as plt
import preprocessing.stack_movie
import scipy.ndimage
import matplotlib.cm

DTYPE = np.float16

#RAW_DATA_PATH = /data2/acarter/faxtor
RAW_DATA_PATH = r'D:\RAW'

def load_tiffs_from_folder(folder, every_nth=1, fraction_of=1, trim_start=0):
    assert fraction_of >= 1
    assert every_nth >= 1

    files = os.listdir(folder)

    files = sorted(files)
    files = files[trim_start:]
    files = files[:len(files)//fraction_of:every_nth]
    print(f'{len(files)} files')

    stack = np.full((len(files), 2048, 2048), np.nan, dtype=DTYPE)
    print('stack arraysize', common.arraysize(stack))

    for file_i, file in enumerate(tqdm.tqdm(files, desc='loading tiffs')):
        tiff = common.load_tif(f'{folder}/{file}').astype(DTYPE)
        assert len(tiff.shape) == 2, f'{file} has shape {tiff.shape}'

        # flip the ting
        # tiff = tiff[::-1, :]
        tiff = np.swapaxes(tiff, 1, 0)

        stack[file_i, :, :] = tiff

    return stack

def go(flat_folder, dark_folder, sample_folder, name, exposure_time, frame_gap, every_nth=1, fraction_of=1, trim_start=0, **kwargs):
    time_step = (exposure_time + frame_gap) * every_nth

    fig, ax = plt.subplots(1, 1)

    if flat_folder:
        flats = load_tiffs_from_folder(flat_folder, every_nth, fraction_of)
        ax.plot(flats.mean(axis=(1, 2)), label='flats')

    if dark_folder:
        darks = load_tiffs_from_folder(dark_folder, every_nth, fraction_of)
        ax.plot(darks.mean(axis=(1, 2)), label='darks')

    samples = load_tiffs_from_folder(sample_folder, every_nth, fraction_of, trim_start)
    ax.plot(samples.mean(axis=(1, 2)), label='samples')

    ax.legend()
    ax.set_xlabel('frame')
    ax.set_ylabel('average intensity')
    common.save_fig(fig, f'preprocessing/figures_png/intensity_time_{name}.png', dpi=300)

    plt.close(fig) # prevent memory warning/leak

    # preprocessing.stack_movie.save_array_movie(flats  [:20, :, :], 1, 1, 'faxtor_test', 'preprocessing/figures_png/movie_faxtor_test_flats.gif', 2048, 2048)
    # preprocessing.stack_movie.save_array_movie(samples[:20, :, :], 1, 1, 'faxtor_test', 'preprocessing/figures_png/movie_faxtor_test_samples.gif', 2048, 2048)

    print('doing processing')

    samples = samples[1:] # the first frame seems to be off sometimes
    print(f'sample mean {samples.mean()} {samples.dtype}')

    if dark_folder:
        darks = darks[1:]
        dark = darks.mean(axis=0)
        print(f'dark mean {darks.mean()} {darks.dtype}')
        assert np.all(samples > dark)
    
    if flat_folder:
        flats = flats[1:]
        flat = flats.mean(axis=0)
        print(f'flat mean {flats.mean()} {flats.dtype}')
        assert np.all(flat > dark)

    # print(f'flat mean  [:100, :100] {flats  [:, :100, :100].mean()} {flats.dtype}')
    # print(f'dark mean  [:100, :100] {darks  [:, :100, :100].mean()} {darks.dtype}')
    # print(f'sample mean[:100, :100] {samples[:, :100, :100].mean()} {samples.dtype}')


    if flat_folder and dark_folder:
        stack = np.log( (samples - dark) / (flat - dark) )
    elif flat_folder:
        stack = np.log( samples / flat )
    elif dark_folder:
        stack = np.log( samples - dark )
    else:
        stack = np.log( samples )


    # stack = np.log(samples / flat)

    # fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4))

    # vmin = min(flats[0, :, :].min(), samples[0, :, :].min())
    # vmax = max(flats[0, :, :].max(), samples[0, :, :].max())
    # ax1.imshow(flats[0, :, :], vmin=vmin, vmax=vmax)
    # ax2.imshow(samples[0, :, :], vmin=vmin, vmax=vmax)
    # ax3.imshow(samples[0, :, :]/flats[0, :, :])
    # common.save_fig(fig, f'preprocessing/figures_png/flat_sample_{name}.png')

    # stack = np.log(flats)

    # preprocessing.stack_movie.save_array_movie(stack[:20, :, :], 1, time_step, name, f'preprocessing/figures_png/movie_{name}_stack.gif', 2048, 2048)
    # preprocessing.stack_movie.save_array_movie(stack[:, :, :], 1, time_step, name, f'preprocessing/figures_png/movie_{name}_stack_highlights.gif', 2048, 2048, hightlights=True)

    # stack_new = np.full((4, 2048, 2048), np.nan, dtype=np.float32)
    # stack_new [0, :, :] = flats[0, :, :]
    # stack_new [1, :, :] = flats[-1, :, :]
    # stack_new [2, :, :] = samples[0, :, :]
    # stack_new [3, :, :] = samples[-1, :, :]
    # preprocessing.stack_movie.save_array_movie(stack_new[:, :, :], 1, 2, name, f'preprocessing/figures_png/movie_{name}_flat_sample.gif', 2048, 2048)

    print(common.nanfrac(stack), 'nanfrac')
    assert np.isfinite(stack).all()

    common.save_data(f'preprocessing/data/stack_{name}.npz',       stack=stack,      pixel_size=0.65, time_step=time_step, **kwargs)
    every = max(stack.shape[0]//50, 1)
    common.save_data(f'preprocessing/data/stack_{name}_small.npz', stack=stack[::every], pixel_size=0.65, time_step=time_step*every, **kwargs)

    print()

# def go_filter(sample_folder, name, exposure_time, frame_gap, **kwargs):
#     time_step = exposure_time + frame_gap

#     fig, ax = plt.subplots(1, 1)

#     samples = load_tiffs_from_folder(sample_folder)
#     ax.plot(samples.mean(axis=(1, 2)), label='samples')

#     ax.legend()
#     ax.set_xlabel('frame')
#     ax.set_ylabel('average intensity')
#     common.save_fig(fig, f'preprocessing/figures_png/intensity_time_{name}.png', dpi=300)

#     print('doing processing')

#     samples = samples[1:]

#     stack = np.log( samples )
#     print('doing filter')

#     fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))

#     stack = stack - stack.mean(axis=0)

#     stack_lowpass = np.full_like(stack, np.nan)
#     for frame in range(stack.shape[0]):
#         stack_lowpass[frame] = scipy.ndimage.gaussian_filter(stack[frame], 3)

#     stack_highpass = stack - stack_lowpass
#     # https://stackoverflow.com/a/6117387/1103752

#     s = stack.shape[1]
#     vmin = stack[0, s//3:2*s//3, s//3:2*s//3].min()
#     vmax = stack[0, s//3:2*s//3, s//3:2*s//3].max()

#     ax1.imshow(stack[0, :, :], vmin=vmin, vmax=vmax, cmap=matplotlib.cm.Greys)
#     ax2.imshow(stack_highpass[0, :, :], vmin=vmin, vmax=vmax, cmap=matplotlib.cm.Greys)
#     common.save_fig(fig, f'preprocessing/figures_png/highpass_{name}.png', dpi=300)


#     assert np.isfinite(stack_highpass).all()
#     common.save_data(f'preprocessing/data/stack_{name}_filter.npz',       stack=stack_highpass,      pixel_size=0.65, time_step=time_step, **kwargs)
#     common.save_data(f'preprocessing/data/stack_{name}_filter_small.npz', stack=stack_highpass[:50], pixel_size=0.65, time_step=time_step, **kwargs)

# go(
#     flat_folder   = '/data2/acarter/faxtor/SAMPLE_0003/MEASUREMENT_0007/PCO_EDGE',
#     dark_folder   = '/data2/acarter/faxtor/SAMPLE_0003/MEASUREMENT_0006/PCO_EDGE',
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0003/MEASUREMENT_0008/PCO_EDGE',
#     name = 'faxtor003',
# )
# go(
#     flat_folder   = '/data2/acarter/faxtor/SAMPLE_0003/MEASUREMENT_0007/PCO_EDGE',
#     dark_folder   = '/data2/acarter/faxtor/SAMPLE_0003/MEASUREMENT_0006/PCO_EDGE',
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0003/MEASUREMENT_0008/PCO_EDGE',
#     name = 'faxtor003',
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0006/EXPERIMENT_0000/MEASUREMENT_0004/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor006a',
#     exposure_time = 0.05,
#     frame_gap = 0.25,
#     particle_diameter = 7.75,
#     particle_material = 'SiO2',
#     trim_start = 5000,
#     # every_nth = 2,
#     fraction_of = 2,
#     NAME = 'Si 120um, Si 7.75um injected',
# )
# go(
#     flat_folder   = '/data2/acarter/faxtor/SAMPLE_0006/EXPERIMENT_0000/MEASUREMENT_0008/PCO_EDGE',
#     dark_folder   = '/data2/acarter/faxtor/SAMPLE_0006/EXPERIMENT_0000/MEASUREMENT_0009/PCO_EDGE',
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0006/EXPERIMENT_0000/MEASUREMENT_0007/PCO_EDGE',
#     name = 'faxtor006b',
#     exposure_time = 0.05,
#     frame_gap = 1,
#     particle_diameter = 8,
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0007/EXPERIMENT_0000/MEASUREMENT_0003/PCO_EDGE',
#     flat_folder   = '/data2/acarter/faxtor/SAMPLE_0007/EXPERIMENT_0000/MEASUREMENT_0004/PCO_EDGE',
#     dark_folder   = '/data2/acarter/faxtor/SAMPLE_0007/EXPERIMENT_0000/MEASUREMENT_0005/PCO_EDGE',
#     name = 'faxtor007',
#     exposure_time = 0.05,
#     frame_gap = 1,
#     particle_diameter = 8,
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0008/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     flat_folder   = '/data2/acarter/faxtor/SAMPLE_0008/EXPERIMENT_0000/MEASUREMENT_0002/PCO_EDGE',
#     dark_folder   = '/data2/acarter/faxtor/SAMPLE_0008/EXPERIMENT_0000/MEASUREMENT_0003/PCO_EDGE',
#     name = 'faxtor008',
#     exposure_time = 0.05,
#     frame_gap = 1,
#     particle_diameter = 8,
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0008/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor_flat',
#     exposure_time = 0.05,
#     frame_gap = 1,
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0012/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     flat_folder = '/data2/acarter/faxtor/SAMPLE_0012/EXPERIMENT_0000/MEASUREMENT_0002/PCO_EDGE',
#     dark_folder = '/data2/acarter/faxtor/SAMPLE_0012/EXPERIMENT_0000/MEASUREMENT_0004/PCO_EDGE',
#     name = 'faxtor012',
#     exposure_time = 0.05,
#     frame_gap = 1,
#     particle_diameter = 4,
#     NAME = '4um PS 0.25g/L'
# )
# go_filter(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0012/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     name = 'faxtor012',
#     exposure_time = 0.05,
#     frame_gap = 1,
#     particle_diameter = 4,
#     NAME = '4um PS 0.25g/L'
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0013/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     flat_folder = '/data2/acarter/faxtor/SAMPLE_0013/EXPERIMENT_0000/MEASUREMENT_0002/PCO_EDGE',
#     dark_folder = '/data2/acarter/faxtor/SAMPLE_0013/EXPERIMENT_0000/MEASUREMENT_0003/PCO_EDGE',
#     name = 'faxtor013',
#     exposure_time = 0.05,
#     frame_gap = 1,
#     particle_diameter = 3,
#     NAME = '3um PS 0.1g/L'
# )
# go_filter(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0013/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     name = 'faxtor013',
#     exposure_time = 0.05,
#     frame_gap = 1,
#     particle_diameter = 3,
#     NAME = '3um PS 0.1g/L'
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0015/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     flat_folder   = '/data2/acarter/faxtor/SAMPLE_0015/EXPERIMENT_0000/MEASUREMENT_0002/PCO_EDGE',
#     dark_folder   = '/data2/acarter/faxtor/SAMPLE_0015/EXPERIMENT_0000/MEASUREMENT_0003/PCO_EDGE',
#     name = 'faxtor015',
#     exposure_time = 0.05,
#     frame_gap = 1,
#     # particle_diameter = 3,
#     # NAME = '3um PS 0.1g/L'
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0013/EXPERIMENT_0000/MEASUREMENT_0004/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor_flat_long',
#     exposure_time = 0.05,
#     frame_gap = 1,
#     # particle_diameter = 3,
#     NAME = 'long flatfield',
#     every_nth = 8,
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0014/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     flat_folder   = '/data2/acarter/faxtor/SAMPLE_0014/EXPERIMENT_0000/MEASUREMENT_0002/PCO_EDGE',
#     dark_folder   = '/data2/acarter/faxtor/SAMPLE_0014/EXPERIMENT_0000/MEASUREMENT_0003/PCO_EDGE',
#     name = 'faxtor014',
#     exposure_time = 0.05,
#     frame_gap = 1,
#     particle_diameter = 3,
#     NAME = '3um PS 50g/L',
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0013/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     flat_folder   = '/data2/acarter/faxtor/SAMPLE_0013/EXPERIMENT_0000/MEASUREMENT_0002/PCO_EDGE',
#     dark_folder   = '/data2/acarter/faxtor/SAMPLE_0013/EXPERIMENT_0000/MEASUREMENT_0003/PCO_EDGE',
#     name = 'faxtor013',
#     exposure_time = 0.05,
#     frame_gap = 1,
#     particle_diameter = 3,
#     NAME = '3um PS 0.1g/L',
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0008/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     # flat_folder   = '/data2/acarter/faxtor/SAMPLE_0008/EXPERIMENT_0000/MEASUREMENT_0002/PCO_EDGE',
#     # dark_folder   = '/data2/acarter/faxtor/SAMPLE_0008/EXPERIMENT_0000/MEASUREMENT_0003/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor008',
#     exposure_time = 0.05,
#     frame_gap = 1,
#     particle_diameter = 8,
#     NAME = '8um PS 2g/L',
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0010/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     # flat_folder   = '/data2/acarter/faxtor/SAMPLE_0010/EXPERIMENT_0000/MEASUREMENT_0002/PCO_EDGE',
#     # dark_folder   = '/data2/acarter/faxtor/SAMPLE_0010/EXPERIMENT_0000/MEASUREMENT_0003/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor010',
#     exposure_time = 0.05,
#     frame_gap = 1,
#     particle_diameter = 8,
#     particle_material = 'PS',
#     NAME = '8um PS 50g/L',
# # )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0006/EXPERIMENT_0000/MEASUREMENT_0004/PCO_EDGE',
#     # flat_folder   = '/data2/acarter/faxtor/SAMPLE_0008/EXPERIMENT_0000/MEASUREMENT_0002/PCO_EDGE',
#     # dark_folder   = '/data2/acarter/faxtor/SAMPLE_0008/EXPERIMENT_0000/MEASUREMENT_0003/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor006a',
#     exposure_time = 0.05,
#     frame_gap = 0.25,
#     particle_diameter = 8,
#     particle_material = 'SiO2',
#     NAME = 'Si 120um, Si 8um injected',
#     # every_nth = 16,
#     # fraction_of = 16,
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0015/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     # flat_folder   = '/data2/acarter/faxtor/SAMPLE_0008/EXPERIMENT_0000/MEASUREMENT_0002/PCO_EDGE',
#     # dark_folder   = '/data2/acarter/faxtor/SAMPLE_0008/EXPERIMENT_0000/MEASUREMENT_0003/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor015a',
#     exposure_time = 0.05,
#     frame_gap = 1,
#     particle_material = 'vermiculite',
#     NAME = 'vermiculite 2-10um bulk',
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0015/EXPERIMENT_0000/MEASUREMENT_0004/PCO_EDGE',
#     # flat_folder   = '/data2/acarter/faxtor/SAMPLE_0008/EXPERIMENT_0000/MEASUREMENT_0002/PCO_EDGE',
#     # dark_folder   = '/data2/acarter/faxtor/SAMPLE_0008/EXPERIMENT_0000/MEASUREMENT_0003/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor015b',
#     exposure_time = 0.05,
#     frame_gap = 1,
#     particle_material = 'vermiculite',
#     NAME = 'vermiculite 2-10um bulk (bottom)',
# )



# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0016/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor016a',
#     exposure_time = 0.05,
#     frame_gap = 1,
#     particle_material = 'PS',
#     particle_diameter = 8,
#     NAME = 'Si 120um, 8um PS (1g/L)',
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0016/EXPERIMENT_0000/MEASUREMENT_0004/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor016b',
#     exposure_time = 0.05,
#     frame_gap = 1,
#     particle_material = 'PS',
#     particle_diameter = 8,
#     NAME = 'Si 120um, 8um PS (1g/L), bottom',
# )

# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0017/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor017a',
#     exposure_time = 0.05,
#     frame_gap = 1,
#     particle_material = 'PS',
#     particle_diameter = 4,
#     NAME = 'Si 120um, 4um PS (0.0625g/L)',
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0017/EXPERIMENT_0000/MEASUREMENT_0004/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor017b',
#     exposure_time = 0.05,
#     frame_gap = 1,
#     particle_material = 'PS',
#     particle_diameter = 4,
#     NAME = 'Si 120um, 4um PS (0.0625g/L), interface',
# )

# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0018/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor018',
#     exposure_time = 0.05,
#     frame_gap = 10,
#     particle_material = 'PS',
#     particle_diameter = 8,
#     NAME = 'Si 120um, 8um PS (0.5g/L)',
#     every_nth = 5,
# )

# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0019/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor019a',
#     exposure_time = 0.05,
#     frame_gap = 1,
#     particle_material = 'PS',
#     particle_diameter = 8,
#     NAME = 'vermiculite 100-200um, 8um PS (2g/L)',
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0019/EXPERIMENT_0000/MEASUREMENT_0002/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor019b',
#     exposure_time = 0.05,
#     frame_gap = 1,
#     particle_material = 'PS',
#     particle_diameter = 8,
#     NAME = 'vermiculite 100-200um, 8um PS (2g/L)',
# )

# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0020/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor020a',
#     exposure_time = 0.05,
#     frame_gap = 10,
#     particle_material = 'PS',
#     particle_diameter = 8,
#     NAME = 'muscovite 100-200um, 8um PS (1g/L)',
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0020/EXPERIMENT_0000/MEASUREMENT_0007/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor020b',
#     exposure_time = 0.05,
#     frame_gap = 1,
#     particle_material = 'PS',
#     particle_diameter = 8,
#     NAME = 'muscovite 100-200um, 8um PS (1g/L)',
# )

# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0021/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor021a',
#     exposure_time = 0.05,
#     frame_gap = 1,
#     particle_material = 'PS',
#     particle_diameter = 3,
#     NAME = 'Si 120um, 3um PS (0.1 g/L)',
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0021/EXPERIMENT_0000/MEASUREMENT_0002/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor021b',
#     exposure_time = 0.05,
#     frame_gap = 1,
#     particle_material = 'PS',
#     particle_diameter = 3,
#     NAME = 'Si 120um, 3um PS (0.1 g/L)',
# )

# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0022/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor022a',
#     exposure_time = 0.05,
#     frame_gap = 1,
#     particle_material = 'PS',
#     particle_diameter = 4,
#     NAME = 'Si 120um, 4um PS (0.125 g/L)',
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0022/EXPERIMENT_0000/MEASUREMENT_0002/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor022b',
#     exposure_time = 0.05,
#     frame_gap = 5,
#     particle_material = 'PS',
#     particle_diameter = 4,
#     NAME = 'Si 120um, 4um PS (0.125 g/L)',
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0022/EXPERIMENT_0000/MEASUREMENT_0003/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor022c',
#     exposure_time = 0.05,
#     frame_gap = 1,
#     particle_material = 'PS',
#     particle_diameter = 4,
#     NAME = 'Si 120um, 4um PS (0.125 g/L)',
# )

# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0023/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor023a',
#     exposure_time = 0.05,
#     frame_gap = 1,
#     particle_material = 'PS',
#     particle_diameter = 4,
#     NAME = 'Si 120um, 4um PS (0.25 g/L)',
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0023/EXPERIMENT_0000/MEASUREMENT_0002/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor023a',
#     exposure_time = 0.05,
#     frame_gap = 5,
#     particle_material = 'PS',
#     particle_diameter = 4,
#     NAME = 'Si 120um, 4um PS (0.25 g/L)',
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0023/EXPERIMENT_0000/MEASUREMENT_0003/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor023c',
#     exposure_time = 0.05,
#     frame_gap = 1,
#     particle_material = 'PS',
#     particle_diameter = 4,
#     NAME = 'Si 120um, 4um PS (0.25 g/L), interface',
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0023/EXPERIMENT_0000/MEASUREMENT_0004/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor023d',
#     exposure_time = 0.05,
#     frame_gap = 5,
#     particle_material = 'PS',
#     particle_diameter = 4,
#     NAME = 'Si 120um, 4um PS (0.25 g/L), interface',
# )

# ## 25 is faxtor_manual.py (manual moving sample out the way for darks)

# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0026/EXPERIMENT_0000/MEASUREMENT_0003/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor026c',
#     exposure_time = 0.05,
#     frame_gap = 1,
#     particle_material = 'PS',
#     particle_diameter = 8,
#     NAME = 'Si 120um, 8um PS injected',
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0026/EXPERIMENT_0000/MEASUREMENT_0004/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor026d',
#     exposure_time = 0.05,
#     frame_gap = 1,
#     particle_material = 'PS',
#     particle_diameter = 8,
#     NAME = 'Si 120um, 8um PS injected, top',
# )

# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0028/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor028a',
#     exposure_time = 0.03,
#     frame_gap = 0.25,
#     particle_material = 'PS',
#     particle_diameter = 8,
#     NAME = 'Si 120um, 8um PS injected',
#     every_nth = 5,
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0028/EXPERIMENT_0000/MEASUREMENT_0005/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor028b',
#     exposure_time = 0.03,
#     frame_gap = 0.25,
#     particle_material = 'PS',
#     particle_diameter = 8,
#     NAME = 'Si 120um, 8um PS injected',
#     every_nth = 5,
# )

# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0029/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor029a',
#     exposure_time = 0.03,
#     frame_gap = 0.25,
#     particle_material = 'vermiculite',
#     NAME = 'Si 120um, vermiculite 2-10um injected',
#     every_nth = 5,
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0029/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor029a_part1',
#     exposure_time = 0.03,
#     frame_gap = 0.25,
#     particle_material = 'vermiculite',
#     NAME = 'Si 120um, vermiculite 2-10um injected',
#     fraction_of = 6,
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0029/EXPERIMENT_0000/MEASUREMENT_0004/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor029b',
#     exposure_time = 0.03,
#     frame_gap = 0.25,
#     particle_material = 'vermiculite',
#     NAME = 'Si 120um, vermiculite 2-10um injected',
#     # fraction_of = 5,
#     every_nth = 25,
# )

# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0030/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor030a',
#     exposure_time = 0.03,
#     frame_gap = 0.25,
#     particle_material = 'SiO2',
#     particle_diameter = 4.3,
#     NAME = 'Si 120um, Si 4.3um injected',
#     # every_nth = 5,
#     trim_start = 3800
# )
go(
    sample_folder = f'{RAW_DATA_PATH}/SAMPLE_0030/EXPERIMENT_0000/MEASUREMENT_0004/PCO_EDGE',
    flat_folder   = None,
    dark_folder   = None,
    name = 'faxtor030b',
    exposure_time = 0.03,
    frame_gap = 0.5,
    particle_material = 'SiO2',
    particle_diameter = 4.3,
    NAME = 'Si 120um, Si 4.3um injected',
    fraction_of=3
)

# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0031/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor031',
#     exposure_time = 0.03,
#     frame_gap = 0.5,
#     particle_material = 'vermiculite',
#     NAME = 'vermiculite 2-10um bulk',
#     every_nth = 5,
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0032/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor032',
#     exposure_time = 0.03,
#     frame_gap = 0.5,
#     particle_material = 'vermiculite',
#     NAME = 'vermiculite 2-10um bulk',
#     every_nth = 5,
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0033/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor033',
#     exposure_time = 0.03,
#     frame_gap = 0.5,
#     particle_material = 'vermiculite',
#     NAME = 'Si 85um, vermiculite 2-10um injected',
#     every_nth = 2,
# )

# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0034/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor034',
#     exposure_time = 0.03,
#     frame_gap = 0.25,
#     particle_material = 'SiO2',
#     particle_diameter = 7.75,
#     NAME = 'Si 85um, Si 7.75um injected',
#     trim_start = 2500,
#     every_nth = 1,
#     fraction_of = 2,
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0034/EXPERIMENT_0000/MEASUREMENT_0002/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor034flat',
#     exposure_time = 0.05,
#     frame_gap = 0,
#     NAME = 'flatfield',
# )

# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0035/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor035',
#     exposure_time = 0.03,
#     frame_gap = 0.5,
#     particle_material = 'SiO2',
#     particle_diameter = 4.3,
#     NAME = 'Si 85um, Si 4.3um injected',
#     # every_nth = 4,
#     trim_start=2700,
#     fraction_of=3,
# )

# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0036/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor036',
#     exposure_time = 0.03,
#     frame_gap = 0.5,
#     particle_material = 'PS',
#     particle_diameter = 8,
#     NAME = 'Si 120um, PS 8um injected',
#     fraction_of = 5,
# )

# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0037/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor037',
#     exposure_time = 0.03,
#     frame_gap = 0.5,
#     particle_material = 'PS',
#     particle_diameter = 8,
#     NAME = 'Si 85um, PS 8um injected',
#     fraction_of = 5,
# )

# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0038/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor038a',
#     exposure_time = 0.03,
#     frame_gap = 0.5,
#     particle_material = 'PS',
#     particle_diameter = 4,
#     NAME = 'Si 120um, PS 4um injected',
#     fraction_of = 5,
# )
# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0038/EXPERIMENT_0000/MEASUREMENT_0004/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor038b',
#     exposure_time = 0.03,
#     frame_gap = 0.5,
#     particle_material = 'PS',
#     particle_diameter = 4,
#     NAME = 'Si 120um, PS 4um injected, scintillator back',
#     fraction_of = 5,
# )

# go(
#     sample_folder = '/data2/acarter/faxtor/SAMPLE_0039/EXPERIMENT_0000/MEASUREMENT_0001/PCO_EDGE',
#     flat_folder   = None,
#     dark_folder   = None,
#     name = 'faxtor039',
#     exposure_time = 0.03,
#     frame_gap = 0.5,
#     particle_material = 'PS',
#     particle_diameter = 4,
#     NAME = 'Si 85um, PS 4um injected',
#     fraction_of = 5,
# )