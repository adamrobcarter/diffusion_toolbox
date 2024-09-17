import os, sys
import tifffile
import numpy as np
import datetime, math
import scipy.fft
import tqdm
import matplotlib.animation, matplotlib.cm
import numba
import scipy.stats
import warnings, time
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import termplotlib
import psutil
import matplotlib.pyplot as plt

old_colors = plt.rcParams['axes.prop_cycle']
plt.style.use('dark_background')
plt.rcParams['axes.prop_cycle'] = old_colors
FIT_COLOR = 'white'

def get_directory_files(directory, extension, file_starts_with=''):
    filenames = []
    all_filenames = os.listdir(directory)
    all_filenames.sort()
    for filename in all_filenames:
        if filename.endswith(extension) and filename.startswith(file_starts_with):
            filenames.append(f'{directory}/{filename}')
    endstring = f' starting with {file_starts_with}' if file_starts_with else ''
    print(f'found {len(filenames)} .{extension}s {endstring}in {directory}')
    return filenames

def load_tif(filename):
    tif = tifffile.imread(filename)
    return tif

def intensity_correlation(data1, data2):
    correlation = np.multiply.outer(data1, data2, dtype='int16')
    # correlation[a][b][c][d] is data1[a][b] * data2[c][d]
    print(correlation.shape, correlation.dtype)
    x = np.arange(0, data1.shape[0], dtype='float16')
    y = np.arange(0, data1.shape[1], dtype='float16')

    # x1, y1 = np.meshgrid(x, y)
    # x2, y2 = np.copy(x1), np.copy(y1)
    x1 = x[:, np.newaxis, np.newaxis, np.newaxis]
    y1 = y[np.newaxis, :, np.newaxis, np.newaxis]
    x2 = x[np.newaxis, np.newaxis, :, np.newaxis]
    y2 = y[np.newaxis, np.newaxis, np.newaxis, :]

    r_sq = (x1 - x2)**2 + (y1 - y2)**2

    print(r_sq.shape)

def load(filename):
    if not filename.endswith('.npz'):
        filename += '.npz'
    modified = datetime.datetime.fromtimestamp(os.path.getmtime(f'{filename}'))
    diff = datetime.datetime.now() - modified
    print(f'loading {filename}, last modified {str(diff)[:-10]} ago')
    # print(f'loading {filename}, last modified ago')

    # try:
    data = np.load(f'{filename}', allow_pickle=True)
    # except FileNotFoundError:
    #     psiche = filename.split('/')[-1]
    #     print(psiche)
    #     data = np.load(f'/media/com-psiche/Sans titre/psiche_export1_npzFiles/{psiche}', allow_pickle=True)

    if filename.endswith('.npz'):
        for key in data.keys():
            if data[key].shape: # array
                print(f'  loaded {key.ljust(20)} dtype={str(data[key].dtype).ljust(12)} shape={data[key].shape}, size={arraysize(data[key])}')
            else: # single value
                print(f'  loaded {key.ljust(20)} dtype={str(data[key].dtype).ljust(12)} value={format_value_for_save_load(data[key])}')
    else:
        print(f'  loaded, dtype={data.dtype}, shape={data.shape}, size={arraysize(data)}')
        
    return data

def format_value_for_save_load(value):
    if hasattr(value, 'dtype'):
        if value.dtype.type is np.str_:
            return f'"{value}"'
    return value

def save_data(filename, quiet=False, **data):
    if not quiet:
        print(f'saving {filename}')
        for key in data.keys():
            if isinstance(data[key], np.ndarray):
                if data[key].shape: # array
                    if data[key].size * data[key].itemsize > 10e9:
                        warnings.warn(f'Saving array of size {arraysize(data[key])}')
                    print(f'  saving {key.ljust(20)} dtype={str(data[key].dtype).ljust(12)} shape={data[key].shape}, size={arraysize(data[key])}')
                else: # single value
                    print(f'  saving {key.ljust(20)} dtype={str(data[key].dtype).ljust(12)} value={format_value_for_save_load(data[key])}')
            else:
                if type(data[key]) in [list,tuple]:
                    nparray = np.array(data[key])
                    print(f'  saving {key.ljust(20)} dtype={str(nparray.dtype).ljust(12)} shape={nparray.shape}, size={arraysize(nparray)}')
                else:
                    print(f'  saving {key.ljust(20)} type={str(type(data[key])).ljust(12)} value={format_value_for_save_load(data[key])}')
    
    np.savez(filename, **data)

def arraysize(arr, mult=1):
    size = arr.size * arr.itemsize * mult
    if size < 1e3:
        # return str(size) + 'B'
        return f'{size}B'
    if size < 1e6:
        # return str(size/1e3) + 'kB'
        return f'{size/1e3:.0f}kB'
    if size < 1e9:
        # return str(size/1e6) + 'MB'
        return f'{size/1e6:.0f}MB'
    else:
        # return str(size/1e9) + 'GB'
        return f'{size/1e9:.0f}GB'

def fourier(t, x):
    N = x.shape[0]
    r_spacing = t[1] - t[0]
    X = scipy.fft.rfft(x)
    f = scipy.fft.rfftfreq(N, r_spacing)
    # S = scipy.fft.fftshift(S)
    # k = scipy.fft.fftshift(k)
    # return k[N//2:], S[N//2:]
    X = np.abs(X)

    return f, X#f[:N//2], X[:N//2]

def fourier_2D(x, spacing, axes):
    assert len(axes) == 2
    assert max(axes) < len(x.shape), f'you have asked for axis {max(axes)}, but your data only has {len(x.shape)} axes'
    N_x = x.shape[axes[0]]
    N_y = x.shape[axes[1]]
    X = scipy.fft.fftn(x, axes=axes, workers=16)
    f_x = scipy.fft.fftfreq(N_x, spacing)
    f_y = scipy.fft.fftfreq(N_y, spacing)
    f_xx, f_yy = np.meshgrid(f_y, f_x) # idk why these need to be the wrong way round
    # f_xx, f_yy = np.meshgrid(f_x, f_y)
    # S = scipy.fft.fftshift(S)
    # k = scipy.fft.fftshift(k)
    # return k[N//2:], S[N//2:]
    # X = np.abs(X)

    return f_xx, f_yy, X#f[:N//2], X[:N//2]

def r_squared(y_data, y_pred):
    # from https://stackoverflow.com/a/37899817/1103752
    residual_sum = np.sum((y_data - y_pred       )**2)
    total_sum    = np.sum((y_data - y_data.mean())**2)
    return 1 - (residual_sum / total_sum)

def add_watermark(fig, only_plot):
    command = '.'.join(sys.argv[0].split('/')[4:]).split('.py')[0] + ' ' + ' '.join(sys.argv[1:])
    time_str = time.strftime('%x %X')
    label = f'{time_str} {command}'

    if not only_plot:
        fig.tight_layout()
        fig.text(0.005, 0.008, label, fontsize=7)

def save_fig(fig, path, dpi=100, only_plot=False, hide_metadata=False):
    # only_plot gets rid of the axes labels and border leaving you with just the chart area (for pictures)

    if '--noplot' in sys.argv:
        return

    fig.tight_layout()
    args = {}

    if hide_metadata:
        args['bbox_inches'] = 'tight'
    else:
        add_watermark(fig, only_plot)

    if only_plot:
        for ax in fig.axes:
            ax.set_axis_off() # hide axes, ticks, etc
        args['bbox_inches'] = 'tight'
        args['pad_inches'] = 0

    path2 = path.replace('*', '')
    # if not path2.startswith('/'):
    #     path2 = f'~/Michot_0624/toolbox/{path}'
    print(f'saving {path2}    {len(fig.axes)} axes')
    fig.savefig(path, dpi=dpi, **args)

def add_scale_bar(ax, pixel_size, color='black'):
    image_width = ax.get_xlim()[1] - ax.get_xlim()[0]
    target_scale_bar_length = image_width * pixel_size / 10
    possible_scale_bar_lengths = (0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)

    takeClosest = lambda num,collection:min(collection,key=lambda x:abs(x-num))
    scale_bar_length = takeClosest(target_scale_bar_length, possible_scale_bar_lengths)

    asb = AnchoredSizeBar(ax.transData, scale_bar_length/pixel_size,
                          rf"${scale_bar_length}\mathrm{{\mu m}}$",
                          loc='lower left', pad=0.1, borderpad=0.5, sep=5,
                          frameon=False, color=color)
    ax.add_artist(asb)

def save_gif(func, frames, fig, file, fps=1, dpi=300):
    if fps > 10:
        warnings.warn(f'fps = {fps:.1f} which seems dangerous')
        fps = 5

    progress = tqdm.tqdm(total=len(frames)+1) # unsure why +1 but it seems to be needed

    def frame(timestep):
        func(timestep)
        progress.update()

    # add_watermark(fig, False)

    ani = matplotlib.animation.FuncAnimation(fig, frame, frames=frames, interval=500, repeat=False)
    # ani.save(file, dpi=dpi, writer=matplotlib.animation.PillowWriter(fps=fps))
    writer = None
    used_fps = fps
    if file.endswith('.gif'):
        writer = matplotlib.animation.PillowWriter(fps=fps)
        fps = None

    fig.tight_layout()
    ani.save(file, dpi=dpi, fps=fps, writer=writer)
    progress.close()
    print(f'saved {file} fps={used_fps}')

    
def N2_nointer_2D(t, D0, N_mean, Lx, Ly=None):
    if Ly == None:
        Ly = Lx
    return 2 * N_mean * (1 - famous_f(4*D0*t/Lx**2) * famous_f(4*D0*t/Ly**2)) # countoscope eq. 2, countoscope overleaf doc

def N2_nointer_2D_drift(t, D0, N0, L, drift_x):
    # assert len(drift) == 2
    
    # Grace's presentation slide 7
    # we compute this in a slightly different way to in the presentation
    # as the denominator of tau_0 is zero when drift is zero, so we get tau_0 = np.inf
    # then when this goes on the denominator again, numpy does 1/inf which gives nan

    tau = 4 * D0 * t / L**2
    tau_n_denom = lambda n : (drift_x*t/L + n)**2
    F_of_tau_n = lambda n : famous_f(tau / tau_n_denom(n)) * np.sqrt(tau_n_denom(n) / tau)
    F_of_tau   = famous_f(tau) / np.sqrt(tau)

    N2  = 2 * N0 * (1 - tau*F_of_tau * (0.5 * (F_of_tau_n(1) + F_of_tau_n(-1)) - F_of_tau_n(0)))

    return N2

def N2_nointer_3D(t, D0, N_mean, Lx, Ly, Lz):
    return 2 * N_mean * (1 - famous_f(4*D0*t/Lx**2) * famous_f(4*D0*t/Ly**2) * famous_f(4*D0*t/Lz**2))

def famous_f(tau):
    # print(tau)
    if np.all(np.isinf(tau)):
        return np.zeros_like(tau)
    elif np.any(np.isinf(tau)):
        raise Exception("I can't yet handle a mix of nan and non-nan")
    
    return np.sqrt(tau / np.pi) * ( np.exp(-1/tau) - 1) + scipy.special.erf(np.sqrt(1/tau)) # countoscope eq. 2

def structure_factor_2d_hard_spheres(k, phi, sigma):
    # sigma is disk diameter
    rho = 4 * phi / (np.pi * sigma**2)

    phi = np.pi/4 * rho * sigma**2
    J0 = lambda x: scipy.special.jv(0, x)
    J1 = lambda x: scipy.special.jv(1, x)

    prefactor = np.pi / ( 6 * ( 1 - phi)**3 * k**2 )
    line1 = -5/4 * (1 - phi)**2 * k**2 * sigma**2 * J0(k * sigma / 2)**2
    line23 = 4 * ( (phi - 20) * phi + 7) + 5/4 * (1 - phi)**2 * k**2 * sigma**2
    line23factor = J1(k * sigma / 2)**2
    line4 = 2 * (phi - 13) * (1 - phi) * k * sigma * J1(k * sigma / 2) * J0(k * sigma / 2)
    c = prefactor * (line1 + line23*line23factor + line4)
    # ^^^ Thorneywork et al 2018
    
    S = 1 / (1 - rho * c) # Hansen & McDonald (3.6.10)

    return S

def N1N2_nointer(t, D0, N0, L, drift):
     # Grace's presentation slide 10
    print(drift)
    assert len(L) == 2
    assert len(drift) == 2

    tau_i = lambda i : 4 * D0 * t / L[i-1]**2 # (i-1) so i in {1, 2}
    tau = lambda j, n, i : tau_i(i) / (drift[j-1]*t/L[i-1] + n)**2
    F = lambda t : famous_f(t) / np.sqrt(t) if not np.all(np.isinf(t)) else np.zeros_like(t)

    N1N2  = N0 * np.sqrt(tau_i(1)*tau_i(2))
    # print(N1N2)
    N1N2 *= F(tau(1,1,1)) - F(tau(1,2,1)) - 0.5*(F(tau(1,2,1)) - F(tau(1,-2,1)))
    # print(F(tau(1,-2,1)))
    # print(N1N2)
    N1N2 *= 0.5*(F(tau(2,1,2)) + F(tau(2,-1,2))) - F(tau(2,0,2))
    print(tau(2,0,2), F(tau(2,0,2)))
    # print(N1N2)

    print()
    return N1N2

def N1N2_square(t, D0, N0, L, drift):
    # Grace's presentation slide 7
    # we compute this in a slightly different way to in the presentation
    # as the denominator of tau_0 is zero when drift is zero, so we get tau_0 = np.inf
    # then when this goes on the denominator again, numpy does 1/inf which gives nan

    tau = 4 * D0 * t / L**2
    tau_n_denom = lambda n : (drift*t/L + n)**2
    F_of_tau_n = lambda n : famous_f(tau / tau_n_denom(n)) * np.sqrt(tau_n_denom(n) / tau)
    F_of_tau   = famous_f(tau) / np.sqrt(tau)
    
    N1N2  = N0 * tau * F_of_tau
    N1N2 *= F_of_tau_n(1) - F_of_tau_n(-1) - 0.5*(F_of_tau_n(2) - F_of_tau_n(-2))

    return N1N2


@numba.njit
def numba_binned_statistic(x, y, bins):
    # print(x.min(), x.max())
    # bin_edges = np.linspace(0, 100, bins+1) # +1 because you need eg. 3 points to define 2 bins
    # this will have broken if u used to supply num_bins
    bin_edges = bins
    num_bins = len(bins)-1
    bin_sums   = np.zeros((num_bins))
    bin_counts = np.zeros((num_bins))

    for value_i in range(x.shape[0]):
        for bin_i in range(num_bins):
            # print(x[value_i], '>', bin_edges[bin_i], x[value_i] > bin_edges[bin_i])
            if x[value_i] < bin_edges[bin_i+1]:
                # print('y[value_i]', y[value_i])
                bin_sums  [bin_i] += y[value_i]
                bin_counts[bin_i] += 1
                # print(x[value_i], bin_i)
                # break

    # print('suns', bin_sums)
    # print('counts', bin_counts)
    # print()

    return bin_sums / bin_counts, bin_edges, None

@numba.njit
def numba_sum_3d_axis01(array):
    # numba compatible version of array.sum(axis=(0, 1))
    sums = np.zeros((array.shape[2]))
    for i in range(array.shape[0]):
        sums[i] = array[:, :, i].sum()
    return sums

@numba.njit
def numba_nanmean_2d_axis0(array):
    nanmeans = np.zeros(array.shape[1])
    for i in range(array.shape[1]):
        nanmeans[i] = np.nanmean(array[:, i])
    return nanmeans

@numba.njit
def numba_mean_2d_axis0(array):
    means = np.zeros(array.shape[1])
    for i in range(array.shape[1]):
        means[i] = np.mean(array[:, i])
    return means

@numba.njit
def numba_nanstd_2d_axis0(array):
    nanstds = np.zeros(array.shape[1])
    for i in range(array.shape[1]):
        nanstds[i] = np.nanstd(array[:, i])
    return nanstds

@numba.njit
def numba_p_assert(condition, message):
    if not condition:
        print('Assertion failed: ', message)

def exponential_integers(min, max, num=50):
    assert min < max, f'I got min={min}, max={max}'
    return np.unique(np.round(10**np.linspace(np.log10(min), np.log10(max), num)).astype('int'))

def add_drift(particles, drift_x, drift_y):
    # drift should be per frame, particles should be rows of x,y,t

    width_before  = particles[:, 0].max()
    height_before = particles[:, 1].max()
    max_t = particles[:, 2].max()
    # print('drift over full time as fraction of window:', drift_x*max_t/width_before, drift_y*max_t/height_before)
    assert drift_x * max_t < width_before,  'drift over full time is greater than window width'
    assert drift_y * max_t < height_before, 'drift over full time is greater than window height'

    particles_drifted = particles.copy()
    for i in range(particles_drifted.shape[0]):
        t = particles_drifted[i, 2]
        particles_drifted[i, 0] += t * drift_x
        particles_drifted[i, 1] += t * drift_y

    # now the problem is the viewing window has shifted, at late times there will be
    # no particles in the bottom left (for +ve drift) of the image
    # so we gotta do a crop
    max_t = particles[:, 2].max()
    crop_x_low = max(0, drift_x * max_t)
    crop_y_low = max(0, drift_y * max_t)
    crop_x_high = min(width_before  + drift_x * max_t, width_before)
    crop_y_high = min(height_before + drift_x * max_t, height_before)

    removed_rows = (
          (particles_drifted[:, 0] < crop_x_low)
        | (particles_drifted[:, 1] < crop_y_low)
        | (particles_drifted[:, 0] > crop_x_high)
        | (particles_drifted[:, 1] > crop_y_high)
    )
    assert removed_rows.sum() < removed_rows.size

    particles_drifted = particles_drifted[~removed_rows, :]
    
    # now we've got a square box but without it's bottom left corner at (0, 0)
    particles_drifted[:, 0] -= particles_drifted[:, 0].min()
    particles_drifted[:, 1] -= particles_drifted[:, 1].min()

    print(f'cropped {1-particles_drifted[:, 0].max()/width_before}, {1-particles_drifted[:, 1].max()/height_before}')

    # after doing this cropping, we may have totally removed some particles
    # so we need to resample the particle numbers so they're continuous again
    # but only if we were provided with particle IDs
    if particles.shape[1] == 4:
        particle_ids = np.unique(particles_drifted[:, 3])
        id_map = {}
        for i in range(len(particle_ids)):
            new_id = i
            old_id = particle_ids[i]
            id_map[old_id] = new_id
        for i in range(particles_drifted.shape[0]):
            particles_drifted[i, 3] = id_map[particles_drifted[i, 3]]

    return particles_drifted

def find_drift(particles):
    # tactic:
    # for each particle, do a linear fit to where it exists
    # then weight the fits by length
    assert particles.shape[1] == 4, 'you should provide linked data to find_drift'
    num_particles = int(particles[:, 3].max())
    assert particles[:, 3].min() == 0, 'particles should be zero-based'

    drift_xs = []
    drift_ys = []
    fit_lens = []

    skip = 1
    if num_particles > 1000:
        skip = 5
    if num_particles > 5000:
        skip = 25
    used_particles = list(range(0, num_particles, skip))

    for particle in tqdm.tqdm(used_particles):
        indexes = particles[:, 3] == particle
        x = particles[indexes, 0]
        y = particles[indexes, 1]
        t = particles[indexes, 2]

        if t.size < 2:
            # can't regress with only one point!
            continue
        
        res_x = scipy.stats.linregress(t, x)
        res_y = scipy.stats.linregress(t, y)

        drift_xs.append(res_x.slope)
        drift_ys.append(res_y.slope)
        fit_lens.append(indexes.sum())

    drift_x = np.average(drift_xs, weights=fit_lens)
    drift_y = np.average(drift_ys, weights=fit_lens)
    assert not np.isnan(drift_x)
    assert not np.isnan(drift_y)

    return drift_x, drift_y

def remove_drift(particles):
    print('starting drift removal')
    drift_x, drift_y = find_drift(particles)
    print('drift', (drift_x, drift_y))

    particles_driftremoved = add_drift(particles, -drift_x, -drift_y)

    residual_drift = find_drift(particles_driftremoved)
    print('residual drift', residual_drift)

    return particles_driftremoved

import fnmatch

def files_from_argv(location, prefix):
    files = sys.argv[1:]

    if len(files) == 1 and files[0].endswith('*'):
        target = prefix + files[0]
        
        files = []
        all_filenames = os.listdir(location)
        all_filenames.sort()
        print(all_filenames)

        filenames = fnmatch.filter(all_filenames, target)
        
        for filename in filenames:
            filename_wo_ext = '.'.join(filename.split('.')[:-1])
            filename_wo_stem = filename_wo_ext.split(prefix)[1]
            files.append(filename_wo_stem)

    assert len(files), f'no files found. searching for {prefix}'
    return files

def format_val_and_unc(val, unc, sigfigs=2, latex=True):
    digits_after_decimal = -math.floor(math.log10(abs(val))) + sigfigs-1
    # print(f'.{digits_after_decimal}f digi')
    if digits_after_decimal < 0:
        digits_after_decimal = 0
    # assert digits_after_decimal
    if latex:
        return f'{val:.{digits_after_decimal}f} \pm {unc:.{digits_after_decimal}f}'
    else:
        return f'{val:.{digits_after_decimal}f}±{unc:.{digits_after_decimal}f}'

def nanfrac(arr):
    return np.isnan(arr).sum() / arr.size

def add_drift_intensity(stack, drift):
    # drift is px/frame
    assert drift != 0
    assert np.mod(drift, 1) == 0 # check integer-ness
    fraction_of_width_kept = 0.6
    new_width = int(stack.shape[1] * fraction_of_width_kept)
    useable_timesteps = (stack.shape[1] - new_width) // drift
    useable_timesteps = min(useable_timesteps, stack.shape[0])
    print(new_width, useable_timesteps)

    output = np.full((useable_timesteps, new_width, new_width), np.nan)

    new_x = np.arange(useable_timesteps) * drift
    # new_x[new_x.size//2:] = new_x[:new_x.size//2-1:-1]

    for t in range(useable_timesteps):
        # print(t, useable_timesteps)
        output[t, :, :]
        output[t, :, :] = stack[t, t*drift:t*drift+new_width, :new_width] # we could actually make the y axis keep it's original height, but for now, we crop to square
        # output[t, :, :] = stack[0, new_x[t]:new_x[t]+new_width, :new_width] # we could actually make the y axis keep it's original height, but for now, we crop to square

    return output

def term_hist(data):
    assert not np.any(np.isnan(data)), 'term_hist: nan found in data'
    counts, bin_edges = np.histogram(data, bins=20)
    term_fig = termplotlib.figure()
    term_fig.hist(counts, bin_edges, force_ascii=False, orientation="horizontal")
    term_fig.show()

    
def print_memory_use(s=''):
    pid = os.getpid()
    python_process = psutil.Process(pid)
    memoryUse = python_process.memory_info().rss / 2.0**30  # memory use in GB...I think
    print(f'memory use: {memoryUse:.1f}GB', s)

    
names = {
    'psiche0': '?',
    'psiche0007': '?',
    'psiche0008': 'solution_mere_10dil_3mu_COOH_NaOH',
    'psiche0009': 'test_flat_long_corr',
    'psiche0010': 'emptyglasscap',
    'psiche0011': 'H2O_glass_cap',
    'psiche0012': 'PS2um_COOH_Au_solution_mere',
    'psiche0013': '',
    'psiche0014': 'vermiculite_1_2um-dil10',
    'psiche0020': 'silice7p75um_sediment_x6p8_z15',
    'psiche0021': 'silice1p7um',
    'psiche0024': 'silice4um_sediment',
    'psiche0026': 'silice1p7um_texp500ms',
    'psiche0027': 'silice1p7um_texp200ms',
    'psiche0028': 'silice1p7um_bottom_texp1s',
    'psiche0029': 'silice0p96_texp500ms',
    'psiche0030': 'silice0p96_texp500ms',
    'psiche0031': 'silice0p96_texp1s',
    'psiche0035': 'PS2um_Au_texp200ms',
    'psiche0036': 'PS2um_Au_texp500ms',
    'psiche0036': 'PS2um_Au_texp1s',    
    'psiche0040': 'silice0p9um_10dil_texp200ms',
    'psiche0041': 'silice0p9um_10dil_texp500ms',
    'psiche0042': 'silice0p9um_10dil_texp1s',
    'psiche0045': 'muscovite_100_200um',
    'psiche0047': 'muscovite_100-200um_time15min-silice4um',
    'psiche0048': 'muscovite 100-200um silica 4um interface',
    'psiche0049': 'muscovite_100-200um_time15min-silice4um',

    # 'psiche0047': 'muscovite_100-200um_time15min-silice4um',
    # 'psiche0048': 'muscovite 100-200um silica 4um interface',
    # 'psiche0049': 'muscovite_100-200um_time15min-silice4um',
}

def name(file):
    if file.startswith('psiche'):
        if file.endswith('_small'):
            file = file[:-6]

        # try:
        #     with open(f'/media/com-psiche/Sans titre/psiche_export_1/{file}/NAME', 'r') as namefile:
        #         name = namefile.read()
        #     return f'{file} {name}'
        # except FileNotFoundError:
        #     pass

        try:
            with open(f'raw_data/psiche/{file}/NAME', 'r') as namefile:
                name = namefile.read()
            return f'{file} {name}'
        except FileNotFoundError:
            pass
    return names.get(file, file)

def colormap(value, min=0, max=1):
    return matplotlib.cm.afmhot(np.interp(value, (min, max), (0.25, 0.85)))

def colormap_cool(value, min=0, max=1):
    return matplotlib.cm.summer(np.interp(value, (min, max), (0, 0.8)))



class DisplayScript:
    figs = []

    def go(self):
        raise NotImplementedError
    
    def run(self, *args):
        self.go(*args)

        for figdata in self.figs:
            path = figdata['path'] + '/' + figdata['name']
            if type(figdata['file_or_files']) is list:
                path += '_' + '_'.join(figdata['file_or_files'])
            else:
                path += '_' + figdata['file_or_files']

            save_fig(figdata['fig']+'.png', path, dpi=figdata['dpi'])

    def fig(self, path, name, file_or_files, subplots=(1, 1), figsize=None, dpi=100):
        fig, ax = plt.subplots(*subplots, figsize=figsize)
        self.figs.append({
            fig: fig,
            path: path,
            name: name,
            dpi: dpi,
            file_or_files: file_or_files,
        })
        return fig, ax