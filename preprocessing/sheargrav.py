"""
updated sheargrav mesu->isis import for dec 2025
"""

import os
import numpy as np
import json
import common
import tqdm
import warnings

def read_binary_file(dir) -> np.ndarray:
    pos_file  = f"{dir}/colloids.bin"
    meta_file = f"{dir}/binary_metadata.json"

    with open(meta_file, "r") as f:
        metadata = json.load(f)

    n_rows   = metadata["n_rows"]
    row_size = metadata["row_size"]
    dtype = np.dtype(metadata["dtype"])

    data = np.fromfile(pos_file, dtype=dtype)
    if data.size != n_rows * row_size:
        warnings.warn(f"data size {data.size} does not match expected {n_rows*row_size} ({n_rows}*{row_size})")
    data = data.reshape((-1, row_size))
    return data

def go(remote_dir):
    if remote_dir.endswith('/'):
        remote_dir = remote_dir[:-1]

    name = remote_dir.split('/')[-1]
    local_dir = f'raw_data/preprocessing/sheargrav/{name}'
    print(f'local_dir: {local_dir}')
    os.makedirs(local_dir, exist_ok=True)

    common.rsync(f'{remote_dir}/colloids.bin', f'{local_dir}/colloids.bin')
    common.rsync(f'{remote_dir}/metadata.json', f'{local_dir}/metadata.json')
    common.rsync(f'{remote_dir}/binary_metadata.json', f'{local_dir}/binary_metadata.json')


    data = read_binary_file(local_dir)
    with open(f"{local_dir}/metadata.json", "r") as f:
        metadata = json.load(f)

    data = read_binary_file(local_dir)
    # data is rows of t, [x, y, z]*n, [qw, qx, qy, qz]*n
    n_particles = metadata["n_bodies"]
    n_frames = data.shape[0]

    particles = data

    # particles = np.full((n_particles*n_frames, 9), np.nan)
    # next_free_row = 0

    # for i in tqdm.trange(n_frames):
    #     frame = data[i, 0]
    #     Xs = data[i, 1:1 + n_particles*3].reshape((n_particles, 3))
    #     Qs = data[i, 1 + n_particles*3:1 + n_particles*7].reshape((n_particles, 4))
        
    #     for particle_id in range(n_particles):
    #         particles[next_free_row, 0:3] = Xs[particle_id, :]
    #         particles[next_free_row, 3] = frame
    #         particles[next_free_row, 4] = particle_id
    #         particles[next_free_row, 5:9] = Qs[particle_id, :]
    #         next_free_row += 1

    # for i in tqdm.trange(n_frames, desc='reshaping'):
    #     t = data[i, 0]
    #     for particle_id in range(n_particles):
    #         base_index = 1 + particle_id * 7
    #         x  = data[i, base_index + 0]
    #         y  = data[i, base_index + 1]
    #         z  = data[i, base_index + 2]
    #         qw = data[i, base_index + 3]
    #         qx = data[i, base_index + 4]
    #         qy = data[i, base_index + 5]
    #         qz = data[i, base_index + 6]

    #         row_index = i*n_particles + particle_id
    #         particles[row_index, 0] = x
    #         particles[row_index, 1] = y
    #         particles[row_index, 2] = z
    #         particles[row_index, 3] = t
    #         particles[row_index, 4] = particle_id
    #         particles[row_index, 5] = qw
    #         particles[row_index, 6] = qx
    #         particles[row_index, 7] = qy
    #         particles[row_index, 8] = qz

    newdata = metadata
    newdata['dimension'] = np.nan # signal to use particles_labels instead
    newdata['source_file'] = local_dir
    newdata['extra_source_file'] = remote_dir

    Lx, Ly = newdata['L'][0], newdata['L'][1]
    del newdata['L']

    print('newdata')
    print(newdata)
    
    common.save_particles(name,
        particles        = particles,
        particles_labels = ['x', 'y', 'z', 't', 'id', 'qw', 'qx', 'qy', 'qz'],
        window_size_x    = Lx,
        window_size_y    = Ly,
        **newdata,
    )

    common.save_trajs(name,
        particles        = particles,
        particles_labels = ['x', 'y', 'z', 't', 'id', 'qw', 'qx', 'qy', 'qz'],
        window_size_x    = Lx,
        window_size_y    = Ly,
        **newdata,
    )
    print(name)
    print()

go('/store/cartera/shear/shear0_T300_theta0_L100_phi0.5_multi_EMmid_nblobs42_dt0.005_tmax60')
go('/store/cartera/shear/shear0.080357_T300_theta10_L100_phi0.5_multi_EMmid_nblobs42_dt0.005_tmax60')
# go('/store/cartera/shear/shear0.080357_T300_theta0_L100_phi0.5_multi_EMmid_nblobs42_dt0.005_tmax60')
# go('/store/cartera/shear/shear0_T300_theta10_L100_phi0.5_multi_EMmid_nblobs42_dt0.005_tmax60')