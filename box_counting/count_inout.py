import countoscope
import numpy as np
import time
import sys
import common

from box_counting.count import calc_and_save

def inout(particles_in, probability, lifetime):
    particles = np.copy(particles_in)
    max_t = int( particles[:, 2].max() )
    lifetime_abs = int(lifetime * max_t)
    print('min t', particles[:, 2].min())
    # assert particles.shape[1] == 4
    # for particle in range(int(particles[:, 3].max())):
    #     particles[:, 3] == particle
    #     row_indexes = np.where(particles[:, 3] == particle)
    #     for i in range(0, len(row_indexes)):
    #         if np.random.rand() < probability:
    #             print(f'deleting {i}', row_indexes[i:])
    #             particles = np.delete(particles, row_indexes[i:], axis=0)
    #             break

    #     row_indexes = np.where(particles[:, 3] == particle)
    #     for i in range(len(row_indexes), 0, -1):
    #         if np.random.rand() < probability:
    #             print(f'deletin2 {i}', row_indexes[:i])
    #             particles = np.delete(particles, row_indexes[:i], axis=0)
    #             break
    # return particles
    for particle in range(int(particles[:, 3].max())):
        if np.random.rand() < probability:
            row_indexes = np.where(particles[:, 3] == particle)[0]
            if row_indexes.shape[0] == 0:
                print(f'skipping {particle}')
                continue

            if particles[row_indexes[0], 2] != particles[:, 2].min():
                print(f'skippin2 {particle}')
                particles = np.delete(particles, row_indexes, axis=0)
                continue

            if particles[row_indexes[-1], 2] != max_t:
                print(f'skippin3 {particle}')
                particles = np.delete(particles, row_indexes, axis=0)
                continue

            in_index = np.random.randint(-lifetime_abs, max_t)
            out_index = in_index + lifetime_abs

            # in_index  = np.random.randint(0, int(row_indexes.shape[0]*(1-lifetime)))
            # in_index  = np.random.randint(0, row_indexes.shape[0])
            # out_index = in_index + lifetime_abs
            if out_index > max_t:
                out_index = max_t
            if in_index < 0:
                in_index = 0
            print(max_t, row_indexes.size, 'taking', in_index, 'to', out_index)
            deleted_indexes = np.concatenate([row_indexes[:in_index], row_indexes[out_index:]])
            particles = np.delete(particles, deleted_indexes, axis=0)
    return particles

for file in common.files_from_argv('box_counting/data', 'counted_'):
    # box_sizes_px = np.array([1,  2,  4,  8,  16,  32,  64])
    # sep_sizes_px = np.array([20, 20, 20, 3,  10, -10, -20])
    box_sizes_px = np.array([8])
    sep_sizes_px = np.array([-3])
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    data = dict(data)
    particles = data['particles']
    l = particles.shape[0]

    # data['particles'] = inout(particles, 1, 0.9)
    # print('rows', data['particles'].shape[0])
    # filename = f'box_counting/data/counted_{file}_inout0.9.npz'
    # calc_and_save(box_sizes_px, sep_sizes_px, data, filename, extra_to_save={'data_fraction':data['particles'].shape[0]/l})

    # data['particles'] = inout(particles, 1, 0.5)
    # print('rows', data['particles'].shape[0])
    # filename = f'box_counting/data/counted_{file}_inout0.5.npz'
    # calc_and_save(box_sizes_px, sep_sizes_px, data, filename, extra_to_save={'data_fraction':data['particles'].shape[0]/l})

    data['particles'] = inout(particles, 1, 0.1)
    print('rows', data['particles'].shape[0])
    filename = f'box_counting/data/counted_{file}_inout0.1.npz'
    calc_and_save(box_sizes_px, sep_sizes_px, data, filename, extra_to_save={'data_fraction':data['particles'].shape[0]/l})

    data['particles'] = inout(particles, 1, 0.01)
    print('rows', data['particles'].shape[0])
    filename = f'box_counting/data/counted_{file}_inout0.01.npz'
    calc_and_save(box_sizes_px, sep_sizes_px, data, filename, extra_to_save={'data_fraction':data['particles'].shape[0]/l})

    data['particles'] = inout(particles, 1, 0.005)
    print('rows', data['particles'].shape[0])
    filename = f'box_counting/data/counted_{file}_inout0.005.npz'
    calc_and_save(box_sizes_px, sep_sizes_px, data, filename, extra_to_save={'data_fraction':data['particles'].shape[0]/l})

    data['particles'] = inout(particles, 1, 1)
    print('rows', data['particles'].shape[0])
    filename = f'box_counting/data/counted_{file}_inout1.npz'
    calc_and_save(box_sizes_px, sep_sizes_px, data, filename, extra_to_save={'data_fraction':data['particles'].shape[0]/l})

# particles = np.array([
#     [0, 0, 1, 1],
#     [0, 0, 2, 1],
#     [0, 0, 3, 1],
#     [0, 0, 4, 1],
#     [0, 0, 5, 1],
#     [0, 0, 6, 1],
#     [0, 0, 7, 1],
#     [0, 0, 8, 1],
#     [0, 0, 1, 0],
#     [0, 0, 2, 0],
#     [0, 0, 3, 0],
#     [0, 0, 4, 0],
#     [0, 0, 5, 0],
#     [0, 0, 6, 0],
#     [0, 0, 7, 0],
#     [0, 0, 8, 0],
# ])
# print(inout(particles, 0.5))
# print(inout(particles, 0.5))
# print(inout(particles, 0.5))
# print(inout(particles, 0.5))
# print(inout(particles, 0.5))