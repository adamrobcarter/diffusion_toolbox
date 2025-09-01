import common
import numpy as np
import tqdm

def calc_vacf(particles, num_dimensions):
    id_column = num_dimensions + 1
    time_column = num_dimensions
    
    num_particles  = int(particles[:, id_column  ].max()) + 1
    num_time_steps = int(particles[:, time_column].max()) + 1

    vacf_sum    = np.zeros((num_time_steps))
    vacf_counts = np.zeros((num_time_steps))

    # for particle in tqdm.trange(num_particles):

    particles = particles[np.lexsort((particles[:, time_column], particles[:, id_column]))] # sort by ID then time

    current_id = 0
    start_row = 0
    for row_id in tqdm.trange(particles.shape[0]):
        if particles[row_id, id_column] != current_id:
            # current_id = particles[row_id, id_column]

            # mask = particles[:, id_column] == particle
            # traj = particles[mask, :]
            traj = particles[start_row:row_id, :]
            # print(traj)
            # print()

            if not np.all(np.diff(traj[:, time_column]) == 1):
                print(f'time steps are not consecutive. np.diff(traj[:, time_column]).max() = {np.diff(traj[:, time_column]).max()}')
                # continue
                start_row = row_id
                current_id = particles[row_id, id_column]
                continue

            assert np.all(traj[:, id_column] == current_id)

            v = np.diff(traj[:, 0:num_dimensions], axis=0)
            # print('v shape:', v.shape)

            # vacf_this = 
            for time_origin in range(traj.shape[0]-1):
                # vacf_this = np.dot(v[time_origin:, :], v[time_origin, :])
                vacf_this = np.sum([v[time_origin:, dimension] * v[time_origin, dimension] for dimension in range(num_dimensions)], axis=0)
                # vacf_this = v[time_origin:, 0] * v[time_origin, 0]

                vacf_sum   [0:vacf_this.shape[0]] += vacf_this
                vacf_counts[0:vacf_this.shape[0]] += 1

            # vacf_this = v[:, 0] * v[0, 0]
            # print(vacf_this)
            # vacf_sum   [0:vacf_this.shape[0]] += vacf_this
            # vacf_counts[0:vacf_this.shape[0]] += 1

            start_row = row_id
            current_id = particles[row_id, id_column]

    vacf = vacf_sum / vacf_counts
    return vacf

for file in common.files_from_argv('particle_linking/data', 'trajs_'):
    data = common.load(f'particle_linking/data/trajs_{file}.npz')
    particles = data['particles']

    vacf = calc_vacf(particles, data.get('dimension', 2))
    assert (nanfrac := common.nanfrac(vacf)) < 0.1, f'vacf nanfrac = {nanfrac:.2g}'
    assert vacf[0] > 0, f'vacf[0] = {vacf[0]:.3g}'

    common.save_data(f'vacf/data/vacf_{file}.npz', vacf=vacf,
                     time_step=data.get('time_step', 1),
                     particle_diameter=data.get('particle_diameter'),
                     particle_mass=data.get('particle_mass'), eta=data.get('eta')
    )