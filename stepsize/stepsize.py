import numpy as np
import tqdm
import common

def calc(particles):
    
    # version that doesn't need reshaping - probably slower but a lot less memory
    assert particles[:, 3].min() == 0

    # need the data are sorted by particle ID
    particles = particles[particles[:, 3].argsort()]
    
    num_particles = int(particles[:, 3].max()) + 1
    num_timesteps = int(particles[:, 2].max()) + 1
    print(f'{num_particles} particles, {num_timesteps} timesteps')

    progress = tqdm.tqdm(total=num_particles)

    start_index = 0
    current_id = 0
    skipped = 0

    all_steps_x = np.full(particles.shape[0] - num_particles, np.nan)
    all_steps_y = np.full(particles.shape[0] - num_particles, np.nan)
    all_steps_next_index = 0

    while start_index < particles.shape[0]-1:
        current_id = particles[start_index, 3]
        end_index = start_index + 1
        # print(start_index, end_index, particles.shape[0])
        # find the index of the last row that has this particle ID
        while end_index < particles.shape[0] and particles[end_index, 3] == current_id:
            end_index += 1

        data_this_particle = particles[start_index:end_index, :]
        assert np.all(data_this_particle[:, 3] == current_id)
        
        # now sort by time
        data_this_particle = data_this_particle[data_this_particle[:, 2].argsort()]

        if data_this_particle.shape[0] == 0:
            pass

        else:
            # assert data_this_particle[0, 2] == 0
            # num_timesteps_this_particle = int(data_this_particle[:, 2].max()) + 1
            num_timesteps_this_particle = int(data_this_particle[-1, 2]) + 1
            # print(num_timesteps_this_particle, data_this_particle.shape)
            # print(data_this_particle)
            # [print(data_this_particle[i, :]) for i in range(data_this_particle.shape[0])]
            # print(num_timesteps_this_particle, data_this_particle.shape[0])
            if num_timesteps_this_particle != data_this_particle.shape[0]:
                # this means that the timestep was non-contiguous
                skipped += 1
            
            else:

                x_steps = data_this_particle[1:, 0] - data_this_particle[:-1, 0]
                y_steps = data_this_particle[1:, 1] - data_this_particle[:-1, 1]

                # steps = np.sqrt(x_steps**2 + y_steps**2)
                # assert not np.any(np.isnan(steps))

                all_steps_x[all_steps_next_index:all_steps_next_index+x_steps.size] = x_steps
                all_steps_y[all_steps_next_index:all_steps_next_index+x_steps.size] = y_steps
                all_steps_next_index += x_steps.size
                
        progress.update()
        start_index = end_index
    progress.close()

    print(f'skipped {skipped/num_particles:.2f}')

    all_steps_x = all_steps_x[:all_steps_next_index]
    all_steps_y = all_steps_y[:all_steps_next_index]

    common.term_hist(all_steps_x)
    common.term_hist(all_steps_y)

    print('x', all_steps_x.mean(), all_steps_x.std())
    print('y', all_steps_y.mean(), all_steps_y.std())