import common

DIV = 16

if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data/', 'particles_'):
        data = common.load(f'particle_detection/data/particles_{file}_div8.npz')

        newdata = dict(data)

        num_timesteps = data['particles'][:, 2].max() + 1
        new_num_timesteps = int(num_timesteps / DIV)
        new_particles = data['particles'][data['particles'][:, 2] < new_num_timesteps]
        newdata['particles'] = new_particles
        newdata['max_time_hours'] = data['max_time_hours'] / DIV

        common.save_data(f'particle_detection/data/particles_{file}_div{8*DIV}.npz', **newdata)