import common

if __name__ == '__main__':
    for file in common.files_from_argv('particle_linking/data', 'particles_'):
        data = common.load(f'particle_linking/data/particles_{file}.npz')

        particles = data['particles']
        window_size_x = data['window_size_x']
        window_size_y = data['window_size_y']

        newdata = common.copy_not_particles(data)

        newdata['particles'] = common.periodic_unwrap(particles, data['dimension'], [window_size_x, window_size_y])

        common.save_data(f'particle_linking/data/particles_{file}_unwrap.npz', **newdata)