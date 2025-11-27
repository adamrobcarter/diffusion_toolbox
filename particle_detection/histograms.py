import common

if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data/', 'particles_'):
        data = common.load(f'particle_detection/data/particles_{file}.npz')
        particles = data['particles']
        
        print('x:')
        common.term_hist(particles[:, 0])
        
        print('y:')
        common.term_hist(particles[:, 1])

        print('time:')
        common.term_hist(particles[:, 2])