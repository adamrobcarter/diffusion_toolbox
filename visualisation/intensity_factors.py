import numpy as np
import common
import sys

def find(file):
    print(file)
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    stack = data['stack']

    # get num particles from stack and known packing fraction
    # print(f'image {stack.shape[1]}px x {stack.shape[2]}px')
    # area = stack.shape[1] * stack.shape[2] * data['pixel_size']**2
    # phi = 0.34
    # # phi = np.pi/4 * rho * stack['particle_diameter']**2
    # rho = 4/np.pi * phi / data['particle_diameter']**2
    # num_particles = rho * area

    # get num particles from tracking
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    particles = data['particles'] # x, y, t
    num_timesteps = particles[:, 2].max() - particles[:, 2].min()
    num_particles = particles.shape[0]/num_timesteps
    print(f'{num_particles:.0f} particles')

    excess = stack[:, :, :] - stack.min(axis=0)
    print(stack.mean(), excess.mean())

    print(excess.sum(axis=(1,2)) / num_particles)
    print(excess.sum(axis=(1,2)).mean() / num_particles)
    print()
    return excess.sum(axis=(1,2)).mean() / num_particles

if __name__ == '__main__':
    for file in sys.argv[1:]:
        find(file)