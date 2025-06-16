import common
import numpy as np

for file in common.files_from_argv('particle_detection/data', 'particles_'):
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    particles = data['particles']
    assert particles.shape[1] == 3
    window_size_x = data['window_size_x']
    window_size_y = data['window_size_y']
    particle_diameter = data['particle_diameter']

    assert np.all(particles[0] >= 0)
    assert np.all(particles[0] <= window_size_x)
    assert np.all(particles[1] >= 0)
    assert np.all(particles[1] <= window_size_y)

    image_flip_x = np.copy(particles)
    image_flip_x[:, 0] = 2*window_size_x - image_flip_x[:, 0] + particle_diameter
    # if you wanted to do this on trajectories you should assign new ID numbers here
    # the ` + particle_diameter` prevents the images of particles overlapping with themselves at the border

    image_flip_y = np.copy(particles)
    image_flip_y[:, 1] = 2*window_size_y - image_flip_y[:, 1] + particle_diameter

    image_flip_xy = np.copy(particles)
    image_flip_xy[:, 0] = 2*window_size_x - image_flip_xy[:, 0] + particle_diameter
    image_flip_xy[:, 1] = 2*window_size_y - image_flip_xy[:, 1] + particle_diameter

    newparticles = np.full((particles.shape[0]*4, 3), np.nan)
    newparticles[0::4, :] = particles
    newparticles[1::4, :] = image_flip_x
    newparticles[2::4, :] = image_flip_y
    newparticles[3::4, :] = image_flip_xy

    newdata = common.copy_not_particles(data)
    newdata['particles'] = newparticles
    newdata['window_size_x'] = 2*window_size_x + particle_diameter
    newdata['window_size_y'] = 2*window_size_y + particle_diameter

    assert np.all(newparticles[:, 0] >= 0)
    assert np.all(newparticles[:, 0] <= newdata['window_size_x']), f'newparticles[:, 0].max() = {newparticles[:, 0].max()}, newdata[window_size_x] = {newdata['window_size_x']}'
    assert np.all(newparticles[:, 1] >= 0)
    assert np.all(newparticles[:, 1] <= newdata['window_size_y'])
    common.save_data(f'particle_detection/data/particles_{file}_mirrortile.npz', **newdata)