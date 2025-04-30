import common
import numpy as np

for file in common.files_from_argv('particle_detection/data', 'particles_'):
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    particles = data['particles']

    particles_t0 = particles[particles[:, 2] == 0, :]
    assert particles_t0.size

    times = np.unique(particles[:, 2])

    particles_out = np.full_like(particles, np.nan)

    for t in times[1:]:
        particles_t = np.copy(particles_t0)
        particles_t[:, 2] = t
        n = particles_t0.shape[0]
        t = int(t)
        particles_out[n*t:n*(t+1), :] = particles_t

    finite_rows = np.isfinite(particles_out).all(axis=1)
    particles_out = particles_out[finite_rows, :]

    newdata = dict(data)
    newdata['particles'] = particles_out
    common.save_data(f'particle_detection/data/particles_{file}_frozen.npz', **newdata)