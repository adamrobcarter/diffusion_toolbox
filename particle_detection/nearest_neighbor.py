import common
import numpy as np
import scipy
import tqdm

for file in common.files_from_argv('particle_detection/data/', 'particles'):
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    particles = data['particles']
    particle_diameter = data['particle_diameter']

    avg_dists = []

    for t in tqdm.trange(0, int(particles[:, 2].max()), 1000):
        particles_t = particles[particles[:, 2] == t, :]
        xy_t = particles_t[:, [0, 1]]

        tree = scipy.spatial.cKDTree(xy_t)

        all_distances, all_indexes = tree.query(xy_t, k=2)
        # find the nearest neighbors distances from the supplied points to the points in the tree
        # k is the number of nearest neighbors to find. k=2 b/c k=1 means finding yourself at zero distance
        print(all_distances)
        distances = all_distances[:, 1] # all_distances[:, 0] is the self distances, = 0
        avg_dists.append(distances.mean())

    print(f'avg nn dist = {common.format_val_and_unc(np.mean(avg_dists)/particle_diameter, np.std(avg_dists)/particle_diameter, latex=False)}σ')
