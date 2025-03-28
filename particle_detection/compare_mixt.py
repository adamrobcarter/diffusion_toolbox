import common
import numpy as np
import tqdm

file_mixt    = 'sim_nohydro_011_L320_test_mixt'
file_singlet = 'sim_nohydro_011_L320_test_singlet'

data_mixt    = common.load(f'particle_detection/data/particles_{file_mixt}.npz')
data_singlet = common.load(f'particle_detection/data/particles_{file_singlet}.npz')

particles_mixt    = data_mixt   ['particles']
particles_singlet = data_singlet['particles']

singlet_timestep = np.unique(particles_singlet[:, 2])[1]
print('singlet timestep', singlet_timestep)

mixt_ongrid = (particles_mixt[:, 2] % singlet_timestep) == 0

print('particles singlet mixt', particles_singlet.shape, particles_mixt.shape)
print('times', np.unique(particles_singlet[:, 2])[:5], np.unique(particles_mixt[:, 2])[:5])
print(mixt_ongrid.shape, particles_mixt.shape)
particles_mixt_ongrid = particles_mixt[mixt_ongrid, :]
assert 0.45 < particles_mixt_ongrid.size/particles_mixt.size < 0.55, particles_mixt_ongrid.size/particles_mixt.size
print(particles_singlet.shape, particles_mixt_ongrid.shape)

# assert np.all(particles_singlet == particles_mixt_ongrid)
for i in tqdm.trange(min(particles_singlet.shape[0], particles_mixt_ongrid.shape[0])):
    assert np.isclose(particles_singlet[i, :], particles_mixt_ongrid[i, :]).all()