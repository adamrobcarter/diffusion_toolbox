import common
import stepsize.stepsize

for file in common.files_from_argv('particle_linking/data/', 'trajs_'):
    data = common.load(f'particle_linking/data/trajs_{file}.npz')

    particles = data['particles']

    res = stepsize.stepsize.calc(particles)