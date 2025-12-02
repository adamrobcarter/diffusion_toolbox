import common
import MSD.MSD
import numpy as np

if __name__ == '__main__':
    for file in common.files_from_argv('particle_linking/data/', 'trajs_'):
        data = common.load(f'particle_linking/data/trajs_{file}.npz')
        particles = data['particles']
        msds = MSD.MSD.calc_individuals(particles)

        common.save_data(f'MSD/data/msd_individuals_{file}', msds=msds, time_step=data['time_step'])