import common
import numpy as np

if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data/', 'particles_'):
        data = common.load(f'particle_detection/data/particles_{file}.npz')
        newdata = dict(data)

        pack_frac_calced = common.calc_pack_frac(data['particles'], data['particle_diameter'], data['window_size_x'], data['window_size_y'])
        assert np.isclose(pack_frac_calced, data['pack_frac_given'], rtol=0.1), f'pack frac calced {pack_frac_calced}, given' + str(data['pack_frac_given'])

        newdata['pack_frac'] = pack_frac_calced
        common.save_data(f'particle_detection/data/particles_{file}.npz', **newdata)
