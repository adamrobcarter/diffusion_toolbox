import common

if __name__ == '__main__':
    for file in common.files_from_argv('particle_linking/data', 'trajs_'):
        data = common.load(f'particle_linking/data/trajs_{file}.npz')
        time_step = data['time_step']

        drift, drift_std = common.find_drift(data['particles'], data.get('dimension', 2))

        for dimension in range(data.get('dimension', 2)):
            name = 'xyz'[dimension]
            print(f'drift {name} = {common.format_val_and_unc(drift[dimension]/time_step, drift_std[dimension]/time_step, latex=False)} um/s')