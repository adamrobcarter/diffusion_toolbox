import common

if __name__ == '__main__':
    output = ''

    for file in common.files_from_argv('particle_linking/data', 'trajs_'):
        data = common.load(f'particle_linking/data/trajs_{file}.npz')
        time_step = data['time_step']

        drift, drift_std = common.find_drift(data)

        output += f'{file}\n'

        for dimension in range(common.particles_num_dimensions(data)):
            name = 'xyz'[dimension]
            s = f'drift {name} = {common.format_val_and_unc(drift[dimension]/time_step, drift_std[dimension]/time_step, latex=False)} um/s'
            output += s + '\n'
        
        output += '\n'
    
    print(output)