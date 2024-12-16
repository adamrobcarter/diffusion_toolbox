import common
import numpy as np

datas = []

SOURCES = [
    'timescaleint_nofit_cropped_var',
    'boxcounting_collective_var',
    'f_first_first',
    # 'f_long',
    'timescaleint_fixexponent_var'
]
for SOURCE in SOURCES:
    K_CROSSOVER = 0.6
    L_CROSSOVER = 3
    print('L crossover', L_CROSSOVER, L_CROSSOVER/3)

    files = common.files_from_argv('box_counting/data/', 'counted_')
    assert len(files) == 2

    data_long  = common.load(f'visualisation/data/Ds_from_{SOURCE}_{files[1]}.npz')
    data_short = common.load(f'visualisation/data/Ds_from_{SOURCE}_{files[0]}.npz')

    if SOURCE.startswith('f'):
        k_long  = data_long ['ks']
        k_short = data_short['ks']

        assert k_long.min()  < K_CROSSOVER
        assert k_long.max()  > K_CROSSOVER
        assert k_short.min() < K_CROSSOVER
        assert k_short.max() > K_CROSSOVER

        assert data_long['max_time_hours'] > data_short['max_time_hours']
        assert (particle_diameter := data_long['particle_diameter']) == data_short['particle_diameter']

        data_long_end_index    = np.argmax(k_long  >  K_CROSSOVER)
        data_short_start_index = np.argmax(k_short >= K_CROSSOVER)
        print(k_long)
        # print
        print(data_long_end_index, data_short_start_index)
        assert data_short_start_index > 0
        assert data_long_end_index    > 0
        assert data_short_start_index < k_short.size
        assert data_long_end_index    < k_long.size

        def concat(arr0, arr1):
            if len(arr0.shape) == 2:
                d0 = arr0[:, :data_long_end_index]
                d1 = arr1[:, data_short_start_index:]
                return np.concatenate((d0, d1), axis=1)
            else:
                d0 = arr0[:data_long_end_index]
                d1 = arr1[data_short_start_index:]
                return np.concatenate((d0, d1), axis=0)

        ks     = concat(k_long,              k_short)
        Ds     = concat(data_long['Ds'],     data_short['Ds'])
        D_uncs = concat(data_long['D_uncs'], data_short['D_uncs'])

        assert ks.shape == Ds.shape
        assert Ds.shape == D_uncs.shape

        dataout = dict(
            Ds                = Ds,
            D_uncs            = D_uncs,
            ks                = ks,
            particle_diameter = particle_diameter,
            max_time_hours    = data_short['max_time_hours']
        )
        
        common.save_data(f'visualisation/data/Ds_from_{SOURCE}_{files[1]}_mergedD.npz', **dataout)

    else:
        L_long  = data_long ['Ls']
        L_short = data_short['Ls']

        assert L_long.min()  < K_CROSSOVER
        assert L_long.max()  > K_CROSSOVER
        assert L_short.min() < K_CROSSOVER
        assert L_short.max() > K_CROSSOVER

        assert data_long['max_time_hours'] > data_short['max_time_hours']
        assert (particle_diameter := data_long['particle_diameter']) == data_short['particle_diameter']

        data_long_start_index = np.argmax(L_long  >  L_CROSSOVER)
        data_short_end_index  = np.argmax(L_short >= L_CROSSOVER)
        # print(k_long)
        # print
        # print(data_long_end_index, data_short_start_index)
        assert data_short_end_index  > 0
        assert data_long_start_index > 0
        assert data_short_end_index  < L_short.size
        assert data_long_start_index < L_long.size

        def concat(arr_long, arr_short):
            if len(arr_short.shape) == 2:
                d_short = arr_short[:, :data_short_end_index]
                d_long  = arr_long [:, data_long_start_index:]
                return np.concatenate((d_short, d_long), axis=1)
            else:
                d_short = arr_short[:data_short_end_index]
                d_long  = arr_long [data_long_start_index:]
                return np.concatenate((d_short, d_long), axis=0)

        Ls     = concat(L_long,              L_short)
        Ds     = concat(data_long['Ds'],     data_short['Ds'])
        D_uncs = concat(data_long['D_uncs'], data_short['D_uncs'])

        assert Ls.shape == Ds.shape

        dataout = dict(
            Ds                = Ds,
            D_uncs            = D_uncs,
            Ls                = Ls,
            particle_diameter = particle_diameter,
            max_time_hours    = data_short['max_time_hours']
        )

        common.save_data(f'visualisation/data/Ds_from_{SOURCE}_{files[1]}_mergedD.npz', **dataout)