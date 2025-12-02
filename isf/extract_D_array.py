import common
import numpy as np

def save_Ds(t, Ds, D_uncs, k, file, metadata):
    if t % 1 == 0:
        t = int(t) # save "t64" not "t64.0"

    common.save_data(f'visualisation/data/Ds_from_f_t{t}_array{file}',
        Ds=Ds, D_uncs=D_uncs, ks=k,
        **metadata
    )

if __name__ == '__main__':
    for file in common.files_from_argv('isf/data/', 'F_array'):
        data = common.load(f'isf/data/F_array{file}.npz')
        F = data['F_array']
        k = data['k']
        t = data['t']
        
        metadata = dict(
            particle_diameter=data['particle_diameter'], max_time_origins=data['max_time_origins'],
            pixel_size=data.get('pixel_size'), max_time_hours=data.get('max_time_hours'),
            channel=data.get('channel'), NAME=data.get('NAME'), computation_time=data['computation_time']
        )

        F_avg = F.mean(axis=0)
        F_std = F.std (axis=0)

        f = F_avg / F_avg[0, :]
        f_unc = np.sqrt((F_std / F_avg[0, :])**2 + (F_avg * F_std[0, :] / F_avg[0, :])**2)

        ts_nonzero = range(1, F_avg.shape[0])
        assert len(ts_nonzero)
        for t_index in ts_nonzero:
            D = - np.log(f[t_index, :]) / (k**2 * t[t_index])
            D_unc = f_unc[t_index, :] / (f[t_index] * k**2 * t[t_index])

            save_Ds(t[t_index], D, D_unc, k, file, metadata)