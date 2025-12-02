import common
import numpy as np

if __name__ == '__main__':
    data = np.loadtxt('raw_data/0.02_alldata.dat')
    print(data.shape)

    nmsd = data[:, 1:]

    nmsd = np.swapaxes(nmsd, 1, 0)

    common.save_data('box_counting/data/counted_sophiecountedhack.npz',
        N2_mean=nmsd,
        N2_std=np.zeros_like(nmsd),
        box_sizes=[32, 16, 8, 4, 3.2, 2.7, 2, 1],
        sep_sizes=np.zeros_like([32, 16, 8, 4, 3.2, 2.7, 2, 1]),
        # num_boxes=results.num_boxes,
        # N_mean=results.N_mean,
        # N_var=results.N_var,
        # N_var_mod=results.N_var_mod,
        # N_var_mod_std=results.N_var_mod_std,
        # N_mean_std=results.N_mean_std,
        time_step=0.5,
        pack_frac=0.01,
        particle_diameter=2.82,
        # particle_diameter_calced=particle_diameter_calced,
        # computation_time=time.time()-t0,
        # depth_of_field=depth_of_field,
        # old_overlap=use_old_overlap,
        # box_coords=results.box_coords,
        pack_frac_given=0.01,
        # window_size_x=window_size_x,
        # window_size_y=window_size_y,
        # pixel_size=data.get('pixel_size'),
        # **extra_to_save
    )