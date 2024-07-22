import numpy as np
import common

# eleanorlong is
# x min max 0.6982900415999999 293.91458112
# y min max 0.7108168319999999 367.6492512
# num_timesteps 71997.0

data = np.loadtxt('raw_data/brennan/spec_softetakt_long_run_dtau_0.025_nsave_2.suspension_phi_0.34_L_320_modified.txt')

data[:, 2] -= data[:, 2].min()
assert data[:, 2].min() == 0

not_too_long = data[:, 2] < 72000
print(f'keeping {not_too_long.sum()/not_too_long.size:.2f}')
data = data[not_too_long, :]

num_timesteps = int(data[:, 2].max() + 1)

common.save_data(f'particle_detection/data/particles_brennan_hydro_034.npz', particles=data,
            time_step=0.5, particle_diameter=2.79, pack_frac_given=0.34,
            num_timesteps=num_timesteps)



data = np.loadtxt('raw_data/brennan/noHydro2D_Leim_run_dt_0.0125_nsave_40.suspension_phi_0.34_L_640_modified.txt')

data[:, 2] -= data[:, 2].min()
assert data[:, 2].min() == 0

inside_crop = (data[:, 0] < 320) & (data[:, 1] < 320)
print(f'keeping {inside_crop.sum()/inside_crop.size:.2f}')
data = data[inside_crop, :]

not_too_long = data[:, 2] < 72000
print(f'keeping {not_too_long.sum()/not_too_long.size:.2f}')
data = data[not_too_long, :]

num_timesteps = int(data[:, 2].max() + 1)

common.save_data(f'particle_detection/data/particles_brennan_nohydro_034.npz', particles=data,
            time_step=0.5, particle_diameter=2.79, pack_frac_given=0.34,
            num_timesteps=num_timesteps)