import numpy as np
import common
import time

# eleanorlong is
# x min max 0.6982900415999999 293.91458112
# y min max 0.7108168319999999 367.6492512
# num_timesteps 71997.0
def go(infile, outfile, orig_width, out_width):
    print('loading')
    t0 = time.time()
    data = np.loadtxt(infile)
    t1 = time.time()
    print(f'loaded in {t1-t0:.0f}s. shape', data.shape, common.arraysize(data))

    data[:, 2] -= data[:, 2].min()
    assert data[:, 2].min() == 0

    num_timesteps_before = int(data[:, 2].max() + 1)
    print('num timesteps before crop', num_timesteps_before)

    not_too_long = data[:, 2] < 26477
    print(f'keeping {not_too_long.sum()/not_too_long.size:.2f} (too long)')
    data = data[not_too_long, :]

    print('raw histogram:')
    common.term_hist(data[:, 1])
    data[:, 0] = data[:, 0] % orig_width
    data[:, 1] = data[:, 1] % orig_width

    print('modded into window:')
    common.term_hist(data[:, 1])
    keep = ( data[:, 0] > 0 ) & ( data[:, 0] < out_width ) &  ( data[:, 1] > 0 ) & ( data[:, 1] < out_width )
    assert keep.sum() == keep.size

    if orig_width > out_width:
        num_timesteps = int(data[:, 2].max() + 1)
        common.save_data(f'particle_detection/data/particles_{outfile}_nocrop.npz', particles=data,
                    time_step=0.5, particle_diameter=2.79, pack_frac_given=0.34,
                    window_size_x=orig_width, window_size_y=orig_width,
                    num_timesteps=num_timesteps)

    # print('cropped into window:')
    # common.term_hist(data[keep, 1])
    # print(f'keeping {keep.sum()/keep.size:.2f} (inside crop)')
    # data = data[keep, :]

    num_timesteps = int(data[:, 2].max() + 1)

    common.save_data(f'particle_detection/data/particles_{outfile}.npz', particles=data,
                time_step=0.5, particle_diameter=2.79, pack_frac_given=0.34,
                window_size_x=out_width, window_size_y=out_width,
                num_timesteps=num_timesteps)
    
    
# 0.34
# go('raw_data/brennan/spec_softetakt_long_run_dtau_0.025_nsave_2.suspension_phi_0.34_L_320_modified.txt', 'brennan_hydro_034', 320, 320)
# go('raw_data/brennan/noHydro2D_Leim_run_dt_0.0125_nsave_40.suspension_phi_0.34_L_640_modified.txt', 'brennan_nohydro_034', 640, 320)
# 0.66
# go('/data2/acarter/Spectral_Sophie_Boxes/data/spec_softetakt_long_run_dtau_0.025_nsave_4.suspension_phi_0.66_L_288_modified.txt', 'brennan_hydro_066', 288, 288)
# go('/data2/acarter/Spectral_Sophie_Boxes/data/noHydro2D_Leim_run_dt_9.765625e5_nsave_256_long.suspension_phi_0.66_L_640_eq_modified.txt', 'brennan_nohydro_066', 640, 288)
go('/data2/acarter/sim/RigidMultiblobsWall/Lubrication/Lubrication_Examples/Monolayer/data/noHydro2D_Leim_run_dt_1.25e2_nsave_40_long_D0p04.suspension_phi_0.34_L_1280_modified.txt', 'sim_nohydro_034', 1280, 1280)


"""
#################################
#################################
##########  Hydro Data ##########
#################################
#################################

#################################
# this is the longest and best run that I did for phi = 0.34 *with* hydro
# The naming is sort of wierd, but I used a timestep of 0.25(s) and saved every other step to a file,
# so the effective timestep (the time betwen data frames in the file) is dt=0.5(s).  
# it's quated in the file name as a dimensionless timestep dtau=0.025. This is a stupid way to name files and I
# stopped doing it later and just quoted dt
#################################
Lx = 320
Ly = 320
infile = "./data/spec_softetakt_long_run_dtau_0.025_nsave_2.suspension_phi_0.34_L_320_modified.txt"
outfile = "./Count_Data_Cpp/New_Py_Test_phi_0.34" #
Nframes = 26478 # number of data frames (this might be less than there actually are)
a = 1.395 #radius of particles


#################################
# this is the longest and best run that I did for phi = 0.66 *with* hydro
# The naming is sort of wierd, but I used a timestep of 0.125(s) and saved every 4th step to a file,
# so the effective timestep (the time betwen data frames in the file) is dt=0.5(s).  
# it's quated in the file name as a dimensionless timestep dtau=0.025 *which is wrong in this case and I copy/pasted wrong* 
# it should've been dtau=0.0125. but dt=0.5(s) is right. 
#################################
Lx = 288
Ly = 288
infile = "./data/spec_softetakt_long_run_dtau_0.025_nsave_4.suspension_phi_0.66_L_288_modified.txt"
outfile = "./Count_Data_Cpp/Py_Test_phi_0.66"
Nframes = 8828 # number of data frames


#################################
#################################
########  No Hydro Data #########
#################################
#################################

# this is no-hydro data with an effective dt = (sim. timestep)*(save frequency) = 0.0125*40
# with packing fraction phi=0.34 and Lx=Ly=640
# to see how many frames are in the file just look at the third column of the last row
# $ tail -1 noHydro2D_Leim_run_dt_0.0125_nsave_40.suspension_phi_0.34_L_640_modified.txt
# 445.3523336427782 390.71188213483447 40000
infile = noHydro2D_Leim_run_dt_0.0125_nsave_40.suspension_phi_0.34_L_640_modified.txt


# this is no-hydro data with an effective dt = (sim. timestep)*(save frequency) = 9.765625e-5*32
# with packing fraction phi=0.34 and Lx=Ly=640
infile = noHydro2D_Leim_run_dt_9.765625e5_nsave_32.suspension_phi_0.34_L_640_modified.txt

# this is no-hydro data with an effective dt = (sim. timestep)*(save frequency) = 9.765625e-5*256
# with packing fraction phi=0.66 and Lx=Ly=640
noHydro2D_Leim_run_dt_9.765625e5_nsave_256_long.suspension_phi_0.66_L_640_eq_modified.txt
"""