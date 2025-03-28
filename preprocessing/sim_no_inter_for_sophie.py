import numpy as np
import matplotlib.pyplot as plt

L =  # window size
phi = # packing fraction
sigma = # particle diameter
dt = # timestep between frames
D = # diffusion coefficiant
num_timesteps = 

num_particles = # function of L, phi, sigma

rng = np.random.default_rng()
# this is a `numpy.random.Generator`, the proper new way of doing random in numpy (it changed fairly recently)
# you can provide a seed to it which will mean it produces the same output every time you run the code (don't need to bother with that at the moment)
# following code should use this `rng` object to get random numbers

stepsize = # variance of the guassian distribution of the stepsize of a particle in one frame
           # function of D and dt

# now create the trajectories.
# generate random starting points for each particle
# an array of shape (num_particles)
# you can treat x and y separately
# suggest to use rng.uniform(): https://numpy.org/doc/2.2/reference/random/generated/numpy.random.Generator.uniform.html

# generate the gaussian steps and put the trajectories in an array
# an array of shape (num_particles x num_timesteps)
# you can still treat x and y separately
# the easiest way to do this would be to loop through each timestep
# but that is not the quickest way... we can talk about that
# suggest to use rng.normal(): https://numpy.org/doc/2.2/reference/random/generated/numpy.random.Generator.normal.html

# some of your coordinates might now be outside 0 < x (or y) < L
# wrap those back into the window so we have periodic boundary conditions


# now display the trajectories
fig, ax = plt.subplots(1, 1) # generate a figure with a single axes

# for each particle
    ax.plot(x, y) # x and y should be 1D arrays

fig.show()