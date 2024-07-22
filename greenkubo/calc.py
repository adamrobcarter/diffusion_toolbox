import common
import MSD.MSD
import numpy as np
import scipy.integrate
import tqdm

for file in common.files_from_argv('particle_linking/data/', 'trajs_'):
    data = common.load(f'particle_linking/data/trajs_{file}.npz')

    particles_ = data['particles']
    time_step = data['time_step']
    
    particles = MSD.MSD.reshape(particles_) # (num particles) x (num timesteps) x (x, y)
    num_timesteps = particles.shape[1]

    # first compute v
    v = (particles[:, 1:, :] - particles[:, :-1, :]) / time_step

    J_of_t = lambda t: v[:, t, :].sum(axis=0)

    all_t = np.arange(num_timesteps)

    Dcm = np.full(all_t.size, np.nan)

    # Falck et al 2004, eq 10

    for t in tqdm.tqdm(all_t):
        time_origin = 0

        # for time_origin in time_origins:

        integrand = lambda t_prime : np.dot(J_of_t(0), J_of_t(t_prime))
        ts = range(t)
        integrand_data = [integrand(t_prime) for t_prime in ts]
        integral = scipy.integrate.trapezoid(integrand_data, ts)

        Dcm[t] = integral

    common.save_data(f'greenkubo/data/Dc_{file}.npz', Dcm=Dcm)