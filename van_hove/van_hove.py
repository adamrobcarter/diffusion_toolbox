import numpy as np
import scipy.stats

def density_correlation(particles_t0, particles_t1, max_r, num_r_bins, crop_border):
    r_bin_edges = np.linspace(0, max_r, num_r_bins+1)
    r_bin_width = r_bin_edges[1] - r_bin_edges[0]
    densities = np.zeros((num_r_bins,))

    assert particles_t0.shape[0] > 0
    assert particles_t1.shape[0] > 0


    x = particles_t0[:, 0]
    y = particles_t0[:, 1]
    width  = max(particles_t0[:, 0].max(), particles_t1[:, 0].max())
    height = max(particles_t0[:, 1].max(), particles_t1[:, 1].max())

    out_of_crop = (x < crop_border) | (x > (width -crop_border)) | (y < crop_border) | (y > (height-crop_border))
    cropped_fraction = out_of_crop.sum() / x.size
    assert cropped_fraction < 1, 'all particles were cropped out'
    # used_locations = np.delete(locations, out_of_crop, axis=0)
    used_locations_t0 = particles_t0[~out_of_crop, :]
    num_used_particles_t0 = used_locations_t0.shape[0]

    # precalculate distance in one go
    used_locations_x_t0 = used_locations_t0[:, 0]
    used_locations_y_t0 = used_locations_t0[:, 1]
    locations_x_t1 = particles_t1[:, 0]
    locations_y_t1 = particles_t1[:, 1]
    dx = used_locations_x_t0[:, np.newaxis] - locations_x_t1[np.newaxis, :]
    dy = used_locations_y_t0[:, np.newaxis] - locations_y_t1[np.newaxis, :]
    r = np.sqrt(dx**2 + dy**2) #  r is num_particles x num_all_particles  where elements are the distance

    # find number of interparticle distances in donuts
    num_particles_in_donuts = np.histogram(r, r_bin_edges)[0]
    num_donuts = num_used_particles_t0
    assert num_donuts > 0
    avg_particles_per_donut = num_particles_in_donuts / num_donuts
    donut_areas = np.pi * (r_bin_edges[1:]**2 - r_bin_edges[:-1]**2) # approx = (left_edge + r_bin_width/2) * r_bin_width
    densities = avg_particles_per_donut / donut_areas

    r = (r_bin_edges[1:] + r_bin_edges[:-1])/2

    avg_density = num_used_particles_t0 / ( (width - 2*crop_border) * (height - 2*crop_border) )
    densities = densities / avg_density
    return r, densities, avg_density
