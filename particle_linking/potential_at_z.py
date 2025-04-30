import numpy as np
import tqdm
import scipy
import matplotlib.pyplot as plt

def load_struct(file_name):
    with open(file_name, 'r') as f:
        lines = f.readlines()
        # s should be a float
        s = float(lines[0].split()[1])  
        Cfg = np.array([[float(j) for j in i.split()] for i in lines[1:]])
    return s, Cfg

struct_files = {
    42: 'particle_linking/structures/shell_N_42_Rg_0_8913_Rh_1.vertex',
    162: 'particle_linking/structures/shell_N_162_Rg_0_9497_Rh_1.vertex',
    642: 'particle_linking/structures/shell_N_642_Rg_0_9767_Rh_1.vertex',
    2562: 'particle_linking/structures/shell_N_2562_Rg_0_9888_Rh_1.vertex',
}

def wall_potential_blobs(r_vectors, a, debye_length, repulsion_strength):
    # reshape r_vecrors to be (N,3)
    r_vectors = r_vectors.reshape(-1, 3)
    #
    z = r_vectors[:,2]
    U = np.full_like(z, np.nan)
    lr_mask = z > a
    sr_mask = z <= a
    U[lr_mask] = repulsion_strength * np.exp(-(z[lr_mask]-a)/debye_length)
    U[sr_mask] = repulsion_strength + repulsion_strength / debye_length * (a - z[sr_mask])
    
    return U

fig_U, ax_U = plt.subplots(1, 1)
fig_P, ax_P = plt.subplots(1, 1)

def potential_at_z(z, hydrodynamic_radius, debye_length, repulsion_strength, nblobs, T, bare_radius):
    U = np.full_like(z, np.nan)
    
    atto = 1e-18
    k_B = scipy.constants.k / atto
    kT = k_B * T # in aJ

    # calculate magnitude of gravity force
    body_volume = 4/3*np.pi*bare_radius**3 # um^3
    rho_particles = 1510 # kg/m^3 (collective diffusion paper SI)
    rho_water = 970 # kg/m^3 (Eleanor in slack)
    delta_rho = rho_particles - rho_water # kg/m^3
    delta_rho *= (1e-6)**3 # kg/um^3
    effective_mass = delta_rho * body_volume # kg
    g = scipy.constants.g # ms^-2
    f_gravity_mag = effective_mass * g # N
    f_gravity_mag *= 1e12 # pN. approx 16kT

    struct_file = struct_files[int(nblobs)]
    s, Cfg = load_struct(struct_file)
    # Cfg is rows of x,y,z

    s *= hydrodynamic_radius
    Cfg *= hydrodynamic_radius
    blob_radius = 0.5*s

    for i in tqdm.trange(z.size):
        blob_coords = np.copy(Cfg)
        blob_coords[:, 2] += z[i]

        U_wall = wall_potential_blobs(blob_coords, blob_radius, debye_length, repulsion_strength)
        U_wall = U_wall.sum()
        U_grav = f_gravity_mag * z[i]
        U[i] = U_wall + U_grav

    return U

if __name__ == '__main__' and False:

    z = np.linspace(1, 2.5, 1000)

    hydrodynamic_diameter = 2.972 # um
    bare_diameter = 2.82 # um
    hydrodynamic_radius = hydrodynamic_diameter / 2
    bare_radius = bare_diameter / 2

    for nblobs in [42, 162, 642, 2562]:
        potential_at_z(z, U, hydrodynamic_radius, debye_length, repulsion_strength, f_gravity_mag, wall_potential_blobs, nblobs)

        ax_U.plot(z, U, label=f'nblobs={nblobs}')
        
        probability = np.exp(-U / kT)
        probability /= scipy.integrate.trapezoid(probability, z)
        ax_P.plot(z, probability, label=f'nblobs={nblobs}')

    ax_U.legend(fontsize=8)
    ax_U.set_xlabel('$z$')
    ax_U.set_ylabel('$U(z)$')
    ax_U.set_ylim(-1, 10)
    fig_U.savefig('U.png')

    ax_P.legend(fontsize=8)
    ax_P.set_xlabel('$z$')
    ax_P.set_ylabel('probability density')
    ax_U.set_ylim(-1, 10)
    fig_P.savefig('P.png')