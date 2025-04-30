import common
import particle_detection.add_drift_periodic

# V_PROFILE = 'cos'
# V_PROFILE = 'cos2'
# V_PROFILE = 'step'
V_PROFILE = 'const'
# V_PROFILE = 'step2'

def go(file, v_profile, velocity_multiplier=1):
    particle_detection.add_drift_periodic.go_internal(
        infile = f'particle_linking/data/trajs_{file}.npz',
        outfile = f'particle_linking/data/trajs_{file}_drifted_{v_profile}_v{velocity_multiplier}.npz',
        v_profile = v_profile,
        velocity_multiplier = velocity_multiplier,
    )

if __name__ == '__main__':
    for file in common.files_from_argv('particle_linking/data', 'trajs_'):
        go(file, V_PROFILE)