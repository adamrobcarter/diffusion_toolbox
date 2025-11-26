import isf.calc_f
import common

if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data', 'particles_'):
        isf.calc_f.go(file, window=True)