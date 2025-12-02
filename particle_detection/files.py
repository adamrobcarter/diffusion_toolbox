# this file is just for seeing which files turn up when you use a wildcard 
import common

if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data/', 'particles_'):
        print(file)