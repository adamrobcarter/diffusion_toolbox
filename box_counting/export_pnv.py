import common
import matplotlib.pyplot as plt
import matplotlib.cm
import numpy as np
import csv

def go(file):
    data = common.load(f'box_counting/data/pnv_{file}.npz')
    

    N1N2 = data['N1N2']
    dt = data['time_step']
    
    # for box_size_index in range(len(data['box_sizes_x'])):
    L1 = data['box_sizes_x']
    L2 = data['box_sizes_y']
    N_mean = data['N_mean']
    N_var = data['N_var']

        
    csv_filename = f'box_counting/data/pnv_{file}.csv'
    with open(csv_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        writer.writerow(['L1'] + list(L1))
        writer.writerow(['L2'] + list(L2))
        writer.writerow(['N_mean'] + list(N_mean))
        writer.writerow(['N_var'] + list(N_var))


        for frame in range(N1N2.shape[1]):
            writer.writerow([frame * dt] + list(N1N2[:, frame]))
        
    print(f'wrote {csv_filename}')

if __name__ == '__main__':
    for file in common.files_from_argv('box_counting/data/', 'pnv_'):
        go(file)