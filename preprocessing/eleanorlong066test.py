import numpy as np
import common
import os, time
import tqdm

if __name__ == '__main__':
    t0 = time.time()
    a = np.load(f'particle_detection/data/particles_eleanorlong066.npy')
    t1 = time.time()
    print(f'took {t1-t0:.0f}s')

    a = None

    t0 = time.time()
    a = np.load(f'particle_detection/data/particles_eleanorlong066.npz')
    b=a['particles']
    t1 = time.time()
    print(f'took {t1-t0:.0f}s')


    """
    This should contain x,y,t for the duration of the high phi experiment. It is 
    saved in 155 separate numpy arrays (columns x,y,t) each spanning 1000 frames 
    of the experiment (this is an artefact of the way we save videos but I have 
    found helpful for dealing with larger data sets piecewise rather than all at 
    once). The slightly annoying thing is that to cut it down to 16 bit precision 
    the t coordinate ranges from 0-999 in each file, but it should be okay to add 
    the right multiple of 1000 based on the index in the file name.

    I don't know how much ram you have - but if you can't get it into a format 
    where all the coordinates are in a single file I can send you the hacky code 
    I use for counting in boxes one file at a time and then stitching it all 
    together at the end

    Sorry this is slightly hacky! Hope it's not too irritating to use. 21hrs of 
    high packing fraction is just a lot of coordinates haha
    """