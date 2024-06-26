import czifile
import common
import numpy as np
import json


def do_file(i, file, kind):
    with czifile.CziFile(file) as czi:
        image = czi.asarray()
        stack = image.squeeze()

        metadata = czi.metadata(raw=False)["ImageDocument"]["Metadata"]
        pixel_x = metadata["Scaling"]["Items"]["Distance"][0]["Value"]
        pixel_y = metadata["Scaling"]["Items"]["Distance"][0]["Value"]
        assert pixel_x == pixel_y
        pixel_x *= 1e6

        lensNA = metadata["Information"]["Instrument"]["Objectives"]["Objective"]["LensNA"]
        depth_of_field = 0.55 / lensNA**2 # https://www.physics.hmc.edu/~gerbode/wppriv/wp-content/uploads/2012/06/DS_ZISRAW-FileFormat_Rel_ZEN2011.pdf
        print(depth_of_field)

        time_step = metadata["Information"]["Image"]["Dimensions"]["T"]["Positions"]["Interval"]["Increment"]

        print(f'{file} is marine{kind}{i}, px={pixel_x:.4f}, t={time_step:.3f}, dof={depth_of_field:.3f}')
        np.savez(f'preprocessing/data/stack_marine{kind}{i}', stack=stack, pixel_size=pixel_x, time_step=time_step,
                 particle_diameter=np.nan, depth_of_field=depth_of_field)

    # with open('preprocessing/data/czi.json', 'w') as o:
    #     json.dump(metadata, o, sort_keys=False, indent=4)
    # break
        
        
for i, file in enumerate(common.get_directory_files('raw_data/marine', 'czi', 'DDTP_TL')):
    do_file(i, file, 'A')
        
for i, file in enumerate(common.get_directory_files('raw_data/marine', 'czi', 'PPR_TL')):
    do_file(i, file, 'B')