import numpy as np
import nd2
import common
import json

if __name__ == '__main__':
    for i, file in enumerate(common.get_directory_files('raw_data/marine3', 'nd2')):
        with nd2.ND2File(file) as f:
            # print(f.shape)
            full_stack = f.asarray()
            # print(f.attributes)
            # print(f.metadata)
            # print(f.unstructured_metadata())
            # with open('metadata.json', 'w') as jf:
            #     print(json.dump(f.unstructured_metadata(), jf, sort_keys=True, indent=4))
            time_step = f.unstructured_metadata()['ImageMetadataLV']['SLxExperiment']['uLoopPars']['dAvgPeriodDiff'] / 1000
            
            
            for channel in f.metadata.channels:
                
                dataset_num = int(file[-8:-4])

                assert channel.volume.axesCalibration[0] == channel.volume.axesCalibration[1]
                pixel_size = channel.volume.axesCalibration[0]

                if len(f.metadata.channels) > 1:
                    stack = full_stack[:, channel.channel.index, : :]
                    
                    if channel.channel.color.as_hex() == '#ff0000':
                        color = 'red'
                        filename = f'marine3_{i}r'
                    elif channel.channel.color.as_hex() == '#00ff00':
                        color = 'green'
                        filename = f'marine3_{i}g'
                    else:
                        raise Exception(f'color was {channel.channel.color.as_hex()}')
                    
                    total_color = 'red&green'

                else:
                    stack = full_stack
                    filename = f'marine3_{i}'

                    if 1024 <= dataset_num <= 1029:
                        color = 'green'
                    elif dataset_num == 1030:
                        color = 'red'
                    total_color = color

                if 1024 <= dataset_num <= 1029:
                    particle_type = 'DDTP'
                elif dataset_num == 1030:
                    particle_type = 'PPR'
                elif dataset_num >= 1031:
                    particle_type = 'DDTP&PPR'

                description = f'{particle_type} {total_color}'

                common.save_data(f'preprocessing/data/stack_{filename}',
                    stack=stack,
                    pixel_size=pixel_size,
                    time_step=time_step,
                    particle_diameter=np.nan,
                    NAME=description,
                    channel=color
                    #, depth_of_field=depth_of_field
                )
                # print(file, description, color)
            # break

            
                
            print()

            # break